// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file multimat.cpp
 *
 * \brief Implementation of the MultiMat class
 */

#include "axom/multimat/multimat.hpp"

#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "axom/slic/interface/slic.hpp"
#include <cassert>

#include "multimat.hpp"


using namespace std;
using namespace axom::multimat;


MultiMat::MultiMat(DataLayout d, SparsityLayout s) :
  m_ncells(0),
  m_nmats(0),
  m_dataLayout(d),
  m_sparsityLayout(s),
  m_cellMatRel(nullptr),
  m_cellMatRelDyn(nullptr),
  m_cellMatNZSet(nullptr),
  m_cellMatProdSet(nullptr),
  m_dynamic_mode(false)
{}

MultiMat::~MultiMat()
{
  for (auto mapPtr : m_mapVec)
  {
    delete mapPtr;
  }
  delete m_cellMatProdSet;
  delete m_cellMatNZSet;
  delete m_cellMatRel;
  delete m_cellMatRelDyn;
}

template<typename T>
MultiMat::MapBaseType* MultiMat::helper_copyField(const MultiMat& mm, int map_i)
{
  MapBaseType* other_map_ptr = mm.m_mapVec[map_i];
  if (mm.getFieldMapping(map_i) == FieldMapping::PER_CELL_MAT)
  {
    BivariateSetType* biSetPtr = get_mapped_biSet();
    Field2D<T>* typed_ptr = dynamic_cast<Field2D<T>*>(other_map_ptr);
    Field2D<T>* new_ptr = new Field2D<T>(biSetPtr, T(), typed_ptr->stride());
    new_ptr->copy(typed_ptr->getMap()->data().data());
    return new_ptr;
  }
  else
  {
    SetType* setPtr = get_mapped_set(getFieldMapping(map_i));
    Field1D<T>* typed_ptr = dynamic_cast<Field1D<T>*>(other_map_ptr);
    Field1D<T>* new_ptr = new Field1D<T>(setPtr, T(), typed_ptr->stride());
    new_ptr->copy(*typed_ptr);
    return new_ptr;
  }
}


MultiMat::MultiMat(const MultiMat& other) :
  m_ncells(other.m_ncells),
  m_nmats(other.m_nmats),
  m_dataLayout( other.m_dataLayout),
  m_sparsityLayout( other.m_sparsityLayout),
  m_cellSet(0,other.m_ncells),
  m_matSet(0, other.m_nmats),
  m_cellMatRel_beginsVec(other.m_cellMatRel_beginsVec),
  m_cellMatRel_indicesVec(other.m_cellMatRel_indicesVec),
  m_cellMatRel(nullptr),
  m_cellMatRelDyn(nullptr),
  m_cellMatNZSet(nullptr),
  m_cellMatProdSet(nullptr),
  m_arrNameVec(other.m_arrNameVec),
  m_fieldMappingVec(other.m_fieldMappingVec),
  m_dataTypeVec(other.m_dataTypeVec),
  m_dynamic_mode(false)
{
  RangeSetType& set1 = (isCellDom() ? m_cellSet : m_matSet);
  RangeSetType& set2 = (isCellDom() ? m_matSet : m_cellSet);
  m_cellMatRel = new StaticVariableRelationType(&set1, &set2);
  m_cellMatRel->bindBeginOffsets(set1.size(), &m_cellMatRel_beginsVec);
  m_cellMatRel->bindIndices(m_cellMatRel_indicesVec.size(),
                            &m_cellMatRel_indicesVec);
  m_cellMatNZSet = new RelationSetType(m_cellMatRel);
  m_cellMatProdSet = new ProductSetType(&set1, &set2);

  for (unsigned int map_i = 0 ; map_i < other.m_mapVec.size() ; ++map_i)
  {
    MapBaseType* new_map_ptr = nullptr;
    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble)
    {
      new_map_ptr = helper_copyField<double>(other, map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat)
    {
      new_map_ptr = helper_copyField<float>(other, map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt)
    {
      new_map_ptr = helper_copyField<int>(other, map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar)
    {
      new_map_ptr = helper_copyField<unsigned char>(other, map_i);
    }
    else
      SLIC_ASSERT_MSG(false,
                      "\t*MultiMat copy constructor : Unsupported Datatype");

    SLIC_ASSERT(new_map_ptr != nullptr);
    m_mapVec.push_back(new_map_ptr);
  }
}

MultiMat& MultiMat::operator=(const MultiMat& other)
{
  if (this == &other)
    return *this;

  throw axom::slam::NotImplementedException();

  return *this;
}


void MultiMat::setNumberOfMaterials(int n)
{
  SLIC_ASSERT(n > 0);
  m_nmats = n;

  m_matSet = RangeSetType(0, m_nmats);
  SLIC_ASSERT(m_matSet.isValid());
}

void MultiMat::setNumberOfCells(int c)
{
  SLIC_ASSERT(c > 0);
  m_ncells = c;

  m_cellSet = RangeSetType(0, m_ncells);
  SLIC_ASSERT(m_cellSet.isValid());
}

void MultiMat::setCellMatRel(vector<bool>& vecarr)
{
  //Setup the SLAM cell to mat relation
  //This step is necessary if the volfrac field is sparse

  SLIC_ASSERT(vecarr.size() == m_ncells * m_nmats); //Check it's dense
  SLIC_ASSERT(m_cellMatRel == nullptr); //cellmatRel has not been set
  SLIC_ASSERT(m_cellMatRelDyn == nullptr);

  RangeSetType& set1 = (isCellDom() ? m_cellSet : m_matSet);
  RangeSetType& set2 = (isCellDom() ? m_matSet : m_cellSet);

  //count the non-zeros
  int nz_count = 0;
  for (bool b : vecarr)
    nz_count += b;

  //Set-up the cell/mat relation
  m_cellMatRel_beginsVec.resize(set1.size() + 1, -1);
  m_cellMatRel_indicesVec.resize(nz_count);

  SetPosType curIdx = SetPosType();
  for (SetPosType i = 0 ; i < set1.size() ; ++i)
  {
    m_cellMatRel_beginsVec[i] = curIdx;
    for (SetPosType j = 0 ; j < set2.size() ; ++j)
    {
      if (vecarr[i*set2.size() + j])
      {
        m_cellMatRel_indicesVec[curIdx] = j;
        ++curIdx;
      }
    }
  }
  m_cellMatRel_beginsVec[set1.size()] = curIdx;

  m_cellMatRel = new StaticVariableRelationType(&set1, &set2);
  m_cellMatRel->bindBeginOffsets(set1.size(), &m_cellMatRel_beginsVec);
  m_cellMatRel->bindIndices(m_cellMatRel_indicesVec.size(),
                            &m_cellMatRel_indicesVec);

  SLIC_ASSERT(m_cellMatRel->isValid());

  //Set-up both dense and sparse BivariateSets.
  m_cellMatNZSet = new RelationSetType(m_cellMatRel);
  m_cellMatProdSet = new ProductSetType(&set1, &set2);

  //Create a field for VolFrac as the 0th field
  m_mapVec.push_back(nullptr);
  m_arrNameVec.push_back("Volfrac");
  m_fieldMappingVec.push_back(FieldMapping::PER_CELL_MAT);
  m_dataTypeVec.push_back(DataTypeSupported::TypeDouble);
  SLIC_ASSERT(m_mapVec.size() == 1);
  SLIC_ASSERT(m_arrNameVec.size() == 1);
  SLIC_ASSERT(m_fieldMappingVec.size() == 1);
  SLIC_ASSERT(m_dataTypeVec.size() == 1);
}


int MultiMat::setVolfracField(double* arr)
{
  //m_mapVec[0] should already be a volfrac map. This functions add a new map,
  //with the input arr, then swap the new map with the 0th map,
  //and delete the new map.

  //Volfrac map is a CellxMat mapping, named "Volfrac", and is stride 1.
  int arr_i = addFieldArray_impl<double>("Volfrac",
                                         FieldMapping::PER_CELL_MAT, arr, 1);

  //move the data to the first one (index 0) in the list
  std::iter_swap(m_mapVec.begin(), m_mapVec.begin() + arr_i);
  std::iter_swap(m_dataTypeVec.begin(), m_dataTypeVec.begin() + arr_i);

  //remove the new entry...
  int nfield = m_mapVec.size() - 1;
  m_mapVec.resize(nfield);
  m_fieldMappingVec.resize(nfield);
  m_arrNameVec.resize(nfield);
  m_dataTypeVec.resize(nfield);

  return 0;
}


MultiMat::Field2D<double>& MultiMat::getVolfracField()
{
  return *dynamic_cast<Field2D<double>*>(m_mapVec[0]);
}


int MultiMat::getFieldIdx(const std::string& field_name) const
{
  for (unsigned int i = 0 ; i < m_arrNameVec.size() ; i++)
  {
    if (m_arrNameVec[i] == field_name)
      return i;
  }

  return -1;
}


MultiMat::IdSet MultiMat::getMatInCell(int c)
{
  SLIC_ASSERT(m_dataLayout == DataLayout::CELL_CENTRIC);
  SLIC_ASSERT(m_cellMatRel != nullptr);

  return (*m_cellMatRel)[c];
}


MultiMat::IndexSet MultiMat::getIndexingSetOfCell(int c)
{
  SLIC_ASSERT(m_dataLayout == DataLayout::CELL_CENTRIC);
  SLIC_ASSERT(0 <= c && c < (int)m_ncells);

  if (m_sparsityLayout == SparsityLayout::SPARSE)
  {
    int start_idx = m_cellMatRel_beginsVec[c];
    int end_idx = m_cellMatRel_beginsVec[c + 1];
    return RangeSetType::SetBuilder().range(start_idx, end_idx);
  }
  else
  {
    SLIC_ASSERT(m_sparsityLayout == SparsityLayout::DENSE);
    int size2 = m_cellMatProdSet->secondSetSize();
    return RangeSetType::SetBuilder().range(c*size2, (c + 1)*size2 - 1);
  }
}

void MultiMat::convertToDynamic()
{
  if (m_dynamic_mode)
    return;
  SLIC_ASSERT(m_cellMatRelDyn == nullptr);

  // Save what the current layout is for later
  m_static_layout = Layout { m_dataLayout, m_sparsityLayout };

  // For now, handle dynamic by changing maps to dense,
  // and relation to DynamicRelation
  if (isSparse())
    convertLayoutToDense();

  //create the dynamic relation
  SetType* set1 = m_cellMatRel->fromSet();
  SetType* set2 = m_cellMatRel->toSet();

  m_cellMatRelDyn = new DynamicVariableRelationType(set1, set2);
  for (int i = 0 ; i < m_cellMatRel->fromSetSize() ; i++)
  {
    auto&& rel_vec = (*m_cellMatRel)[i];
    for (int j = 0 ; j < rel_vec.size() ; j++)
    {
      m_cellMatRelDyn->insert(i, rel_vec[j]);
    }
  }

  SLIC_ASSERT(m_cellMatRelDyn->totalSize() == m_cellMatRel->totalSize());

  delete m_cellMatRel;
  m_cellMatRel = nullptr;
  delete m_cellMatNZSet;
  m_cellMatNZSet = nullptr;

  m_dynamic_mode = true;
}


void MultiMat::convertToStatic()
{
  if (!m_dynamic_mode)
    return;

  // Change dynamicRelation back to staticRelation
  // change the layout to previously stored static layout

  RangeSetType* set1 = isCellDom() ? &m_cellSet : &m_matSet;
  RangeSetType* set2 = isCellDom() ? &m_matSet : &m_cellSet;

  SLIC_ASSERT(m_cellMatRel == nullptr);
  SLIC_ASSERT(SetPosType( m_cellMatRel_beginsVec.size()) == set1->size() + 1);
  int rel_data_size = 0;
  for (int i = 0 ; i < m_cellMatRelDyn->fromSetSize() ; i++)
  {
    auto& rel_vec = (*m_cellMatRelDyn)[i];
    m_cellMatRel_beginsVec[i] = rel_data_size;
    rel_data_size += rel_vec.size();
  }
  m_cellMatRel_beginsVec.back() = rel_data_size;
  m_cellMatRel_indicesVec.resize(rel_data_size);
  int idx = 0;
  for (int i = 0 ; i < m_cellMatRelDyn->fromSetSize() ; i++)
  {
    auto& rel_vec = (*m_cellMatRelDyn)[i];
    for (unsigned int j = 0 ; j < rel_vec.size() ; j++)
    {
      m_cellMatRel_indicesVec[idx++] = rel_vec[j];
    }
  }
  SLIC_ASSERT(idx == m_cellMatRelDyn->totalSize());

  m_cellMatRel = new StaticVariableRelationType(set1, set2);
  m_cellMatRel->bindBeginOffsets(set1->size(), &m_cellMatRel_beginsVec);
  m_cellMatRel->bindIndices(m_cellMatRel_indicesVec.size(),
                            &m_cellMatRel_indicesVec);

  SLIC_ASSERT(m_cellMatRel->isValid());

  SLIC_ASSERT(m_cellMatNZSet == nullptr);
  m_cellMatNZSet = new RelationSetType(m_cellMatRel);

  m_dynamic_mode = false;

  if (m_static_layout.sparsity_layout == SparsityLayout::SPARSE)
    convertLayoutToSparse();

  delete m_cellMatRelDyn;
  m_cellMatRelDyn = nullptr;
}

bool MultiMat::addEntry(int idx1, int idx2)
{
  SLIC_ASSERT(m_dynamic_mode);

  //right now this is implemented by converting the data layout to dense so that
  // adding an entry only involve changing the relation data.

  auto& rel_vec = m_cellMatRelDyn->data(idx1);
  auto&& found_iter = std::find(rel_vec.begin(), rel_vec.end(), idx2);
  if (found_iter != rel_vec.end())
  {
    SLIC_ASSERT_MSG(false, "MultiMat::addEntry() -- entry already exists.");
    return false;
  }

  m_cellMatRelDyn->insert(idx1, idx2);

  return true;
}

bool MultiMat::removeEntry(int idx1, int idx2)
{
  SLIC_ASSERT(m_dynamic_mode);

  auto& rel_vec = m_cellMatRelDyn->data(idx1);
  auto&& found_iter = std::find(rel_vec.begin(), rel_vec.end(), idx2);
  if (found_iter == rel_vec.end())
  {
    SLIC_ASSERT_MSG(false, "MultiMat::removeEntry() -- entry not found");
    return false;
  }

  rel_vec.erase(found_iter);

  //TODO make all map value zero?

  return true;
}


void MultiMat::convertLayoutToCellDominant()
{
  if (m_dataLayout == DataLayout::CELL_CENTRIC)
    return;
  transposeData();
}

void MultiMat::convertLayoutToMaterialDominant()
{
  if (m_dataLayout == DataLayout::MAT_CENTRIC)
    return;

  transposeData();
}

void MultiMat::transposeData()
{
  //create the new relation to pass into helper...
  RangeSetType& set1 = *(m_cellMatRel->fromSet());
  RangeSetType& set2 = *(m_cellMatRel->toSet());

  if(isCellDom())
    SLIC_ASSERT(&set1 == &m_cellSet);
  else
    SLIC_ASSERT(&set1 == &m_matSet);

  auto nz_count = m_cellMatRel_indicesVec.size();
  StaticVariableRelationType* new_cellMatRel = nullptr;
  std::vector<SetPosType> new_cellMatRel_beginsVec(set2.size() + 1, 0);
  //initialized to 0 becuase it will be used for counting
  std::vector<SetPosType> new_cellMatRel_indicesVec(nz_count, -1);
  RelationSetType* new_cellMatNZSet = nullptr;
  ProductSetType* new_cellMatProdSet = nullptr;
  std::vector<SetPosType> move_indices(nz_count, -1); //map from old to new loc

  //construct the new transposed relation
  //count the non-zero in each rows
  for (auto idx1 = 0 ; idx1 < m_cellMatRel->fromSetSize() ; ++idx1)
  {
    IdSet relSubset = (*m_cellMatRel)[idx1];
    for (auto j = 0 ; j < relSubset.size() ; ++j)
    {
      auto idx2 = relSubset[j];
      new_cellMatRel_beginsVec[idx2] += 1;
    }
  }
  //add them to make this the end index
  {
    unsigned int i;
    for( i = 1 ; i < new_cellMatRel_beginsVec.size() - 1 ; i++)
    {
      new_cellMatRel_beginsVec[i] += new_cellMatRel_beginsVec[i - 1];
    }
    new_cellMatRel_beginsVec[i] = new_cellMatRel_beginsVec[i - 1];
  }
  //fill in the indicesVec and the move_indices backward
  for (auto idx1 = m_cellMatRel->fromSetSize() - 1 ; idx1 >= 0 ; --idx1)
  {
    IdSet relSubset = (*m_cellMatRel)[idx1];
    for (auto j = relSubset.size()-1 ; j >= 0 ; --j)
    {
      auto idx2 = relSubset[j];
      auto compress_idx = --new_cellMatRel_beginsVec[idx2];
      new_cellMatRel_indicesVec[compress_idx] = idx1;
      move_indices[ m_cellMatRel_beginsVec[idx1] + j ] = compress_idx;
    }
  }

  new_cellMatRel = new StaticVariableRelationType(&set2, &set1);
  new_cellMatRel->bindBeginOffsets(set2.size(), &new_cellMatRel_beginsVec);
  new_cellMatRel->bindIndices(new_cellMatRel_indicesVec.size(),
                              &new_cellMatRel_indicesVec);

  new_cellMatNZSet = new RelationSetType(new_cellMatRel);
  new_cellMatProdSet = new ProductSetType(&set2, &set1);

  for (unsigned int map_i = 0 ; map_i < m_fieldMappingVec.size() ; map_i++)
  {
    //no conversion needed unless the field is PER_CELL_MAT
    if (m_fieldMappingVec[map_i] != FieldMapping::PER_CELL_MAT)
      continue;

    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble)
    {
      transposeData_helper<double>(map_i, new_cellMatNZSet, new_cellMatProdSet,
                                   move_indices);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat)
    {
      transposeData_helper<float>(map_i, new_cellMatNZSet, new_cellMatProdSet,
                                  move_indices);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt)
    {
      transposeData_helper<int>(map_i, new_cellMatNZSet, new_cellMatProdSet,
                                move_indices);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar)
    {
      transposeData_helper<unsigned char>(map_i, new_cellMatNZSet,
                                          new_cellMatProdSet, move_indices);
    }
    else
      SLIC_ASSERT(false);    //TODO
  }

  if(m_dataLayout == DataLayout::MAT_CENTRIC)
    m_dataLayout = DataLayout::CELL_CENTRIC;
  else
    m_dataLayout = DataLayout::MAT_CENTRIC;

  //switch vectors and rebind relation begin offset
  new_cellMatRel_beginsVec.swap(m_cellMatRel_beginsVec);
  new_cellMatRel_indicesVec.swap(m_cellMatRel_indicesVec);
  new_cellMatRel->bindBeginOffsets(set2.size(), &m_cellMatRel_beginsVec);
  new_cellMatRel->bindIndices(m_cellMatRel_indicesVec.size(),
                              &m_cellMatRel_indicesVec);

  //delete old relation and biSets
  delete m_cellMatNZSet;
  delete m_cellMatRel;
  delete m_cellMatProdSet;
  m_cellMatRel = new_cellMatRel;
  m_cellMatNZSet = new_cellMatNZSet;
  m_cellMatProdSet = new_cellMatProdSet;
}



template<typename DataType>
void MultiMat::transposeData_helper(int map_i,
                                    RelationSetType* new_cellMatNZSet,
                                    ProductSetType* new_cellMatProdSet,
                                    std::vector<SetPosType>& move_indices)
{
  MapBaseType* mapPtr = m_mapVec[map_i];

  //Skip if no volume fraction array is set-up
  if (map_i == 0 && mapPtr == nullptr)
    return;

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data;

  RangeSetType& set1 = *(m_cellMatRel->fromSet());
  RangeSetType& set2 = *(m_cellMatRel->toSet());

  int set1Size = set1.size();
  int set2Size = set2.size();

  Field2D<DataType>* new_map = nullptr;
  if (m_sparsityLayout == SparsityLayout::SPARSE)
  {
    arr_data.resize(m_cellMatRel->totalSize()*stride);
    for (unsigned int i = 0 ; i < move_indices.size() ; ++i)
    {
      for (int c = 0 ; c < stride ; ++c)
      {
        arr_data[move_indices[i] * stride + c] = (*old_map.getMap())(i, c);
      }
    }
    new_map = new Field2D<DataType>(new_cellMatNZSet, DataType(), stride);
  }
  else //dense
  {
    arr_data.resize(set1Size*set2Size*stride);
    for (int i = 0 ; i < set1Size ; ++i)
    {
      for (auto iter = old_map.begin(i) ; iter != old_map.end(i) ; ++iter)
      {
        int elem_idx = iter.index() * set1Size + i;
        for (int c = 0 ; c < stride ; ++c)
        {
          arr_data[elem_idx*stride + c] = iter.value(c);
        }
      }
    }
    new_map = new Field2D<DataType>(new_cellMatProdSet, DataType(), stride);
  }

  new_map->copy(arr_data.data());
  m_mapVec[map_i] = new_map;
  delete mapPtr;
}


void MultiMat::convertLayoutToSparse()
{
  if (m_sparsityLayout == SparsityLayout::SPARSE)
    return;

  for (unsigned int map_i = 0 ; map_i < m_fieldMappingVec.size() ; map_i++)
  {
    //no conversion needed unless the field is PER_CELL_MAT
    if (m_fieldMappingVec[map_i] != FieldMapping::PER_CELL_MAT)
      continue;

    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble)
    {
      convertToSparse_helper<double>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat)
    {
      convertToSparse_helper<float>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt)
    {
      convertToSparse_helper<int>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar)
    {
      convertToSparse_helper<unsigned char>(map_i);
    }
    else
      SLIC_ASSERT(false);    //TODO
  }
  m_sparsityLayout = SparsityLayout::SPARSE;
}


void MultiMat::convertLayoutToDense()
{
  if(m_sparsityLayout == SparsityLayout::DENSE)
    return;

  for (unsigned int map_i = 0 ; map_i < m_fieldMappingVec.size() ; map_i++)
  {
    //no conversion needed unless the field is PER_CELL_MAT
    if (m_fieldMappingVec[map_i] != FieldMapping::PER_CELL_MAT)
      continue;

    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble)
    {
      convertToDense_helper<double>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat)
    {
      convertToDense_helper<float>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt)
    {
      convertToDense_helper<int>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar)
    {
      convertToDense_helper<unsigned char>(map_i);
    }
    else
      SLIC_ASSERT(false);    //TODO
  }
  m_sparsityLayout = SparsityLayout::DENSE;
}


void MultiMat::convertLayout(DataLayout new_layout, SparsityLayout new_sparsity)
{
  if (new_layout == m_dataLayout && new_sparsity == m_sparsityLayout)
    return;

  //sparse/dense conversion
  if (m_sparsityLayout == SparsityLayout::DENSE
      && new_sparsity == SparsityLayout::SPARSE)
  {
    convertLayoutToSparse();
  }
  else if(m_sparsityLayout == SparsityLayout::SPARSE
          && new_sparsity == SparsityLayout::DENSE)
  {
    convertLayoutToDense();
  }

  //cell/mat centric conversion
  if (m_dataLayout == DataLayout::CELL_CENTRIC
      && new_layout == DataLayout::MAT_CENTRIC)
  {
    convertLayoutToMaterialDominant();
  }
  else if(m_dataLayout == DataLayout::MAT_CENTRIC
          && new_layout == DataLayout::CELL_CENTRIC)
  {
    convertLayoutToCellDominant();
  }
}


std::string MultiMat::getDataLayoutAsString() const
{
  if(isCellDom())
    return "Cell-Centric";
  else if(isMatDom())
    return "Material-Centric";
  else
    SLIC_ASSERT(false);
  return "";
}

std::string MultiMat::getSparsityLayoutAsString() const
{
  if(isSparse())
    return "Sparse";
  else if(isDense())
    return "Dense";
  else
    SLIC_ASSERT(false);
  return "";
}


void MultiMat::print() const
{
  std::stringstream sstr;

  sstr << "  Multimat Object Details:";
  sstr << "\nNumber of materials: " << m_nmats;
  sstr << "\nNumber of cells:     "<< m_ncells;
  sstr << "\nData layout:     " << getDataLayoutAsString();
  sstr << "\nSparsity layout: " << getSparsityLayoutAsString();

  sstr << "\n\n Number of fields: " << m_mapVec.size() << "\n";
  for (unsigned int i = 0 ; i < m_mapVec.size() ; i++)
  {
    sstr << "Field " << i << ": " << m_arrNameVec[i].c_str();
    sstr << "  Mapping per ";
    switch (m_fieldMappingVec[i])
    {
    case FieldMapping::PER_CELL: sstr << "cell"; break;
    case FieldMapping::PER_MAT: sstr << "material"; break;
    case FieldMapping::PER_CELL_MAT: sstr << "cellXmaterial"; break;
    }
  }
  sstr << "\n\n";

  cout << sstr.str() << endl;
}


bool MultiMat::isValid(bool verboseOutput) const
{
  bool bValid = true;
  std::stringstream errStr;

  if (m_cellSet.size() > 0 && m_matSet.size() > 0)
  {
    //It's a non-empty MultiMat object.
    //Check Volfrac exists
    if (m_mapVec[0] == nullptr)
    {
      errStr << "\n\t*No Volfrac field added.";
      bValid = false;
    }
    else
    {
      const Field2D<double>& volfrac_map =
        *dynamic_cast<Field2D<double>*>(m_mapVec[0]);

      //Check Volfrac values match the relation value (if dense)
      if (isDense() && !m_dynamic_mode)
      {
        for (int i = 0 ; i < volfrac_map.firstSetSize() ; ++i)
        {
          auto rel_iter = m_cellMatRel->begin(i);
          for (int j = 0 ; j < volfrac_map.secondSetSize() ; ++j)
          {
            double volfrac = volfrac_map[i* volfrac_map.secondSetSize() + j];
            bool zero_volfrac = volfrac == 0.0; //exact comp?
            if (rel_iter != m_cellMatRel->end(i) && *rel_iter == j)
            { //cellmat rel is present
              if (zero_volfrac)
              {
                errStr << "\n\t*Volume fraction is zero for a material "
                       << " that exists in a cell";
                bValid = false;
              }
              ++rel_iter;
            }
            else
            {
              if (!zero_volfrac)
              {
                errStr << "\n\t*Volume fraction is non-zero for a material "
                       << "not presented in a cell.";
                bValid = false;
              }
            }
          }
        }
      }

      //Check Volfrac sums to up 1.0
      std::vector<double> volfrac_sum(m_cellSet.size(), 0.0);
      for (int i = 0 ; i < volfrac_map.firstSetSize() ; ++i)
      {
        auto& submap = volfrac_map(i);
        for (int j = 0 ; j < submap.size() ; ++j)
        {
          if (isCellDom())
            volfrac_sum[i] += submap.value(j);
          else
            volfrac_sum[submap.index(j)] += submap.value(j);
        }
      }
      for (unsigned int i = 0 ; i < volfrac_sum.size() ; ++i)
      {
        if (abs(volfrac_sum[i] - 1.0) > 10e-9)
        {
          errStr << "\n\t*Volfrac does not sum to 1.0 in cell " << i;
          bValid = false;
        }
      }

    }
  }

  if (verboseOutput)
  {
    if (bValid)
    {
      errStr << "\n\t*MultiMat data was valid";
    }

    std::cout << errStr.str() << std::endl;
  }

  return bValid;
}

MultiMat::SetType* MultiMat::get_mapped_set(FieldMapping fm)
{
  SetType* set_ptr = nullptr;
  switch (fm)
  {
  case FieldMapping::PER_CELL:
    set_ptr = &m_cellSet;
    break;
  case FieldMapping::PER_MAT:
    set_ptr = &m_matSet;
    break;
  case FieldMapping::PER_CELL_MAT:
    set_ptr = dynamic_cast<SetType*>(get_mapped_biSet());
    break;
  default:
    SLIC_ASSERT(false);
    return nullptr;
  }
  return set_ptr;
}

MultiMat::BivariateSetType* MultiMat::get_mapped_biSet()
{
  BivariateSetType* set_ptr = nullptr;
  if (m_sparsityLayout == SparsityLayout::SPARSE)
    set_ptr = m_cellMatNZSet;
  else if (m_sparsityLayout == SparsityLayout::DENSE)
    set_ptr = m_cellMatProdSet;

  SLIC_ASSERT(set_ptr != nullptr);
  return set_ptr;
}
