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
#include "axom/slic.hpp"

#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>

#include <cassert>


using namespace std;
using namespace axom::multimat;


MultiMat::MultiMat(DataLayout d, SparsityLayout s) :
  m_ncells(0),
  m_nmats(0),
  m_cellMatRel(nullptr),
  m_cellMatRelDyn(nullptr),
  m_matCellRel(nullptr),
  m_matCellRelDyn(nullptr),
  m_cellMatNZSet(nullptr),
  m_cellMatProdSet(nullptr),
  m_matCellNZSet(nullptr),
  m_matCellProdSet(nullptr),
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
  delete m_matCellProdSet;
  delete m_matCellNZSet;
  delete m_cellMatRel;
  delete m_cellMatRelDyn;
  delete m_matCellRel;
  delete m_matCellRelDyn;
}

template<typename T>
MultiMat::MapBaseType* MultiMat::helper_copyField(const MultiMat& mm, int map_i)
{
  MapBaseType* other_map_ptr = mm.m_mapVec[map_i];
  if (mm.getFieldMapping(map_i) == FieldMapping::PER_CELL_MAT)
  {
    //BivariateSetType* biSetPtr = get_mapped_biSet(map_i);
    Field2D<T>* typed_ptr = dynamic_cast<Field2D<T>*>(other_map_ptr);
    
    //old field2d 
    //Field2D<T>* new_ptr = new Field2D<T>(biSetPtr, T(), typed_ptr->stride());
    //new_ptr->copy(typed_ptr->getMap()->data().data());

    Field2D<T>* new_ptr = new Field2D<T>(*typed_ptr);

    return new_ptr;
  }
  else
  {
    SetType* setPtr = get_mapped_set(map_i);
    Field1D<T>* typed_ptr = dynamic_cast<Field1D<T>*>(other_map_ptr);
    Field1D<T>* new_ptr = new Field1D<T>(setPtr, T(), typed_ptr->stride());
    new_ptr->copy(*typed_ptr);
    return new_ptr;
  }
}

// Copy constructor
MultiMat::MultiMat(const MultiMat& other) :
  m_ncells(other.m_ncells),
  m_nmats(other.m_nmats),
  m_cellSet(0,other.m_ncells),
  m_matSet(0, other.m_nmats),
  m_cellMatRel_beginsVec(other.m_cellMatRel_beginsVec),
  m_cellMatRel_indicesVec(other.m_cellMatRel_indicesVec),
  m_cellMatRel(nullptr),
  m_cellMatRelDyn(nullptr),
  m_matCellRel_beginsVec(other.m_matCellRel_beginsVec),
  m_matCellRel_indicesVec(other.m_matCellRel_indicesVec),
  m_matCellRel(nullptr),
  m_matCellRelDyn(nullptr),
  m_cellMatNZSet(nullptr),
  m_cellMatProdSet(nullptr),
  m_matCellNZSet(nullptr),
  m_matCellProdSet(nullptr),
  m_fieldNameVec(other.m_fieldNameVec),
  m_fieldMappingVec(other.m_fieldMappingVec),
  m_dataTypeVec(other.m_dataTypeVec),
  m_fieldDataLayoutVec(other.m_fieldDataLayoutVec),
  m_fieldSparsityLayoutVec(other.m_fieldSparsityLayoutVec),
  m_dynamic_mode(false)
{
  
  if (other.m_cellMatRel) {
    m_cellMatRel = new StaticVariableRelationType(&m_cellSet, &m_matSet);
    m_cellMatRel->bindBeginOffsets(m_cellSet.size(), &m_cellMatRel_beginsVec);
    m_cellMatRel->bindIndices(m_cellMatRel_indicesVec.size(),
                              &m_cellMatRel_indicesVec);
    m_cellMatNZSet = new RelationSetType(m_cellMatRel);
    m_cellMatProdSet = new ProductSetType(&m_cellSet, &m_matSet);
  }
  if (other.m_matCellRel) {
    m_matCellRel = new StaticVariableRelationType(&m_matSet, &m_cellSet);
    m_matCellRel->bindBeginOffsets(m_matSet.size(), &m_matCellRel_beginsVec);
    m_matCellRel->bindIndices(m_matCellRel_indicesVec.size(),
                              &m_matCellRel_indicesVec);
    m_matCellNZSet = new RelationSetType(m_matCellRel);
    m_matCellProdSet = new ProductSetType(&m_matSet, &m_cellSet);
  }

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

  //return *this;
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

void MultiMat::setCellMatRel(vector<bool>& vecarr, DataLayout layout)
{
  //Setup the SLAM cell to mat relation
  //This step is necessary if the volfrac field is sparse

  SLIC_ASSERT(vecarr.size() == m_ncells * m_nmats); //Check it's dense

  StaticVariableRelationType** Rel_ptr_ptr;
  std::vector<SetPosType> * Rel_beginsVec_ptr;
  std::vector<SetPosType> * Rel_indicesVec_ptr;
  if (layout == DataLayout::CELL_DOM)
  {
    Rel_ptr_ptr = &m_cellMatRel;
    Rel_beginsVec_ptr = &m_cellMatRel_beginsVec;
    Rel_indicesVec_ptr = &m_cellMatRel_indicesVec;
  }
  else
  {
    Rel_ptr_ptr = &m_matCellRel;
    Rel_beginsVec_ptr = &m_matCellRel_beginsVec;
    Rel_indicesVec_ptr = &m_matCellRel_indicesVec;
  }
  // why am I doing this? beacuse I want to use reference,
  // but don't know how to initialize them with an if-statement
  StaticVariableRelationType*& Rel_ptr = *Rel_ptr_ptr;
  std::vector<SetPosType>& Rel_beginsVec = *Rel_beginsVec_ptr;
  std::vector<SetPosType>& Rel_indicesVec = *Rel_indicesVec_ptr;

  SLIC_ASSERT(Rel_ptr == nullptr); //cellmatRel has not been set
  
  RangeSetType& set1 = (layout == DataLayout::CELL_DOM ? m_cellSet : m_matSet);
  RangeSetType& set2 = (layout == DataLayout::CELL_DOM ? m_matSet : m_cellSet);

  //count the non-zeros
  int nz_count = 0;
  for (bool b : vecarr)
    nz_count += b;

  //Set-up the cell/mat relation
  Rel_beginsVec.resize(set1.size() + 1, -1);
  Rel_indicesVec.resize(nz_count);

  SetPosType curIdx = SetPosType();
  for (SetPosType i = 0 ; i < set1.size() ; ++i)
  {
    Rel_beginsVec[i] = curIdx;
    for (SetPosType j = 0 ; j < set2.size() ; ++j)
    {
      if (vecarr[i*set2.size() + j])
      {
        Rel_indicesVec[curIdx] = j;
        ++curIdx;
      }
    }
  }
  Rel_beginsVec[set1.size()] = curIdx;

  Rel_ptr = new StaticVariableRelationType(&set1, &set2);
  Rel_ptr->bindBeginOffsets(set1.size(), &Rel_beginsVec);
  Rel_ptr->bindIndices(Rel_indicesVec.size(),
                            &Rel_indicesVec);

  SLIC_ASSERT(Rel_ptr->isValid());

  //Set-up both dense and sparse BivariateSets.
  if (layout == DataLayout::CELL_DOM)
  {
    m_cellMatNZSet = new RelationSetType(Rel_ptr);
    m_cellMatProdSet = new ProductSetType(&set1, &set2);
  }
  else
  {
    m_matCellNZSet = new RelationSetType(Rel_ptr);
    m_matCellProdSet = new ProductSetType(&set1, &set2);

  }

  //Create a field for VolFrac as the 0th field
  m_mapVec.push_back(nullptr);
  m_fieldNameVec.push_back("Volfrac");
  m_fieldMappingVec.push_back(FieldMapping::PER_CELL_MAT);
  m_dataTypeVec.push_back(DataTypeSupported::TypeDouble);
  m_fieldDataLayoutVec.push_back(DataLayout::CELL_DOM);
  m_fieldSparsityLayoutVec.push_back(SparsityLayout::SPARSE);

  SLIC_ASSERT(m_mapVec.size() == 1);
  SLIC_ASSERT(m_fieldNameVec.size() == 1);
  SLIC_ASSERT(m_fieldMappingVec.size() == 1);
  SLIC_ASSERT(m_dataTypeVec.size() == 1);
  SLIC_ASSERT(m_fieldDataLayoutVec.size() == 1);
  SLIC_ASSERT(m_fieldSparsityLayoutVec.size() == 1);

  //not efficient to have both rel, but here for debug
  makeOtherRelation();
}


int MultiMat::setVolfracField(double* arr, DataLayout layout, SparsityLayout sparsity)
{
  //m_mapVec[0] should already be a volfrac map. This functions add a new map,
  //with the input arr, then swap the new map with the 0th map,
  //and delete the new map.

  //Volfrac map is a CellxMat mapping, named "Volfrac", and is stride 1.
  int arr_i = addFieldArray_impl<double>("Volfrac", FieldMapping::PER_CELL_MAT, 
    layout, sparsity, arr, 1);

  //move the data to the first one (index 0) in the list
  std::iter_swap(m_mapVec.begin(), m_mapVec.begin() + arr_i);
  std::iter_swap(m_dataTypeVec.begin(), m_dataTypeVec.begin() + arr_i);
  m_fieldDataLayoutVec[0] = layout;
  m_fieldSparsityLayoutVec[0] = sparsity;

  //remove the new entry...
  int nfield = m_mapVec.size() - 1;
  m_mapVec.resize(nfield);
  m_fieldMappingVec.resize(nfield);
  m_fieldNameVec.resize(nfield);
  m_dataTypeVec.resize(nfield);
  m_fieldDataLayoutVec.resize(nfield);
  m_fieldSparsityLayoutVec.resize(nfield);

  return 0;
}


MultiMat::Field2D<double>& MultiMat::getVolfracField()
{
  return *dynamic_cast<Field2D<double>*>(m_mapVec[0]);
}


int MultiMat::getFieldIdx(const std::string& field_name) const
{
  for (unsigned int i = 0 ; i < m_fieldNameVec.size() ; i++)
  {
    if (m_fieldNameVec[i] == field_name)
      return i;
  }

  return -1;
}


MultiMat::IdSet MultiMat::getMatInCell(int c)
{
  //SLIC_ASSERT(m_dataLayout == DataLayout::CELL_DOM);
  SLIC_ASSERT(m_cellMatRel != nullptr);

  return (*m_cellMatRel)[c];
}


MultiMat::IdSet MultiMat::getCellContainingMat(int m)
{
  SLIC_ASSERT(m_matCellRel != nullptr);

  return (*m_matCellRel)[m];
}


MultiMat::IndexSet MultiMat::getSubfieldIndexingSet(int idx, 
  DataLayout layout, SparsityLayout sparsity)
{
  if (layout == DataLayout::CELL_DOM)
    return getIndexingSetOfCell(idx, sparsity);
  else
    return getIndexingSetOfMat(idx, sparsity);
}

MultiMat::IndexSet MultiMat::getIndexingSetOfCell(int c, SparsityLayout sparsity)
{
  SLIC_ASSERT(0 <= c && c < (int)m_ncells);

  if (sparsity == SparsityLayout::SPARSE)
  {
    int start_idx = m_cellMatRel_beginsVec[c];
    int end_idx = m_cellMatRel_beginsVec[c + 1];
    return RangeSetType(start_idx, end_idx);
  }
  else
  {
    SLIC_ASSERT(sparsity == SparsityLayout::DENSE);
    int size2 = m_cellMatProdSet->secondSetSize();
    return RangeSetType(c*size2, (c + 1)*size2);
  }
}

MultiMat::IndexSet MultiMat::getIndexingSetOfMat(int m, SparsityLayout sparsity)
{
  //SLIC_ASSERT(m_dataLayout == DataLayout::MAT_DOM);
  SLIC_ASSERT(0 <= m && m < (int)m_nmats);

  if (sparsity == SparsityLayout::SPARSE)
  {
    int start_idx = m_matCellRel_beginsVec[m];
    int end_idx = m_matCellRel_beginsVec[m + 1];
    return RangeSetType::SetBuilder().range(start_idx, end_idx);
  }
  else
  {
    SLIC_ASSERT(sparsity == SparsityLayout::DENSE);
    int size2 = m_matCellProdSet->secondSetSize();
    return RangeSetType::SetBuilder().range(m * size2, (m + 1) * size2 - 1);
  }
}


void MultiMat::convertToDynamic()
{
  SLIC_ASSERT(false); //need to fix
                      /*
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
  */
}


void MultiMat::convertToStatic()
{
  SLIC_ASSERT(false); //need to fix
  /*
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
  */
}

bool MultiMat::addEntry(int idx1, int idx2)
{
  SLIC_ASSERT(false); //need to fix
                      /*
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
  */
  return true;
}

bool MultiMat::removeEntry(int idx1, int idx2)
{
  SLIC_ASSERT(false); //need to fix
                      /*
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
  */
  return true;
}

void MultiMat::makeOtherRelation()
{
  StaticVariableRelationType** oldRelptr, **newRelptr;
  std::vector<SetPosType>* newBeginVec_ptr, *newIndicesVec_ptr;
  if (m_cellMatRel == nullptr && m_matCellRel != nullptr)
  {
    oldRelptr = &m_matCellRel;
    newRelptr = &m_cellMatRel;
    newBeginVec_ptr = &m_cellMatRel_beginsVec;
    newIndicesVec_ptr = &m_cellMatRel_indicesVec;
  }
  else if (m_matCellRel == nullptr && m_cellMatRel != nullptr)
  {
    oldRelptr = &m_cellMatRel;
    newRelptr = &m_matCellRel;
    newBeginVec_ptr = &m_matCellRel_beginsVec;
    newIndicesVec_ptr = &m_matCellRel_indicesVec;
  }
  else SLIC_ASSERT(false);

  StaticVariableRelationType *& oldRel = *oldRelptr;
  StaticVariableRelationType *& newRel = *newRelptr;
  std::vector<SetPosType>& newBeginVec = *newBeginVec_ptr;
  std::vector<SetPosType>& newIndicesVec = *newIndicesVec_ptr;

  RangeSetType& set1 = *(oldRel->fromSet());
  RangeSetType& set2 = *(oldRel->toSet());
  
  auto nz_count = oldRel->totalSize();
  //auto nz_count = m_cellMatRel_indicesVec.size();

  newBeginVec.resize(set2.size() + 1, 0);
  newIndicesVec.resize(nz_count, -1);

  //construct the new transposed relation
                                                      
  //count the non-zero in each rows
  for (auto idx1 = 0; idx1 < oldRel->fromSetSize(); ++idx1)
  {
    IdSet relSubset = (*oldRel)[idx1];
    for (auto j = 0; j < relSubset.size(); ++j)
    {
      auto idx2 = relSubset[j];
      newBeginVec[idx2] += 1;
    }
  }

  //add them to make this the end index
  {
    unsigned int i;
    for (i = 1; i < newBeginVec.size() - 1; i++)
    {
      newBeginVec[i] += newBeginVec[i - 1];
    }
    newBeginVec[i] = newBeginVec[i - 1];
  }

  //fill in the indicesVec and the move_indices backward
  for (auto idx1 = oldRel->fromSetSize() - 1; idx1 >= 0; --idx1)
  {
    IdSet relSubset = (*oldRel)[idx1];
    for (auto j = relSubset.size() - 1; j >= 0; --j)
    {
      auto idx2 = relSubset[j];
      auto compress_idx = --newBeginVec[idx2];
      newIndicesVec[compress_idx] = idx1;
    }
  }

  newRel = new StaticVariableRelationType(&set2, &set1);
  newRel->bindBeginOffsets(set2.size(), &newBeginVec);
  newRel->bindIndices(newIndicesVec.size(), &newIndicesVec);

  if (newRelptr == &m_matCellRel)
  {
    m_matCellNZSet = new RelationSetType(newRel);
    m_matCellProdSet = new ProductSetType(&set2, &set1);
  }
  else
  {
    m_cellMatNZSet = new RelationSetType(newRel);
    m_cellMatProdSet = new ProductSetType(&set2, &set1);
  }

}

void MultiMat::convertLayoutToCellDominant()
{
  for (unsigned int i = 0; i < m_mapVec.size(); ++i)
  {
    convertFieldToCellDom(i);
  }
}

void MultiMat::convertLayoutToMaterialDominant()
{
  for (unsigned int i = 0; i < m_mapVec.size(); ++i)
  {
    convertFieldToMatDom(i);
  }
}

void MultiMat::convertLayoutToSparse()
{
  for (unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    convertFieldToSparse(map_i);
  }
}

void MultiMat::convertLayoutToDense()
{
  for (unsigned int map_i = 0 ; map_i < m_fieldMappingVec.size() ; map_i++)
  {
    convertFieldToDense(map_i);
  }
}

void MultiMat::convertLayout(DataLayout new_layout, SparsityLayout new_sparsity)
{

  //sparse/dense conversion
  if ( new_sparsity == SparsityLayout::SPARSE)
  {
    convertLayoutToSparse();
  }
  else if( new_sparsity == SparsityLayout::DENSE)
  {
    convertLayoutToDense();
  }

  //cell/mat centric conversion
  if ( new_layout == DataLayout::MAT_DOM)
  {
    convertLayoutToMaterialDominant();
  }
  else if( new_layout == DataLayout::CELL_DOM)
  {
    convertLayoutToCellDominant();
  }
}

void MultiMat::convertFieldLayout(int field_idx, SparsityLayout new_sparsity, DataLayout new_layout)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < m_mapVec.size());

  DataLayout field_data_layout = m_fieldDataLayoutVec[field_idx];
  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  if (new_layout == field_data_layout && new_sparsity == field_sparsity_layout)
    return;

  //sparse/dense conversion
  if (field_sparsity_layout == SparsityLayout::DENSE
    && new_sparsity == SparsityLayout::SPARSE)
  {
    convertFieldToSparse(field_idx);
  }
  else if (field_sparsity_layout == SparsityLayout::SPARSE
    && new_sparsity == SparsityLayout::DENSE)
  {
    convertFieldToDense(field_idx);
  }

  //cell/mat centric conversion
  if (field_data_layout == DataLayout::CELL_DOM
    && new_layout == DataLayout::MAT_DOM)
  {
    convertFieldToMatDom(field_idx);
  }
  else if (field_data_layout == DataLayout::MAT_DOM
    && new_layout == DataLayout::CELL_DOM)
  {
    convertFieldToCellDom(field_idx);
  }

  SLIC_ASSERT(false);
}

void axom::multimat::MultiMat::convertFieldToSparse(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < m_mapVec.size());

  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  //no conversion needed unless the field is PER_CELL_MAT
  if (field_sparsity_layout == SparsityLayout::SPARSE ||
    m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT) {
    return;
  }

  if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeDouble)
  {
    convertToSparse_helper<double>(field_idx);
  }
  else if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeFloat)
  {
    convertToSparse_helper<float>(field_idx);
  }
  else if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeInt)
  {
    convertToSparse_helper<int>(field_idx);
  }
  else if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeUnsignChar)
  {
    convertToSparse_helper<unsigned char>(field_idx);
  }
  else
    SLIC_ASSERT(false);    //TODO

  m_fieldSparsityLayoutVec[field_idx] = SparsityLayout::SPARSE;
}

void axom::multimat::MultiMat::convertFieldToDense(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < m_mapVec.size());

  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  //no conversion needed unless the field is PER_CELL_MAT
  if (field_sparsity_layout == SparsityLayout::DENSE ||
    m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT) {
    return;
  }

  if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeDouble)
  {
    convertToDense_helper<double>(field_idx);
  }
  else if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeFloat)
  {
    convertToDense_helper<float>(field_idx);
  }
  else if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeInt)
  {
    convertToDense_helper<int>(field_idx);
  }
  else if (m_dataTypeVec[field_idx] == DataTypeSupported::TypeUnsignChar)
  {
    convertToDense_helper<unsigned char>(field_idx);
  }
  else
    SLIC_ASSERT(false);    //TODO

  m_fieldSparsityLayoutVec[field_idx] = SparsityLayout::DENSE;
}


template<typename DataType>
void MultiMat::convertToSparse_helper(int map_i)
{
  SLIC_ASSERT(m_fieldSparsityLayoutVec[map_i] != SparsityLayout::SPARSE);

  MapBaseType* mapPtr = m_mapVec[map_i];

  //Skip if no volume fraction array is set-up
  if (map_i == 0 && mapPtr == nullptr)
    return;

  StaticVariableRelationType* Rel = getRel(map_i);

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(Rel->totalSize()*stride);
  int idx = 0;
  for (int i = 0; i < Rel->fromSetSize(); ++i)
  {
    auto relset = (*Rel)[i];
    auto submap = old_map(i);
    for (int j = 0; j < relset.size(); ++j)
    {
      for (int s = 0; s < stride; ++s)
      {
        arr_data[idx++] = submap[relset[j] * stride + s];
      }
    }
  }
  SLIC_ASSERT(idx == Rel->totalSize()*stride);

  RelationSetType* nz_set = m_fieldDataLayoutVec[map_i] == DataLayout::CELL_DOM ?
    m_cellMatNZSet : m_matCellNZSet;
  //old field2d
  //Field2D<DataType>* new_field = new Field2D<DataType>(nz_set, DataType(), stride);
  //new_field->copy(arr_data.data());
  Field2D<DataType>* new_field = new Field2D<DataType>(*this, nz_set,
    old_map.getName(), arr_data.data(), stride);

  delete m_mapVec[map_i];
  m_mapVec[map_i] = new_field;
}


template<typename DataType>
void MultiMat::convertToDense_helper(int map_i)
{
  SLIC_ASSERT(m_fieldSparsityLayoutVec[map_i] != SparsityLayout::DENSE);

  MapBaseType* mapPtr = m_mapVec[map_i];

  //Skip if no volume fraction array is set-up
  if (map_i == 0 && mapPtr == nullptr)
    return;

  ProductSetType* prod_set = m_fieldDataLayoutVec[map_i] == DataLayout::CELL_DOM ?
    m_cellMatProdSet : m_matCellProdSet;

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(prod_set->size()*stride);
  for (int i = 0; i < old_map.firstSetSize(); ++i)
  {
    for (auto iter = old_map.begin(i); iter != old_map.end(i); ++iter)
    {
      int elem_idx = i * old_map.secondSetSize() + iter.index();
      for (int c = 0; c < stride; ++c)
      {
        arr_data[elem_idx*stride + c] = iter.value(c);
      }
    }
  }

  //old field2d
  //Field2D<DataType>* new_field = new Field2D<DataType>(prod_set, DataType(), stride);
  //new_field->copy(&arr_data[0]);
  Field2D<DataType>* new_field = new Field2D<DataType>(*this, prod_set,
    old_map.getName(), arr_data.data(), stride);

  delete m_mapVec[map_i];
  m_mapVec[map_i] = new_field;
}


DataLayout MultiMat::getFieldDataLayout(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < m_mapVec.size());

  return m_fieldDataLayoutVec[field_idx]; 
}

SparsityLayout MultiMat::getFieldSparsityLayout(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < m_mapVec.size());

  return m_fieldSparsityLayoutVec[field_idx];
}


template<typename DataType>
void MultiMat::transposeField_helper(int field_idx)
{
  SLIC_ASSERT(m_fieldMappingVec[field_idx] == FieldMapping::PER_CELL_MAT);

  MapBaseType* oldMapPtr = m_mapVec[field_idx];

  //Skip if no volume fraction array is set-up
  if (field_idx == 0 && oldMapPtr == nullptr)
    return;

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(oldMapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data;
  
  DataLayout oldDataLayout = getFieldDataLayout(field_idx);
  StaticVariableRelationType* oldRel = nullptr;
  StaticVariableRelationType* newRel = nullptr;
  RelationSetType* newNZSet = nullptr;
  ProductSetType* newProdSet = nullptr;
  DataLayout new_layout;
  if (oldDataLayout == DataLayout::CELL_DOM)
  {
    oldRel = m_cellMatRel;
    if (!m_matCellRel) makeOtherRelation();
    newRel = m_matCellRel;
    newNZSet = m_matCellNZSet;
    newProdSet = m_matCellProdSet;
    new_layout = DataLayout::MAT_DOM;
  }
  else
  {
    oldRel = m_matCellRel;
    if (!m_cellMatRel) makeOtherRelation();
    newRel = m_cellMatRel;
    newNZSet = m_cellMatNZSet;
    newProdSet = m_cellMatProdSet;
    new_layout = DataLayout::CELL_DOM;
  }

  SLIC_ASSERT(oldRel != nullptr);
  SLIC_ASSERT(newRel != nullptr);
  SLIC_ASSERT((newNZSet != nullptr) || (newProdSet != nullptr));

  auto& set1 = *(oldRel->fromSet());
  auto& set2 = *(oldRel->toSet());

  int set1Size = set1.size();
  int set2Size = set2.size();

  Field2D<DataType>* new_map = nullptr;
  if (m_fieldSparsityLayoutVec[field_idx] == SparsityLayout::SPARSE)
  {
    //copy begin vector for moving
    std::vector<SetPosType> vec_idx;  //a copy of the beginVec to keep track
    const std::vector<SetPosType>& indicesVec = *oldRel->relationData();
    if (oldDataLayout == DataLayout::CELL_DOM)
      vec_idx = m_matCellRel_beginsVec;
    else
      vec_idx = m_cellMatRel_beginsVec;

    arr_data.resize(oldRel->totalSize() * stride);
    for (int i = 0; i < oldRel->totalSize(); ++i)
    {
      int col = indicesVec[i];
      for (int c = 0; c < stride; ++c)
      {
        arr_data[vec_idx[col] * stride + c] = (*old_map.getMap())(i, c);
      }
      ++vec_idx[col];
    }

    //old
    //new_map = new Field2D<DataType>(newNZSet, DataType(), stride);
    new_map = new Field2D<DataType>(*this, newNZSet,
      old_map.getName(), arr_data.data(), stride);
  }
  else //dense
  {
    arr_data.resize(set1Size * set2Size * stride);
    for (int i = 0; i < set1Size; ++i)
    {
      for (auto iter = old_map.begin(i); iter != old_map.end(i); ++iter)
      {
        int elem_idx = iter.index() * set1Size + i;
        for (int c = 0; c < stride; ++c)
        {
          arr_data[elem_idx * stride + c] = iter.value(c);
        }
      }
    }
    //new_map = new Field2D<DataType>(newProdSet, DataType(), stride);
    new_map = new Field2D<DataType>(*this, newProdSet, old_map.getName(),
      arr_data.data(), stride);
  }

  //new_map->copy(arr_data.data());
  m_mapVec[field_idx] = new_map;
  delete oldMapPtr;
  m_fieldDataLayoutVec[field_idx] = new_layout;
  
}



void axom::multimat::MultiMat::transposeField(int field_idx)
{
  if (m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT)
    return;

  switch (m_dataTypeVec[field_idx])
  {
  case DataTypeSupported::TypeDouble:
    transposeField_helper<double>(field_idx);
    break;
  case DataTypeSupported::TypeFloat:
    transposeField_helper<float>(field_idx);
    break;
  case DataTypeSupported::TypeInt:
    transposeField_helper<int>(field_idx);
    break;
  case DataTypeSupported::TypeUnsignChar:
    transposeField_helper<unsigned char>(field_idx);
    break;
  default:
    SLIC_ASSERT(false);
  }
}

void MultiMat::convertFieldToMatDom(int field_idx)
{
  if (getFieldDataLayout(field_idx) == DataLayout::MAT_DOM)
  {
    return;
  }
  transposeField(field_idx);
}

void MultiMat::convertFieldToCellDom(int field_idx)
{
  if (getFieldDataLayout(field_idx) == DataLayout::CELL_DOM)
  {
    return;
  }
  transposeField(field_idx);
}


std::string MultiMat::getFieldDataLayoutAsString(int field_i) const
{
  if (m_fieldDataLayoutVec[field_i] == DataLayout::CELL_DOM)
    return "Cell-Centric";
  else if (m_fieldDataLayoutVec[field_i] == DataLayout::MAT_DOM)
    return "Material-Centric";
  else
    SLIC_ASSERT(false);
  return "";
}

std::string MultiMat::getFieldSparsityLayoutAsString(int field_i) const
{
  if (m_fieldSparsityLayoutVec[field_i] == SparsityLayout::SPARSE)
    return "Sparse";
  else if (m_fieldSparsityLayoutVec[field_i] == SparsityLayout::DENSE)
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

  sstr << "\n\n Number of fields: " << m_mapVec.size() << "\n";
  for (unsigned int i = 0 ; i < m_mapVec.size() ; i++)
  {
    sstr << "Field " << i << ": " << m_fieldNameVec[i].c_str();
    sstr << "  Mapping per ";
    switch (m_fieldMappingVec[i])
    {
    case FieldMapping::PER_CELL: sstr << "cell"; break;
    case FieldMapping::PER_MAT: sstr << "material"; break;
    case FieldMapping::PER_CELL_MAT: 
      sstr << "cellXmaterial"; 
      sstr << "\n  Data layout: " << getFieldDataLayoutAsString(i);
      sstr << "\n  Sparsity layout: " << getFieldSparsityLayoutAsString(i);
      break;
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
      const auto& volfrac_map = *dynamic_cast<Field2D<double>*>(m_mapVec[0]);
      auto volfrac_sparsity = m_fieldSparsityLayoutVec[0];
      auto volfrac_layout = m_fieldDataLayoutVec[0];

      //Check Volfrac values match the relation value (if dense)
      if (volfrac_sparsity==SparsityLayout::DENSE && !m_dynamic_mode)
      {
        const StaticVariableRelationType* relPtr = (volfrac_layout ==
          DataLayout::CELL_DOM ? m_cellMatRel : m_matCellRel);

        for (int i = 0 ; i < volfrac_map.firstSetSize() ; ++i)
        {
          auto rel_iter = relPtr->begin(i);
          for (int j = 0 ; j < volfrac_map.secondSetSize() ; ++j)
          {
            double volfrac = volfrac_map[i* volfrac_map.secondSetSize() + j];
            bool zero_volfrac = volfrac == 0.0; //exact comp?
            if (rel_iter != relPtr->end(i) && *rel_iter == j)
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
        const auto& submap = volfrac_map(i);
        for (int j = 0; j < submap.size(); ++j)
        {
          if (volfrac_layout == DataLayout::CELL_DOM)
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
    SLIC_ASSERT(false);
    //dynamic_cast<SetType*>(get_mapped_biSet(m_dataLayout, m_sparsityLayout));
    //todo fix this haha
    break;
  default:
    SLIC_ASSERT(false);
    return nullptr;
  }
  return set_ptr;
}

MultiMat::SetType * MultiMat::get_mapped_set(int field_idx)
{
  SetType* set_ptr = nullptr;
  switch (m_fieldMappingVec[field_idx])
  {
  case FieldMapping::PER_CELL:
    set_ptr = &m_cellSet;
    break;
  case FieldMapping::PER_MAT:
    set_ptr = &m_matSet;
    break;
  case FieldMapping::PER_CELL_MAT:
    set_ptr = dynamic_cast<SetType*>(get_mapped_biSet(field_idx));
    break;
  default:
    SLIC_ASSERT(false);
    return nullptr;
  }
  return set_ptr;
}

MultiMat::BivariateSetType* MultiMat::get_mapped_biSet(
  DataLayout layout, SparsityLayout sparsity)
{
  BivariateSetType* set_ptr = nullptr;

  if (sparsity == SparsityLayout::SPARSE)
    if (layout == DataLayout::CELL_DOM)
      set_ptr = m_cellMatNZSet;
    else
      set_ptr = m_matCellNZSet;
  else if (sparsity == SparsityLayout::DENSE)
    if (layout == DataLayout::CELL_DOM)
      set_ptr = m_cellMatProdSet;
    else
      set_ptr = m_matCellProdSet;

  SLIC_ASSERT(set_ptr != nullptr);
  return set_ptr;
}

MultiMat::BivariateSetType * MultiMat::get_mapped_biSet(int field_idx)
{
  SLIC_ASSERT(m_fieldMappingVec[field_idx] == FieldMapping::PER_CELL_MAT);
  SLIC_ASSERT(0 <= field_idx && field_idx < m_fieldNameVec.size());

  DataLayout layout = m_fieldDataLayoutVec[field_idx];
  SparsityLayout sparsity = m_fieldSparsityLayoutVec[field_idx];

  return get_mapped_biSet(layout, sparsity);
}

MultiMat::StaticVariableRelationType * MultiMat::getRel(int field_idx)
{
  //get the sparse relation vector
  MultiMat::StaticVariableRelationType * rel_ptr;
  if (m_fieldDataLayoutVec[field_idx] == DataLayout::CELL_DOM)
    rel_ptr = m_cellMatRel;
  else if (m_fieldDataLayoutVec[field_idx] == DataLayout::MAT_DOM)
    rel_ptr = m_matCellRel;

  SLIC_ASSERT(rel_ptr != nullptr);
  return rel_ptr;
}

