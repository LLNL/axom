// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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

MultiMat::MultiMat(DataLayout AXOM_UNUSED_PARAM(d),
                   SparsityLayout AXOM_UNUSED_PARAM(s))
  : m_ncells(0)
  , m_nmats(0)
  , m_sets(2)
  , m_staticRelations(2)
  , m_dynamicRelations(2)
  , m_sparseBivarSet(2)
  , m_denseBivarSet(2)
  , m_dynamic_mode(false)
{ }

template <typename T>
MultiMat::MapUniquePtr MultiMat::helper_copyField(const MultiMat& mm, int map_i)
{
  MapBaseType* other_map_ptr = mm.m_mapVec[map_i].get();
  if(mm.getFieldMapping(map_i) == FieldMapping::PER_CELL_MAT)
  {
    //BivariateSetType* biSetPtr = get_mapped_biSet(map_i);
    Field2D<T>* typed_ptr = dynamic_cast<Field2D<T>*>(other_map_ptr);

    //old field2d with variable stride
    //Field2D<T>* new_ptr = new Field2D<T>(biSetPtr, T(), typed_ptr->stride());
    //new_ptr->copy(typed_ptr->getMap()->data().data());

    //fixed strideOne (in the definition of Field2D)
    Field2D<T>* new_ptr = new Field2D<T>(*typed_ptr);

    return MapUniquePtr {new_ptr};
  }
  else
  {
    const RangeSetType& setPtr =
      *static_cast<RangeSetType*>(get_mapped_set(map_i));
    Field1D<T>* typed_ptr = dynamic_cast<Field1D<T>*>(other_map_ptr);
    Field1D<T>* new_ptr = new Field1D<T>(setPtr, T(), typed_ptr->stride());
    new_ptr->copy(*typed_ptr);
    return MapUniquePtr {new_ptr};
  }
}

MultiMat::IndBufferType& MultiMat::relBeginVec(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? m_cellMatRel_beginsVec
                                          : m_matCellRel_beginsVec;
}

MultiMat::IndBufferType& MultiMat::relIndVec(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? m_cellMatRel_indicesVec
                                          : m_matCellRel_indicesVec;
}

MultiMat::StaticVariableRelationType& MultiMat::relStatic(DataLayout layout)
{
  return m_staticRelations[(int)layout];
}

MultiMat::DynamicVariableRelationType& MultiMat::relDynamic(DataLayout layout)
{
  return m_dynamicRelations[(int)layout];
}

MultiMat::RangeSetType& MultiMat::relDominantSet(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? getCellSet() : getMatSet();
}

MultiMat::RangeSetType& MultiMat::relSecondarySet(DataLayout layout)
{
  return (layout == DataLayout::CELL_DOM) ? getMatSet() : getCellSet();
}

MultiMat::RelationSetType& MultiMat::relSparseSet(DataLayout layout)
{
  return m_sparseBivarSet[(int)layout];
}

MultiMat::ProductSetType& MultiMat::relDenseSet(DataLayout layout)
{
  return m_denseBivarSet[(int)layout];
}

bool MultiMat::hasValidStaticRelation(DataLayout layout) const
{
  const StaticVariableRelationType& rel = m_staticRelations[(int)layout];

  return rel.hasFromSet() && rel.hasToSet();
}

bool MultiMat::hasValidDynamicRelation(DataLayout layout) const
{
  const DynamicVariableRelationType& rel = m_dynamicRelations[(int)layout];

  return rel.hasFromSet() && rel.hasToSet();
}

// Copy constructor
MultiMat::MultiMat(const MultiMat& other)
  : m_ncells(other.m_ncells)
  , m_nmats(other.m_nmats)
  , m_sets(other.m_sets)
  , m_cellMatRel_beginsVec(other.m_cellMatRel_beginsVec)
  , m_cellMatRel_indicesVec(other.m_cellMatRel_indicesVec)
  , m_matCellRel_beginsVec(other.m_matCellRel_beginsVec)
  , m_matCellRel_indicesVec(other.m_matCellRel_indicesVec)
  , m_staticRelations(other.m_staticRelations)
  , m_dynamicRelations(other.m_dynamicRelations)
  , m_sparseBivarSet(other.m_sparseBivarSet)
  , m_denseBivarSet(other.m_denseBivarSet)
  , m_fieldNameVec(other.m_fieldNameVec)
  , m_fieldMappingVec(other.m_fieldMappingVec)
  , m_dataTypeVec(other.m_dataTypeVec)
  , m_fieldDataLayoutVec(other.m_fieldDataLayoutVec)
  , m_fieldSparsityLayoutVec(other.m_fieldSparsityLayoutVec)
  , m_dynamic_mode(false)
{
  if(other.hasValidStaticRelation(DataLayout::CELL_DOM))
  {
    StaticVariableRelationType& cellMatRel = relStatic(DataLayout::CELL_DOM);

    cellMatRel = StaticVariableRelationType(&getCellSet(), &getMatSet());
    cellMatRel.bindBeginOffsets(getCellSet().size(), &m_cellMatRel_beginsVec);
    cellMatRel.bindIndices(m_cellMatRel_indicesVec.size(),
                           &m_cellMatRel_indicesVec);
    relSparseSet(DataLayout::CELL_DOM) = RelationSetType(&cellMatRel);
    relDenseSet(DataLayout::CELL_DOM) =
      ProductSetType(&getCellSet(), &getMatSet());
  }
  if(other.hasValidStaticRelation(DataLayout::MAT_DOM))
  {
    StaticVariableRelationType& matCellRel = relStatic(DataLayout::MAT_DOM);

    matCellRel = StaticVariableRelationType(&getMatSet(), &getCellSet());
    matCellRel.bindBeginOffsets(getMatSet().size(), &m_matCellRel_beginsVec);
    matCellRel.bindIndices(m_matCellRel_indicesVec.size(),
                           &m_matCellRel_indicesVec);
    relSparseSet(DataLayout::MAT_DOM) = RelationSetType(&matCellRel);
    relDenseSet(DataLayout::MAT_DOM) =
      ProductSetType(&getMatSet(), &getCellSet());
  }

  for(unsigned int map_i = 0; map_i < other.m_mapVec.size(); ++map_i)
  {
    MapUniquePtr new_map_ptr {nullptr};
    if(m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble)
    {
      new_map_ptr = helper_copyField<double>(other, map_i);
    }
    else if(m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat)
    {
      new_map_ptr = helper_copyField<float>(other, map_i);
    }
    else if(m_dataTypeVec[map_i] == DataTypeSupported::TypeInt)
    {
      new_map_ptr = helper_copyField<int>(other, map_i);
    }
    else if(m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar)
    {
      new_map_ptr = helper_copyField<unsigned char>(other, map_i);
    }
    else
      SLIC_ASSERT_MSG(false,
                      "\t*MultiMat copy constructor : Unsupported Datatype");

    SLIC_ASSERT(new_map_ptr != nullptr);
    m_mapVec.push_back(std::move(new_map_ptr));
  }
}

void MultiMat::setNumberOfMaterials(int n)
{
  SLIC_ASSERT(n > 0);
  m_nmats = n;

  getMatSet() = RangeSetType(0, m_nmats);
  SLIC_ASSERT(getMatSet().isValid());
}

void MultiMat::setNumberOfCells(int c)
{
  SLIC_ASSERT(c > 0);
  m_ncells = c;

  getCellSet() = RangeSetType(0, m_ncells);
  SLIC_ASSERT(getCellSet().isValid());
}

void MultiMat::setCellMatRel(vector<bool>& vecarr, DataLayout layout)
{
  //Setup the SLAM cell to mat relation
  //This step is necessary if the volfrac field is sparse

  SLIC_ASSERT(vecarr.size() == m_ncells * m_nmats);  //Check it's dense

  StaticVariableRelationType& Rel_ptr = relStatic(layout);
  IndBufferType& Rel_beginsVec = relBeginVec(layout);
  IndBufferType& Rel_indicesVec = relIndVec(layout);

  SLIC_ASSERT(!hasValidStaticRelation(layout));

  RangeSetType& set1 = relDominantSet(layout);
  RangeSetType& set2 = relSecondarySet(layout);

  //count the non-zeros
  int nz_count = 0;
  for(bool b : vecarr) nz_count += b;

  //Set-up the cell/mat relation
  Rel_beginsVec.resize(set1.size() + 1, -1);
  Rel_indicesVec.resize(nz_count);

  SetPosType curIdx = SetPosType();
  for(SetPosType i = 0; i < set1.size(); ++i)
  {
    Rel_beginsVec[i] = curIdx;
    for(SetPosType j = 0; j < set2.size(); ++j)
    {
      if(vecarr[i * set2.size() + j])
      {
        Rel_indicesVec[curIdx] = j;
        ++curIdx;
      }
    }
  }
  Rel_beginsVec[set1.size()] = curIdx;

  Rel_ptr = StaticVariableRelationType(&set1, &set2);
  Rel_ptr.bindBeginOffsets(set1.size(), &Rel_beginsVec);
  Rel_ptr.bindIndices(Rel_indicesVec.size(), &Rel_indicesVec);

  SLIC_ASSERT(Rel_ptr.isValid());

  //Set-up both dense and sparse BivariateSets.
  relSparseSet(layout) = RelationSetType(&Rel_ptr);
  relDenseSet(layout) = ProductSetType(&set1, &set2);

  //Create a field for VolFrac as the 0th field
  m_mapVec.emplace_back(nullptr);
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
}

int MultiMat::setVolfracField(double* arr,
                              DataLayout layout,
                              SparsityLayout sparsity)
{
  //m_mapVec[0] should already be a volfrac map. This functions add a new map,
  //with the input arr, then swap the new map with the 0th map,
  //and delete the new map.

  //Volfrac map is a CellxMat mapping, named "Volfrac", and is stride 1.
  int arr_i = addFieldArray_impl<double>("Volfrac",
                                         FieldMapping::PER_CELL_MAT,
                                         layout,
                                         sparsity,
                                         arr,
                                         1);

  //move the data to the first one (index 0) in the list
  std::iter_swap(m_mapVec.begin(), m_mapVec.begin() + arr_i);
  std::iter_swap(m_dataTypeVec.begin(), m_dataTypeVec.begin() + arr_i);
  m_fieldDataLayoutVec[0] = layout;
  m_fieldSparsityLayoutVec[0] = sparsity;

  //remove the new entry...
  int nfield = m_mapVec.size() - 1;
  m_mapVec.erase(m_mapVec.begin() + nfield);
  m_fieldMappingVec.resize(nfield);
  m_fieldNameVec.resize(nfield);
  m_dataTypeVec.resize(nfield);
  m_fieldDataLayoutVec.resize(nfield);
  m_fieldSparsityLayoutVec.resize(nfield);

  return 0;
}

MultiMat::Field2D<double>& MultiMat::getVolfracField()
{
  return *dynamic_cast<Field2D<double>*>(m_mapVec[0].get());
}

int MultiMat::getFieldIdx(const std::string& field_name) const
{
  for(unsigned int i = 0; i < m_fieldNameVec.size(); i++)
  {
    if(m_fieldNameVec[i] == field_name) return i;
  }

  return -1;
}

MultiMat::IdSet MultiMat::getMatInCell(int c)
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::CELL_DOM));

  return relStatic(DataLayout::CELL_DOM)[c];
}

MultiMat::IdSet MultiMat::getCellContainingMat(int m)
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::MAT_DOM));

  return relStatic(DataLayout::MAT_DOM)[m];
}

MultiMat::IndexSet MultiMat::getSubfieldIndexingSet(int idx,
                                                    DataLayout layout,
                                                    SparsityLayout sparsity)
{
  if(layout == DataLayout::CELL_DOM)
    return getIndexingSetOfCell(idx, sparsity);
  else
    return getIndexingSetOfMat(idx, sparsity);
}

MultiMat::IndexSet MultiMat::getIndexingSetOfCell(int c, SparsityLayout sparsity)
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::CELL_DOM));
  SLIC_ASSERT(0 <= c && c < (int)m_ncells);

  if(sparsity == SparsityLayout::SPARSE)
  {
    int start_idx = m_cellMatRel_beginsVec[c];
    int end_idx = m_cellMatRel_beginsVec[c + 1];
    return RangeSetType(start_idx, end_idx);
  }
  else
  {
    SLIC_ASSERT(sparsity == SparsityLayout::DENSE);
    int size2 = relDenseSet(DataLayout::CELL_DOM).secondSetSize();
    return RangeSetType(c * size2, (c + 1) * size2);
  }
}

MultiMat::IndexSet MultiMat::getIndexingSetOfMat(int m, SparsityLayout sparsity)
{
  SLIC_ASSERT(hasValidStaticRelation(DataLayout::MAT_DOM));
  SLIC_ASSERT(0 <= m && m < (int)m_nmats);

  if(sparsity == SparsityLayout::SPARSE)
  {
    int start_idx = m_matCellRel_beginsVec[m];
    int end_idx = m_matCellRel_beginsVec[m + 1];
    return RangeSetType::SetBuilder().range(start_idx, end_idx);
  }
  else
  {
    SLIC_ASSERT(sparsity == SparsityLayout::DENSE);
    int size2 = relDenseSet(DataLayout::MAT_DOM).secondSetSize();
    return RangeSetType::SetBuilder().range(m * size2, (m + 1) * size2 - 1);
  }
}

void MultiMat::convertToDynamic()
{
  if(m_dynamic_mode)
  {
    return;
  }
  // If we weren't in dynamic mode before, we shouldn't already have dynamic
  // relations
  SLIC_ASSERT(!hasValidDynamicRelation(DataLayout::CELL_DOM));
  SLIC_ASSERT(!hasValidDynamicRelation(DataLayout::MAT_DOM));

  // Save what the current layout is for later
  //m_static_layout = Layout { m_dataLayout, m_sparsityLayout };  //old version with single layout for all MM
  m_layout_when_static.resize(m_fieldDataLayoutVec.size());
  const int SZ = m_fieldDataLayoutVec.size();
  for(auto i = 0; i < SZ; i++)
  {
    m_layout_when_static[i].data_layout = m_fieldDataLayoutVec[i];
    m_layout_when_static[i].sparsity_layout = m_fieldSparsityLayoutVec[i];

    // For now, handle dynamic by changing maps to dense,
    // and relation to DynamicRelation
    // if(isSparse()) convertLayoutToDense(); //old version with single layout for all MM
    if(m_fieldSparsityLayoutVec[i] == SparsityLayout::SPARSE)
    {
      convertFieldToDense(i);
    }
  }

  //create the dynamic relation
  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidStaticRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }
    StaticVariableRelationType& rel = relStatic(layout);

    RangeSetType* set1 = &relDominantSet(layout);
    RangeSetType* set2 = &relSecondarySet(layout);

    DynamicVariableRelationType relDyn(set1, set2);
    for(int i = 0; i < rel.fromSetSize(); i++)
    {
      auto&& rel_vec = rel[i];
      for(int j = 0; j < rel_vec.size(); j++)
      {
        relDyn.insert(i, rel_vec[j]);
      }
    }

    relDynamic(layout) = relDyn;

    SLIC_ASSERT(relDynamic(layout).isValid());
    SLIC_ASSERT(relDynamic(layout).totalSize() == relStatic(layout).totalSize());
  }

  relStatic(DataLayout::CELL_DOM) = StaticVariableRelationType {};
  relSparseSet(DataLayout::CELL_DOM) = RelationSetType {};
  relStatic(DataLayout::MAT_DOM) = StaticVariableRelationType {};
  relSparseSet(DataLayout::MAT_DOM) = RelationSetType {};

  m_dynamic_mode = true;
}

void MultiMat::convertToStatic()
{
  if(!m_dynamic_mode) return;

  // Change dynamicRelation back to staticRelation
  // change the layout to previously stored static layout

  SLIC_ASSERT(!hasValidStaticRelation(DataLayout::CELL_DOM));
  SLIC_ASSERT(!hasValidStaticRelation(DataLayout::MAT_DOM));

  //Create the static relations
  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidDynamicRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }

    DynamicVariableRelationType& relDyn = relDynamic(layout);

    RangeSetType& set1 = relDominantSet(layout);
    RangeSetType& set2 = relSecondarySet(layout);

    IndBufferType& rel_beginvec = relBeginVec(layout);
    IndBufferType& rel_indicesVec = relIndVec(layout);

    SLIC_ASSERT(SetPosType(rel_beginvec.size()) == set1.size() + 1);
    int rel_data_size = 0;
    for(int i = 0; i < relDyn.fromSetSize(); i++)
    {
      auto& rel_vec = relDyn[i];
      rel_beginvec[i] = rel_data_size;
      rel_data_size += rel_vec.size();
    }
    rel_beginvec.back() = rel_data_size;
    rel_indicesVec.resize(rel_data_size);
    int idx = 0;
    for(int i = 0; i < relDyn.fromSetSize(); i++)
    {
      auto& rel_vec = relDyn[i];
      for(unsigned int j = 0; j < rel_vec.size(); j++)
      {
        rel_indicesVec[idx++] = rel_vec[j];
      }
    }
    SLIC_ASSERT(idx == relDyn.totalSize());

    StaticVariableRelationType& rel = relStatic(layout);

    rel = StaticVariableRelationType(&set1, &set2);
    rel.bindBeginOffsets(set1.size(), &rel_beginvec);
    rel.bindIndices(rel_indicesVec.size(), &rel_indicesVec);

    SLIC_ASSERT(rel.isValid());

    RelationSetType& nzSet = relSparseSet(layout);
    SLIC_ASSERT(nzSet.getRelation() == nullptr);
    nzSet = RelationSetType(&rel);
    SLIC_ASSERT(nzSet.isValid());
  }

  m_dynamic_mode = false;

  //Change each field to their corresponding sparsity
  //if (m_static_layout.sparsity_layout == SparsityLayout::SPARSE) convertLayoutToSparse();
  const int SZ = m_layout_when_static.size();
  for(auto i = 0; i < SZ; i++)
  {
    if(m_layout_when_static[i].sparsity_layout == SparsityLayout::SPARSE)
    {
      convertFieldToSparse(i);
    }
  }
  m_layout_when_static.resize(0);

  // Unset the dynamic relations
  relDynamic(DataLayout::CELL_DOM) = DynamicVariableRelationType {};
  relDynamic(DataLayout::MAT_DOM) = DynamicVariableRelationType {};
}

bool MultiMat::addEntry(int cell_id, int mat_id)
{
  SLIC_ASSERT(m_dynamic_mode);

  //right now this is implemented by converting the data layout to dense so that
  // adding an entry only involve changing the relation data.

  bool searched = false;
  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidDynamicRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }

    DynamicVariableRelationType& relDyn = relDynamic(layout);

    std::pair<int, int> idx;
    if(layout == DataLayout::CELL_DOM)
    {
      idx = {cell_id, mat_id};
    }
    else
    {
      idx = {mat_id, cell_id};
    }

    auto& rel_vec = relDyn.data(idx.first);
    if(!searched)
    {
      // Search through the first valid relation we encounter to check if the
      // entry already exists
      auto&& found_iter = std::find(rel_vec.begin(), rel_vec.end(), idx.second);
      if(found_iter != rel_vec.end())
      {
        SLIC_ASSERT_MSG(false, "MultiMat::addEntry() -- entry already exists.");
        return false;
      }
    }

    relDyn.insert(idx.first, idx.second);
  }

  return true;
}

bool MultiMat::removeEntry(int cell_id, int mat_id)
{
  SLIC_ASSERT(m_dynamic_mode);

  for(DataLayout layout : {DataLayout::CELL_DOM, DataLayout::MAT_DOM})
  {
    if(!hasValidDynamicRelation(layout))
    {
      // We don't have a relation for this layout type
      continue;
    }

    DynamicVariableRelationType& relDyn = relDynamic(layout);

    std::pair<int, int> idx;
    if(layout == DataLayout::CELL_DOM)
    {
      idx = {cell_id, mat_id};
    }
    else
    {
      idx = {mat_id, cell_id};
    }

    auto& rel_vec = relDyn.data(idx.first);
    auto&& found_iter = std::find(rel_vec.begin(), rel_vec.end(), idx.second);
    if(found_iter == rel_vec.end())
    {
      SLIC_ASSERT_MSG(false, "MultiMat::removeEntry() -- entry not found");
      return false;
    }

    rel_vec.erase(found_iter);
  }

  //TODO make all map value zero?

  return true;
}

void MultiMat::makeOtherRelation(DataLayout layout)
{
  DataLayout old_layout =
    (layout == DataLayout::CELL_DOM ? DataLayout::MAT_DOM : DataLayout::CELL_DOM);
  StaticVariableRelationType& oldRel = relStatic(old_layout);
  StaticVariableRelationType& newRel = relStatic(layout);
  IndBufferType& newBeginVec = relBeginVec(layout);
  IndBufferType& newIndicesVec = relIndVec(layout);

  RangeSetType& set1 = *(oldRel.fromSet());
  RangeSetType& set2 = *(oldRel.toSet());

  auto nz_count = oldRel.totalSize();
  //auto nz_count = m_cellMatRel_indicesVec.size();

  newBeginVec.resize(set2.size() + 1, 0);
  newIndicesVec.resize(nz_count, -1);

  //construct the new transposed relation

  //count the non-zero in each rows
  for(auto idx1 = 0; idx1 < oldRel.fromSetSize(); ++idx1)
  {
    IdSet relSubset = oldRel[idx1];
    for(auto j = 0; j < relSubset.size(); ++j)
    {
      auto idx2 = relSubset[j];
      newBeginVec[idx2] += 1;
    }
  }

  //add them to make this the end index
  {
    unsigned int i;
    for(i = 1; i < newBeginVec.size() - 1; i++)
    {
      newBeginVec[i] += newBeginVec[i - 1];
    }
    newBeginVec[i] = newBeginVec[i - 1];
  }

  //fill in the indicesVec and the move_indices backward
  for(auto idx1 = oldRel.fromSetSize() - 1; idx1 >= 0; --idx1)
  {
    IdSet relSubset = oldRel[idx1];
    for(auto j = relSubset.size() - 1; j >= 0; --j)
    {
      auto idx2 = relSubset[j];
      auto compress_idx = --newBeginVec[idx2];
      newIndicesVec[compress_idx] = idx1;
    }
  }

  newRel = StaticVariableRelationType(&set2, &set1);
  newRel.bindBeginOffsets(set2.size(), &newBeginVec);
  newRel.bindIndices(newIndicesVec.size(), &newIndicesVec);

  relSparseSet(layout) = RelationSetType(&newRel);
  relDenseSet(layout) = ProductSetType(&set2, &set1);
}

void MultiMat::convertLayoutToCellDominant()
{
  for(unsigned int i = 0; i < m_mapVec.size(); ++i)
  {
    convertFieldToCellDom(i);
  }
}

void MultiMat::convertLayoutToMaterialDominant()
{
  for(unsigned int i = 0; i < m_mapVec.size(); ++i)
  {
    convertFieldToMatDom(i);
  }
}

void MultiMat::convertLayoutToSparse()
{
  for(unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    convertFieldToSparse(map_i);
  }
}

void MultiMat::convertLayoutToDense()
{
  for(unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    convertFieldToDense(map_i);
  }
}

void MultiMat::convertLayout(DataLayout new_layout, SparsityLayout new_sparsity)
{
  //sparse/dense conversion
  if(new_sparsity == SparsityLayout::SPARSE)
  {
    convertLayoutToSparse();
  }
  else if(new_sparsity == SparsityLayout::DENSE)
  {
    convertLayoutToDense();
  }

  //cell/mat centric conversion
  if(new_layout == DataLayout::MAT_DOM)
  {
    convertLayoutToMaterialDominant();
  }
  else if(new_layout == DataLayout::CELL_DOM)
  {
    convertLayoutToCellDominant();
  }
}

void MultiMat::convertFieldLayout(int field_idx,
                                  SparsityLayout new_sparsity,
                                  DataLayout new_layout)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_mapVec.size()));

  DataLayout field_data_layout = m_fieldDataLayoutVec[field_idx];
  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  if(new_layout == field_data_layout && new_sparsity == field_sparsity_layout)
    return;

  //sparse/dense conversion
  if(field_sparsity_layout == SparsityLayout::DENSE &&
     new_sparsity == SparsityLayout::SPARSE)
  {
    convertFieldToSparse(field_idx);
  }
  else if(field_sparsity_layout == SparsityLayout::SPARSE &&
          new_sparsity == SparsityLayout::DENSE)
  {
    convertFieldToDense(field_idx);
  }

  //cell/mat centric conversion
  if(field_data_layout == DataLayout::CELL_DOM &&
     new_layout == DataLayout::MAT_DOM)
  {
    convertFieldToMatDom(field_idx);
  }
  else if(field_data_layout == DataLayout::MAT_DOM &&
          new_layout == DataLayout::CELL_DOM)
  {
    convertFieldToCellDom(field_idx);
  }

  SLIC_ASSERT(false);
}

void MultiMat::convertFieldToSparse(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_mapVec.size()));

  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  //no conversion needed unless the field is PER_CELL_MAT
  if(field_sparsity_layout == SparsityLayout::SPARSE ||
     m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT)
  {
    return;
  }

  if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeDouble)
  {
    convertToSparse_helper<double>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeFloat)
  {
    convertToSparse_helper<float>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeInt)
  {
    convertToSparse_helper<int>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeUnsignChar)
  {
    convertToSparse_helper<unsigned char>(field_idx);
  }
  else
    SLIC_ASSERT(false);  //TODO

  m_fieldSparsityLayoutVec[field_idx] = SparsityLayout::SPARSE;
}

void MultiMat::convertFieldToDense(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_mapVec.size()));

  SparsityLayout field_sparsity_layout = m_fieldSparsityLayoutVec[field_idx];

  //no conversion needed unless the field is PER_CELL_MAT
  if(field_sparsity_layout == SparsityLayout::DENSE ||
     m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT)
  {
    return;
  }

  if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeDouble)
  {
    convertToDense_helper<double>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeFloat)
  {
    convertToDense_helper<float>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeInt)
  {
    convertToDense_helper<int>(field_idx);
  }
  else if(m_dataTypeVec[field_idx] == DataTypeSupported::TypeUnsignChar)
  {
    convertToDense_helper<unsigned char>(field_idx);
  }
  else
    SLIC_ASSERT(false);  //TODO

  m_fieldSparsityLayoutVec[field_idx] = SparsityLayout::DENSE;
}

template <typename DataType>
void MultiMat::convertToSparse_helper(int map_i)
{
  SLIC_ASSERT(m_fieldSparsityLayoutVec[map_i] != SparsityLayout::SPARSE);

  MapBaseType* mapPtr = m_mapVec[map_i].get();

  //Skip if no volume fraction array is set-up
  if(map_i == 0 && mapPtr == nullptr) return;

  StaticVariableRelationType* Rel = getRel(map_i);

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(Rel->totalSize() * stride);
  int idx = 0;
  for(int i = 0; i < Rel->fromSetSize(); ++i)
  {
    auto relset = (*Rel)[i];
    auto submap = old_map(i);
    for(int j = 0; j < relset.size(); ++j)
    {
      for(int s = 0; s < stride; ++s)
      {
        arr_data[idx++] = submap[relset[j] * stride + s];
      }
    }
  }
  SLIC_ASSERT(idx == Rel->totalSize() * stride);

  RelationSetType* nz_set = &relSparseSet(m_fieldDataLayoutVec[map_i]);
  //old field2d
  //Field2D<DataType>* new_field = new Field2D<DataType>(nz_set, DataType(), stride);
  //new_field->copy(arr_data.data());
  Field2D<DataType>* new_field =
    new Field2D<DataType>(*this, nz_set, old_map.getName(), arr_data.data(), stride);

  m_mapVec[map_i].reset(new_field);
}

template <typename DataType>
void MultiMat::convertToDense_helper(int map_i)
{
  SLIC_ASSERT(m_fieldSparsityLayoutVec[map_i] != SparsityLayout::DENSE);

  MapBaseType* mapPtr = m_mapVec[map_i].get();

  //Skip if no volume fraction array is set-up
  if(map_i == 0 && mapPtr == nullptr) return;

  ProductSetType* prod_set = &relDenseSet(m_fieldDataLayoutVec[map_i]);

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(prod_set->size() * stride);
  for(int i = 0; i < old_map.firstSetSize(); ++i)
  {
    for(auto iter = old_map.begin(i); iter != old_map.end(i); ++iter)
    {
      int elem_idx = i * old_map.secondSetSize() + iter.index();
      for(int c = 0; c < stride; ++c)
      {
        arr_data[elem_idx * stride + c] = iter.value(c);
      }
    }
  }

  //old field2d
  //Field2D<DataType>* new_field = new Field2D<DataType>(prod_set, DataType(), stride);
  //new_field->copy(&arr_data[0]);
  Field2D<DataType>* new_field = new Field2D<DataType>(*this,
                                                       prod_set,
                                                       old_map.getName(),
                                                       arr_data.data(),
                                                       stride);

  m_mapVec[map_i].reset(new_field);
}

DataLayout MultiMat::getFieldDataLayout(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_mapVec.size()));

  return m_fieldDataLayoutVec[field_idx];
}

SparsityLayout MultiMat::getFieldSparsityLayout(int field_idx)
{
  SLIC_ASSERT(0 <= field_idx && field_idx < static_cast<int>(m_mapVec.size()));

  return m_fieldSparsityLayoutVec[field_idx];
}

template <typename DataType>
void MultiMat::transposeField_helper(int field_idx)
{
  SLIC_ASSERT(m_fieldMappingVec[field_idx] == FieldMapping::PER_CELL_MAT);

  MapBaseType* oldMapPtr = m_mapVec[field_idx].get();

  //Skip if no volume fraction array is set-up
  if(field_idx == 0 && oldMapPtr == nullptr) return;

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(oldMapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data;

  DataLayout oldDataLayout = getFieldDataLayout(field_idx);
  StaticVariableRelationType& oldRel = relStatic(oldDataLayout);
  DataLayout new_layout;
  if(oldDataLayout == DataLayout::CELL_DOM)
  {
    new_layout = DataLayout::MAT_DOM;
  }
  else
  {
    new_layout = DataLayout::CELL_DOM;
  }

  if(!hasValidStaticRelation(new_layout))
  {
    makeOtherRelation(new_layout);
  }
  RelationSetType* newNZSet = &relSparseSet(new_layout);
  ProductSetType* newProdSet = &relDenseSet(new_layout);

  auto& set1 = *(oldRel.fromSet());
  auto& set2 = *(oldRel.toSet());

  int set1Size = set1.size();
  int set2Size = set2.size();

  Field2D<DataType>* new_map = nullptr;
  if(m_fieldSparsityLayoutVec[field_idx] == SparsityLayout::SPARSE)
  {
    //copy begin vector for moving
    IndBufferType vec_idx =
      relBeginVec(new_layout);  //a copy of the beginVec to keep track
    const IndBufferType& indicesVec = oldRel.relationData();

    arr_data.resize(oldRel.totalSize() * stride);
    for(int i = 0; i < oldRel.totalSize(); ++i)
    {
      int col = indicesVec[i];
      for(int c = 0; c < stride; ++c)
      {
        arr_data[vec_idx[col] * stride + c] = (*old_map.getMap())(i, c);
      }
      ++vec_idx[col];
    }

    //old
    //new_map = new Field2D<DataType>(newNZSet, DataType(), stride);
    new_map = new Field2D<DataType>(*this,
                                    newNZSet,
                                    old_map.getName(),
                                    arr_data.data(),
                                    stride);
  }
  else  //dense
  {
    arr_data.resize(set1Size * set2Size * stride);
    for(int i = 0; i < set1Size; ++i)
    {
      for(auto iter = old_map.begin(i); iter != old_map.end(i); ++iter)
      {
        int elem_idx = iter.index() * set1Size + i;
        for(int c = 0; c < stride; ++c)
        {
          arr_data[elem_idx * stride + c] = iter.value(c);
        }
      }
    }
    //new_map = new Field2D<DataType>(newProdSet, DataType(), stride);
    new_map = new Field2D<DataType>(*this,
                                    newProdSet,
                                    old_map.getName(),
                                    arr_data.data(),
                                    stride);
  }

  //new_map->copy(arr_data.data());
  m_mapVec[field_idx].reset(new_map);
  m_fieldDataLayoutVec[field_idx] = new_layout;
}

void MultiMat::transposeField(int field_idx)
{
  if(m_fieldMappingVec[field_idx] != FieldMapping::PER_CELL_MAT) return;

  switch(m_dataTypeVec[field_idx])
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
  if(getFieldDataLayout(field_idx) == DataLayout::MAT_DOM)
  {
    return;
  }
  transposeField(field_idx);
}

void MultiMat::convertFieldToCellDom(int field_idx)
{
  if(getFieldDataLayout(field_idx) == DataLayout::CELL_DOM)
  {
    return;
  }
  transposeField(field_idx);
}

std::string MultiMat::getFieldDataLayoutAsString(int field_i) const
{
  if(m_fieldDataLayoutVec[field_i] == DataLayout::CELL_DOM)
    return "Cell-Centric";
  else if(m_fieldDataLayoutVec[field_i] == DataLayout::MAT_DOM)
    return "Material-Centric";
  else
    SLIC_ASSERT(false);
  return "";
}

std::string MultiMat::getFieldSparsityLayoutAsString(int field_i) const
{
  if(m_fieldSparsityLayoutVec[field_i] == SparsityLayout::SPARSE)
    return "Sparse";
  else if(m_fieldSparsityLayoutVec[field_i] == SparsityLayout::DENSE)
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
  sstr << "\nNumber of cells:     " << m_ncells;

  sstr << "\n\n Number of fields: " << m_mapVec.size() << "\n";
  for(unsigned int i = 0; i < m_mapVec.size(); i++)
  {
    sstr << "Field " << i << ": " << m_fieldNameVec[i].c_str();
    sstr << "  Mapping per ";
    switch(m_fieldMappingVec[i])
    {
    case FieldMapping::PER_CELL:
      sstr << "cell";
      break;
    case FieldMapping::PER_MAT:
      sstr << "material";
      break;
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

  if(getCellSet().size() > 0 && getMatSet().size() > 0)
  {
    //It's a non-empty MultiMat object.
    //Check Volfrac exists
    if(m_mapVec[0] == nullptr)
    {
      errStr << "\n\t*No Volfrac field added.";
      bValid = false;
    }
    else
    {
      const auto& volfrac_map =
        *dynamic_cast<Field2D<double>*>(m_mapVec[0].get());
      auto volfrac_sparsity = m_fieldSparsityLayoutVec[0];
      auto volfrac_layout = m_fieldDataLayoutVec[0];

      //Check Volfrac values match the relation value (if dense)
      if(volfrac_sparsity == SparsityLayout::DENSE && !m_dynamic_mode)
      {
        const StaticVariableRelationType& relPtr =
          m_staticRelations[(int)volfrac_layout];

        for(int i = 0; i < volfrac_map.firstSetSize(); ++i)
        {
          auto rel_iter = relPtr.begin(i);
          for(int j = 0; j < volfrac_map.secondSetSize(); ++j)
          {
            double volfrac = volfrac_map[i * volfrac_map.secondSetSize() + j];
            bool zero_volfrac = volfrac == 0.0;  //exact comp?
            if(rel_iter != relPtr.end(i) && *rel_iter == j)
            {  //cellmat rel is present
              if(zero_volfrac)
              {
                errStr << "\n\t*Volume fraction is zero for a material "
                       << " that exists in a cell";
                bValid = false;
              }
              ++rel_iter;
            }
            else
            {
              if(!zero_volfrac)
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
      std::vector<double> volfrac_sum(getCellSet().size(), 0.0);
      for(int i = 0; i < volfrac_map.firstSetSize(); ++i)
      {
        const auto& submap = volfrac_map(i);
        for(int j = 0; j < submap.size(); ++j)
        {
          if(volfrac_layout == DataLayout::CELL_DOM)
            volfrac_sum[i] += submap.value(j);
          else
            volfrac_sum[submap.index(j)] += submap.value(j);
        }
      }
      for(unsigned int i = 0; i < volfrac_sum.size(); ++i)
      {
        if(abs(volfrac_sum[i] - 1.0) > 10e-9)
        {
          errStr << "\n\t*Volfrac does not sum to 1.0 in cell " << i;
          bValid = false;
        }
      }
    }
  }

  if(verboseOutput)
  {
    if(bValid)
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
  switch(fm)
  {
  case FieldMapping::PER_CELL:
    set_ptr = &getCellSet();
    break;
  case FieldMapping::PER_MAT:
    set_ptr = &getMatSet();
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

MultiMat::SetType* MultiMat::get_mapped_set(int field_idx)
{
  SetType* set_ptr = nullptr;
  switch(m_fieldMappingVec[field_idx])
  {
  case FieldMapping::PER_CELL:
    set_ptr = &getCellSet();
    break;
  case FieldMapping::PER_MAT:
    set_ptr = &getMatSet();
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

MultiMat::BivariateSetType* MultiMat::get_mapped_biSet(DataLayout layout,
                                                       SparsityLayout sparsity)
{
  BivariateSetType* set_ptr = nullptr;

  if(sparsity == SparsityLayout::SPARSE)
  {
    set_ptr = &relSparseSet(layout);
  }
  else if(sparsity == SparsityLayout::DENSE)
  {
    set_ptr = &relDenseSet(layout);
  }

  return set_ptr;
}

MultiMat::BivariateSetType* MultiMat::get_mapped_biSet(int field_idx)
{
  SLIC_ASSERT(m_fieldMappingVec[field_idx] == FieldMapping::PER_CELL_MAT);
  SLIC_ASSERT(0 <= field_idx &&
              field_idx < static_cast<int>(m_fieldNameVec.size()));

  DataLayout layout = m_fieldDataLayoutVec[field_idx];
  SparsityLayout sparsity = m_fieldSparsityLayoutVec[field_idx];

  return get_mapped_biSet(layout, sparsity);
}

MultiMat::StaticVariableRelationType* MultiMat::getRel(int field_idx)
{
  SLIC_ASSERT(m_fieldDataLayoutVec[field_idx] == DataLayout::CELL_DOM ||
              m_fieldDataLayoutVec[field_idx] == DataLayout::MAT_DOM);

  //get the sparse relation vector
  DataLayout layout = m_fieldDataLayoutVec[field_idx];
  StaticVariableRelationType* rel_ptr = &relStatic(layout);

  return rel_ptr;
}
