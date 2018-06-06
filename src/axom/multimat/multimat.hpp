#ifndef MULTIMAT_H_
#define MULTIMAT_H_


#include "slam/RangeSet.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/StaticRelation.hpp"
#include "slam/Map.hpp"
#include "slam/MappedRelationSet.hpp"
#include "slam/RelationMap.hpp"
#include <vector>
#include <cassert>

namespace axom
{
namespace multimat
{

namespace policies = slam::policies;


enum FieldMapping { PER_CELL, PER_MAT, PER_CELL_MAT };
enum DataLayout { LAYOUT_CELL_DOM, LAYOUT_MAT_DOM };
enum SparcityLayout { LAYOUT_SPARSE, LAYOUT_DENSE };

class MultiMat
{

  /**  private type def **/
private: 
  // SLAM Set type definitions
  using SetType = slam::RangeSet;
  using RangeSetType   = slam::RangeSet;
  using SetPosType  = SetType::PositionType;
  using SetElemType = SetType::ElementType;
  // SLAM Relation typedef (static variable-cardinality relation)
  using STLIndirection             = policies::STLVectorIndirection<SetPosType, SetPosType>;
  using VariableCardinality = policies::VariableCardinality<SetPosType, STLIndirection>;
  using StaticVariableRelationType = slam::StaticRelation<VariableCardinality, STLIndirection,
    RangeSetType, RangeSetType>;
  using RelationSet = StaticVariableRelationType::RelationSet; 
  // SLAM MappedRelationSet for the set of non-zero cell to mat variables
  using MappedRelationSetType = slam::MappedRelationSet<StaticVariableRelationType>;
  
  // SLAM Map type
  using MapBaseType = slam::MapBase;
  template <typename T>
  using MapType = slam::Map<T>;
  template <typename T>
  using RelationMapType = slam::RelationMap<T, StaticVariableRelationType>;
  template<typename T>
  using SubsetMap = slam::SubsetMap<T>;

  /** Public type def **/
public:
  template <typename T>
  using Field1D = MapType<T>;
  template <typename T>
  using Field2D = RelationMapType<T>;
  template <typename T>
  using SubField = SubsetMap<T>;
  using IndexSet = RangeSetType;
  using IdSet = RelationSet;
  
  /** Public functions **/
  MultiMat();
  MultiMat(DataLayout, SparcityLayout);
  
  //Set-up functions
  void setNumberOfMat(int); 
  void setNumberOfCell(int);
  void setCellMatRel(std::vector<bool>&); //relation information

  //functions related to the field arrays
  template<class T>
  int newFieldArray(std::string arr_name, FieldMapping, T* arr);
 
  template<typename T>
  Field1D<T>& get1dField(std::string field_name);
  template<typename T>
  Field2D<T>& get2dField(std::string field_name);

  IdSet getMatInCell(int c); //Should change the func name so there's no assumption of the layout
  IndexSet getIndexingSetOfCell(int c);

  int getNumberOfMaterials() { return m_nmats; };
  int getNumberOfCells() { return m_ncells; }

  //Data modification functions
  //...

  //Layout modification functions
  void convertLayoutToCellDominant() { convertLayout(LAYOUT_CELL_DOM, m_sparcityLayout); }
  void convertLayoutToMaterialDominant() { convertLayout(LAYOUT_MAT_DOM, m_sparcityLayout); }
  void convertLayout(DataLayout, SparcityLayout);
  DataLayout getDataLayout();
  std::string getLayoutAsString();

  void printSelf();
  bool isValid();
  
private: //private functions
  SetType* get_mapped_set(FieldMapping fm);


private:
  int m_nmats, m_ncells;
  DataLayout m_dataLayout;
  SparcityLayout m_sparcityLayout;

  //slam variables
  SetType m_cellSet;
  SetType m_matSet;
  StaticVariableRelationType m_cell2matRel;
  std::vector<SetPosType> m_cell2matRel_beginsVec; //to store the cell2mat relation
  std::vector<SetPosType> m_cell2matRel_indicesVec; //to store the cell2mat relation
  MappedRelationSetType m_cellMatNZSet; // set of non-zero entries in the cellXmat matrix
  
  std::vector<std::string> m_arrNameVec;
  std::vector<FieldMapping> m_fieldMappingVec;
  std::vector<MapBaseType*> m_mapVec; 

  friend class MultiMatArray;
  
}; //end MultiMat class


//--------------- MultiMat template function definitions -----------------//


template<class T>
int MultiMat::newFieldArray(std::string arr_name, FieldMapping arr_mapping, T* data_arr)
{
  int index_val = m_mapVec.size();
  
  if (arr_mapping == PER_CELL_MAT) {
    MappedRelationSetType* s = &m_cellMatNZSet;
    
    RelationMapType<T>* newMapPtr = new RelationMapType<T>(s, data_arr);

    m_mapVec.push_back(newMapPtr);
  }
  else if (arr_mapping == PER_CELL || arr_mapping == PER_MAT)
  {
    MapType<T>* newMapPtr = new MapType<T>(get_mapped_set(arr_mapping));
    
    int i = 0;
    for (auto iter = newMapPtr->begin(); iter != newMapPtr->end(); iter++) {
      *iter = data_arr[i++];
    }

    m_mapVec.push_back(newMapPtr);
  }
  m_arrNameVec.push_back(arr_name);
  m_fieldMappingVec.push_back(arr_mapping);
  assert(m_arrNameVec.size() == m_mapVec.size() && m_mapVec.size() == m_fieldMappingVec.size());

  return index_val;
}

template<typename T>
MultiMat::Field1D<T>& MultiMat::get1dField(std::string field_name)
{
  for (int i = 0; i < m_arrNameVec.size(); i++)
  {
    if (m_arrNameVec[i] == field_name)
    {
      //assert(m_fieldMappingVec[i] == PER_CELL || m_fieldMappingVec[i] == PER_MAT);
      return *dynamic_cast<Field1D<T>*>(m_mapVec[i]);
    }
  }
  assert(false); //No array with such name found
}

template<typename T>
inline MultiMat::RelationMapType<T>& MultiMat::get2dField(std::string field_name)
{
  for (int i = 0; i < m_arrNameVec.size(); i++)
  {
    if (m_arrNameVec[i] == field_name)
    {
      assert(m_fieldMappingVec[i] == PER_CELL_MAT);
      return *dynamic_cast<RelationMapType<T>*>(m_mapVec[i]);
    }
  }
}

} //end namespace multimat
} //end namespace axom


#endif