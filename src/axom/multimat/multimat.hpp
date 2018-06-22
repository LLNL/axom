#ifndef MULTIMAT_H_
#define MULTIMAT_H_


#include "slam/RangeSet.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/StaticRelation.hpp"
#include "slam/Map.hpp"
#include "slam/MappedRelationSet.hpp"
#include "slam/ProductSet.hpp"
#include "slam/BivariateMap.hpp"
#include <vector>
#include <cassert>

namespace axom
{
namespace multimat
{

namespace policies = slam::policies;


enum class FieldMapping { PER_CELL, PER_MAT, PER_CELL_MAT };
enum class DataLayout { CELL_CENTRIC, MAT_CENTRIC };
enum class SparcityLayout { SPARSE, DENSE };
enum class DataTypeSupported { TypeUnknown, TypeInt, TypeDouble, TypeFloat, TypeUnsignChar };

class MultiMat
{

  /**  private type def **/
private: 
  // SLAM Set type definitions
  using SetType = slam::Set;
  using RangeSetType   = slam::RangeSet;
  using BivariateSetType = slam::BivariateSet;
  using SetPosType  = SetType::PositionType;
  using SetElemType = SetType::ElementType;
  // SLAM Relation typedef (static variable-cardinality relation)
  using STLIndirection             = policies::STLVectorIndirection<SetPosType, SetPosType>;
  using VariableCardinality = policies::VariableCardinality<SetPosType, STLIndirection>;
  using StaticVariableRelationType = slam::StaticRelation<VariableCardinality, STLIndirection,
    RangeSetType, RangeSetType>;
  using OrderedSetType = slam::OrderedSet<
                policies::RuntimeSize<SetType::PositionType>,
                policies::RuntimeOffset<SetType::PositionType>,
                policies::StrideOne<SetType::PositionType>,
                policies::STLVectorIndirection<SetType::PositionType, SetType::ElementType>,
                policies::NoSubset>;
  using ProductSetType = slam::ProductSet;
  // SLAM MappedRelationSet for the set of non-zero cell to mat variables
  using MappedRelationSetType = slam::MappedRelationSet<StaticVariableRelationType>;
  
  //stride-related
  using StrideType = slam::policies::RuntimeStride<SetPosType>;

  // SLAM Map type
  using MapBaseType = slam::MapBase;
  template <typename T>
  using MapType = slam::Map<T, StrideType>;
  template <typename T>
  using BivariateMapType = slam::BivariateMap<T, StrideType>;

  template<typename T>
  using SubMap = slam::SubMap<T, StrideType>;

/** Public type def **/
public:
  template <typename T>
  using Field1D = MapType<T>;
  template <typename T>
  using Field2D = BivariateMapType<T>;
  template <typename T>
  using SubField = SubMap<T>;
  using IndexSet = RangeSetType;
  using IdSet = OrderedSetType;
  
/** Public functions **/
  MultiMat();
  MultiMat(DataLayout, SparcityLayout);
  ~MultiMat();

  //Set-up functions
  void setNumberOfMat(int); 
  void setNumberOfCell(int);
  void setCellMatRel(std::vector<bool>&); //relation information

  //functions related to the field arrays
  template<class T>
  int newFieldArray(std::string arr_name, FieldMapping, T* arr, int stride = 1);
 
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
  void convertLayoutToCellDominant() { convertLayout(DataLayout::CELL_CENTRIC, m_sparcityLayout); }
  void convertLayoutToMaterialDominant() { convertLayout(DataLayout::MAT_CENTRIC, m_sparcityLayout); }
  void convertLayoutToSparse() { convertLayout(m_dataLayout, SparcityLayout::SPARSE); }
  void convertLayoutToDense() { convertLayout(m_dataLayout, SparcityLayout::DENSE); }
  void convertLayout(DataLayout, SparcityLayout);
  DataLayout getDataLayout();
  std::string getLayoutAsString();
  SparcityLayout getSparcityLayout();

  void printSelf() const;
  bool isValid(bool verboseOutput = false) const;
  
private: //private functions
  SetType* get_mapped_set(FieldMapping fm);
  template<typename DataType>
  void convertToSparse_helper(int map_i);

private:
  unsigned int m_nmats, m_ncells;
  DataLayout m_dataLayout;
  SparcityLayout m_sparcityLayout;

  //slam variables
  RangeSetType m_cellSet;
  RangeSetType m_matSet;
  StaticVariableRelationType m_cellMatRel;
  std::vector<SetPosType> m_cellMatRel_beginsVec; //to store the cell2mat relation
  std::vector<SetPosType> m_cellMatRel_indicesVec; //to store the cell2mat relation
  MappedRelationSetType m_cellMatNZSet; // set of non-zero entries in the cellXmat matrix
  ProductSetType m_cellMatProdSet;

  std::vector<std::string> m_arrNameVec;
  std::vector<FieldMapping> m_fieldMappingVec;
  std::vector<MapBaseType*> m_mapVec; 
  std::vector<DataTypeSupported> m_dataTypeVec;

}; //end MultiMat class


//--------------- MultiMat template function definitions -----------------//


template<class T>
int MultiMat::newFieldArray(std::string arr_name, FieldMapping arr_mapping, T* data_arr, int stride)
{
  assert(stride > 0);

  int index_val = m_mapVec.size();
  
  if (arr_mapping == FieldMapping::PER_CELL_MAT) 
  {
    BivariateSetType* s = nullptr;
    if (m_sparcityLayout == SparcityLayout::SPARSE) {
      s = &m_cellMatNZSet;
    }
    else if (m_sparcityLayout == SparcityLayout::DENSE) {
      s = &m_cellMatProdSet;
    }
    else assert(false);

    Field2D<T>* new_map_ptr = new Field2D<T>(s, data_arr, stride);
    m_mapVec.push_back(new_map_ptr);
  }
  else if (arr_mapping == FieldMapping::PER_CELL || arr_mapping == FieldMapping::PER_MAT)
  {
    Field1D<T>* new_map_ptr = new Field1D<T>(get_mapped_set(arr_mapping), T(), stride);
    
    int i = 0;
    for (auto iter = new_map_ptr->begin(); iter != new_map_ptr->end(); iter++) {
      for (auto s = 0; s < stride; ++s) {
        iter(s) = data_arr[i++];
      }
    }

    m_mapVec.push_back(new_map_ptr);
  }
  else assert(false);

  m_arrNameVec.push_back(arr_name);
  m_fieldMappingVec.push_back(arr_mapping);
  assert(m_arrNameVec.size() == m_mapVec.size() && m_mapVec.size() == m_fieldMappingVec.size());

  if (std::is_same<T, int>::value)
    m_dataTypeVec.push_back(DataTypeSupported::TypeInt);
  else if (std::is_same<T, double>::value)
    m_dataTypeVec.push_back(DataTypeSupported::TypeDouble);
  else if (std::is_same<T, float>::value)
    m_dataTypeVec.push_back(DataTypeSupported::TypeFloat);
  else if (std::is_same<T, unsigned char>::value)
    m_dataTypeVec.push_back(DataTypeSupported::TypeUnsignChar);
  else
    m_dataTypeVec.push_back(DataTypeSupported::TypeUnknown);

  return index_val;
}

template<typename T>
MultiMat::Field1D<T>& MultiMat::get1dField(std::string field_name)
{
  for (unsigned int i = 0; i < m_arrNameVec.size(); i++)
  {
    if (m_arrNameVec[i] == field_name)
    {
      //assert(m_fieldMappingVec[i] == FieldMapping::PER_CELL || m_fieldMappingVec[i] == FieldMapping::PER_MAT);
      if (m_fieldMappingVec[i] == FieldMapping::PER_CELL || m_fieldMappingVec[i] == FieldMapping::PER_MAT)
        return *dynamic_cast<Field1D<T>*>(m_mapVec[i]);
      else if (m_fieldMappingVec[i] == FieldMapping::PER_CELL_MAT) {
        Field2D<T> * map_2d = dynamic_cast<Field2D<T>*>(m_mapVec[i]);
        return *(map_2d->getMap());
      }
    }
  }
  assert(false); //No array with such name found
}

template<typename T>
inline MultiMat::Field2D<T>& MultiMat::get2dField(std::string field_name)
{
  for (unsigned int i = 0; i < m_arrNameVec.size(); i++)
  {
    if (m_arrNameVec[i] == field_name)
    {
      assert(m_fieldMappingVec[i] == FieldMapping::PER_CELL_MAT);
      return *dynamic_cast<Field2D<T>*>(m_mapVec[i]);
    }
  }
}


template<typename DataType>
void MultiMat::convertToSparse_helper(int map_i)
{
  Field2D<DataType>& old_ptr = *dynamic_cast<Field2D<DataType>*>(m_mapVec[map_i]);
  int stride = old_ptr.stride();
  std::vector<DataType> arr_data(m_cellMatNZSet.totalSize()*stride);
  int idx = 0;
  for (int i = 0; i < m_cellMatRel.fromSetSize(); ++i) {
    auto relset = m_cellMatRel[i];
    auto submap = old_ptr(i);
    for (int j = 0; j < relset.size(); ++j) {
      for (int s = 0; s < stride; ++s) {
        arr_data[idx++] = submap[relset[j]*stride + s];
      }
    }
  }
  assert(idx == m_cellMatNZSet.totalSize()*stride);

  Field2D<DataType>* new_field = new Field2D<DataType>(&m_cellMatNZSet, &arr_data[0], stride);
  delete m_mapVec[map_i];
  m_mapVec[map_i] = new_field;
}

} //end namespace multimat
} //end namespace axom


#endif
