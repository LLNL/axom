#ifndef MULTIMAT_H_
#define MULTIMAT_H_


#include "slam/RangeSet.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/StaticRelation.hpp"
#include "slam/Map.hpp"
#include "slam/RelationSet.hpp"
#include "slam/ProductSet.hpp"
#include "slam/BivariateMap.hpp"
#include <vector>
#include <cassert>
#include <stdexcept>

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
  using RelationSetType = slam::RelationSet<StaticVariableRelationType>;
  
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
  //Constructors 
  MultiMat();
  MultiMat(DataLayout, SparcityLayout);
  ~MultiMat();
  MultiMat(const MultiMat&);
  MultiMat& operator=(const MultiMat&);

  //Set-up functions
  void setNumberOfMat(int); 
  void setNumberOfCell(int);
  void setCellMatRel(std::vector<bool>&); //relation information

  //functions related to the field arrays
  template<class T>
  int addField(const std::string& arr_name, FieldMapping, T* arr, int stride = 1);
  int setVolfracField(double* arr); //TODO field

  int getFieldIdx(const std::string& arr_name) const;
  template<typename T>
  Field1D<T>& get1dField(const std::string& field_name);
  template<typename T>
  Field2D<T>& get2dField(const std::string& field_name);
  Field2D<double>& getVolfracField();

  IdSet getMatInCell(int c); //Should change the func name so there's no assumption of the layout
  IndexSet getIndexingSetOfCell(int c);

  int getNumberOfMaterials() const { return m_nmats; };
  int getNumberOfCells() const { return m_ncells; }

  //Data modification functions
  //...

  //Layout modification functions
  void convertLayoutToCellDominant();
  void convertLayoutToMaterialDominant();
  void transposeData();
  void convertLayoutToSparse();
  void convertLayoutToDense();
  void convertLayout(DataLayout, SparcityLayout);
  DataLayout getDataLayout() const { return m_dataLayout; }
  SparcityLayout getSparcityLayout() const { return m_sparcityLayout; }
  std::string getDataLayoutAsString() const ;
  std::string getSparcityLayoutAsString() const;
  bool isSparse() const { return getSparcityLayout() == SparcityLayout::SPARSE; }
  bool isDense() const { return getSparcityLayout() == SparcityLayout::DENSE; }
  bool isCellDom() const { return getDataLayout() == DataLayout::CELL_CENTRIC; }
  bool isMatDom() const { return getDataLayout() == DataLayout::MAT_CENTRIC; }
  FieldMapping getFieldMapping(int field_i) const { return m_fieldMappingVec[field_i]; }

  void print() const;
  bool isValid(bool verboseOutput = false) const;
  
private: //private functions
  SetType* get_mapped_set(FieldMapping fm);
  BivariateSetType* get_mapped_biSet();
  template<typename DataType>
  void convertToSparse_helper(int map_i);
  template<typename DataType>
  void convertToDense_helper(int map_i);
  template<typename DataType>
  void transposeData_helper(int map_i, RelationSetType*, ProductSetType*, 
    std::vector<SetPosType>&);
  template<typename T>
  int addFieldArray_impl(const std::string&, FieldMapping, T*, int);
  template<typename T>
  MapBaseType* helperfun_copyField(const MultiMat&, int map_i);

private:
  unsigned int m_ncells, m_nmats;
  DataLayout m_dataLayout;
  SparcityLayout m_sparcityLayout;

  //slam variables
  RangeSetType m_cellSet;
  RangeSetType m_matSet;

  std::vector<SetPosType> m_cellMatRel_beginsVec; //to store the cell2mat relation
  std::vector<SetPosType> m_cellMatRel_indicesVec;
  StaticVariableRelationType* m_cellMatRel;
  RelationSetType* m_cellMatNZSet; // set of non-zero entries in the cellXmat matrix
  ProductSetType* m_cellMatProdSet;
  
  std::vector<std::string> m_arrNameVec;
  std::vector<FieldMapping> m_fieldMappingVec;
  std::vector<MapBaseType*> m_mapVec; 
  std::vector<DataTypeSupported> m_dataTypeVec;

}; //end MultiMat class



//--------------- MultiMat template function definitions -----------------//

template<class T>
int MultiMat::addFieldArray_impl(const std::string& arr_name, FieldMapping arr_mapping,
  T* data_arr, int stride)
{
  int new_arr_idx = m_mapVec.size();

  if (arr_mapping == FieldMapping::PER_CELL_MAT)
  {
    BivariateSetType* s = get_mapped_biSet();
    assert(s != nullptr);

    Field2D<T>* new_map_ptr = new Field2D<T>(s, T(), stride);
    new_map_ptr->copy(data_arr);
    m_mapVec.push_back(new_map_ptr);
  }
  else if (arr_mapping == FieldMapping::PER_CELL ||
    arr_mapping == FieldMapping::PER_MAT)
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
  assert(m_mapVec.size() == m_arrNameVec.size());
  assert(m_mapVec.size() == m_fieldMappingVec.size());

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

  return new_arr_idx;
}

template<class T>
int MultiMat::addField(const std::string& arr_name, FieldMapping arr_mapping,
                            T* data_arr, int stride)
{
  assert(stride > 0);

  //make sure the name does not conflict
  int fieldIdx = getFieldIdx(arr_name);
  if (fieldIdx == 0 && m_mapVec[0] != nullptr) 
  { //this is the vol frac array. call setVolfrac instead
    assert(arr_mapping == FieldMapping::PER_CELL_MAT);
    assert(stride == 1);
    assert(data_arr != nullptr);
    setVolfracField(data_arr);
    return 0;
  }
  else if(fieldIdx > 0) 
  {
    // There is already an array with the current name. And it's not Volfrac.
    // Don't add the new field
    return -1;
  }
  else 
  {
    //No field with this name. Proceed with adding the field
    return addFieldArray_impl<>(arr_name, arr_mapping, data_arr, stride);
  }
}

template<typename T>
MultiMat::Field1D<T>& MultiMat::get1dField(const std::string& field_name)
{
  int fieldIdx = getFieldIdx(field_name);
  if (fieldIdx >= 0)
  {
    //assert(m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL ||
    //       m_fieldMappingVec[fieldIdx] == FieldMapping::PER_MAT);
    if (m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL || 
        m_fieldMappingVec[fieldIdx] == FieldMapping::PER_MAT)
    {
      return *dynamic_cast<Field1D<T>*>(m_mapVec[fieldIdx]);
    }
    else if (m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL_MAT) 
    {
      Field2D<T>* map_2d = dynamic_cast<Field2D<T>*>(m_mapVec[fieldIdx]);
      return *(map_2d->getMap());
    }
  }
  else
    assert(false); //No array with such name found
}


template<typename T>
MultiMat::Field2D<T>& MultiMat::get2dField(const std::string& field_name)
{
  for (unsigned int i = 0; i < m_arrNameVec.size(); i++)
  {
    if (m_arrNameVec[i] == field_name)
    {
      assert(m_fieldMappingVec[i] == FieldMapping::PER_CELL_MAT);
      return *dynamic_cast<Field2D<T>*>(m_mapVec[i]);
    }
  }

  throw std::invalid_argument("Field with name not found");
}


template<typename DataType>
void MultiMat::convertToSparse_helper(int map_i)
{
  SLIC_ASSERT(m_sparcityLayout != SparcityLayout::SPARSE);
  
  MapBaseType* mapPtr = m_mapVec[map_i];
  if (map_i == 0 && mapPtr == nullptr) { return; }

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(m_cellMatRel->totalSize()*stride);
  int idx = 0;
  for (int i = 0; i < m_cellMatRel->fromSetSize(); ++i) {
    auto relset = (*m_cellMatRel)[i];
    auto submap = old_map(i);
    for (int j = 0; j < relset.size(); ++j) {
      for (int s = 0; s < stride; ++s) {
        arr_data[idx++] = submap[relset[j]*stride + s];
      }
    }
  }
  assert(idx == m_cellMatRel->totalSize()*stride);
  Field2D<DataType>* new_field = new Field2D<DataType>(m_cellMatNZSet, DataType(), stride);
  new_field->copy(&arr_data[0]);

  delete m_mapVec[map_i];
  m_mapVec[map_i] = new_field;
}


template<typename DataType>
void MultiMat::convertToDense_helper(int map_i)
{
  SLIC_ASSERT(m_sparcityLayout != SparcityLayout::DENSE);

  MapBaseType* mapPtr = m_mapVec[map_i];
  if (map_i == 0 && mapPtr == nullptr) { return; }

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(m_cellMatProdSet->size()*stride);
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

  Field2D<DataType>* new_field = new Field2D<DataType>(m_cellMatProdSet, DataType(), stride);
  new_field->copy(&arr_data[0]);

  delete m_mapVec[map_i];
  m_mapVec[map_i] = new_field;
}

template<typename DataType>
void MultiMat::transposeData_helper(int map_i, RelationSetType* new_cellMatNZSet,
ProductSetType* new_cellMatProdSet, std::vector<SetPosType>& move_indices)
{
  MapBaseType* mapPtr = m_mapVec[map_i];
  if (map_i == 0 && mapPtr == nullptr) { return; }
  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data;

  RangeSetType& set1 = *(m_cellMatRel->fromSet());
  RangeSetType& set2 = *(m_cellMatRel->toSet());

  int set1Size = set1.size();
  int set2Size = set2.size();

  Field2D<DataType>* new_map = nullptr;
  if (m_sparcityLayout == SparcityLayout::SPARSE) {
    arr_data.resize(m_cellMatRel->totalSize()*stride);
    for (unsigned int i = 0; i < move_indices.size(); ++i)
    {
      for (int c = 0; c < stride; ++c)
      {
        arr_data[move_indices[i] * stride + c] = (*old_map.getMap())(i, c);
      }
    }
    new_map = new Field2D<DataType>(new_cellMatNZSet, DataType(), stride);
  }
  else //dense
  {
    arr_data.resize(set1Size*set2Size*stride);
    for (int i = 0; i < set1Size; ++i)
    {
      for (auto iter = old_map.begin(i); iter != old_map.end(i); ++iter)
      {
        int elem_idx = iter.index() * set1Size + i;
        for (int c = 0; c < stride; ++c)
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



} //end namespace multimat
} //end namespace axom


#endif
