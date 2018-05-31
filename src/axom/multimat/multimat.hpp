#ifndef MULTIMAT_H_
#define MULTIMAT_H_


#include "slam/RangeSet.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/StaticRelation.hpp"
#include "slam/Map.hpp"
#include "slam/MappedRelationSet.hpp"
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

//forward declaration of classes
class MultiMatArray;
template<class T> class MultiMatTypedArray;


class MultiMat
{
public:
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
  //using IntMapType = MapType<int>;
  //using doubleMapType = MapType<double>;

  
  MultiMat();
  MultiMat(DataLayout, SparcityLayout);
  
  //Set-up functions
  void setNumberOfMat(int); 
  void setNumberOfCell(int);
  void setCellMatRel(std::vector<bool>&); //relation information

  RelationSet getMatInCell(int c); //Should change it so there's no assumption of the layout

  //functions related to the field arrays
  template<class T>
  MultiMatTypedArray<T>* newFieldArray(std::string arr_name, FieldMapping, T* arr); //get the array
 
  MultiMatArray* getFieldArray(std::string arr_name); //get the array
  template<typename T>
  MultiMatTypedArray<T>* getFieldArray(std::string arr_name); //get the array
  MultiMatArray* getFieldArray(int arr_idx); //get the array
  
  //accessing functions
  template<class T>
  T getFieldValue(std::string field_name, int cell_i, int mat_i); //can be slow
  template<class T>
  T getFieldValue(std::string field_name, int i);

  template<typename T>
  const MapType<T>& getMap(std::string field_name);

  SetType* get_mapped_set(FieldMapping fm);

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
  
public: //private:
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

  std::vector<MultiMatArray*> m_fieldArrayVec; //list of MM-Arrays

  friend class MultiMatArray;
  
}; //end MultiMat class


class MultiMatArray
{
public:
  using SetType = MultiMat::SetType;
  using RangeSetType = MultiMat::RangeSetType;
  using MapBaseType = axom::slam::MapBase;
  template<typename T>
  using MapType = axom::slam::Map<T>;

  MultiMatArray() {}
  MultiMatArray(MultiMat* m, std::string name, FieldMapping f);

  //MultiMatAbstractArray(MultiMat* m, std::string name, SetType* s, FieldMapping f);
  template<typename T>
  MultiMatArray(MultiMat* m, std::string name, FieldMapping f, T* data_arr);

  ~MultiMatArray() {};
  
  void setName(std::string arr_name) { m_arrayName = arr_name; }
  std::string getName() { return m_arrayName; }

  FieldMapping getFieldMapping() { return m_fieldMapping; };

  template<typename T> T getValue(int c);
  template<typename T> T getValue(int c, int m);

  template<typename T> void setValue(int c, T val);
  template<typename T> void setValue(int c, int m, T val);

  const MapBaseType* getMap() { return m_mapPtr; }
  template<typename T>
  const MapType<T>& getMap() { return *dynamic_cast<MapType<T>*>(m_mapPtr); }
  template<typename T> const MapType<T>& getSubsetMap(int c);


  template<typename T>
  class iterator : public std::iterator<std::input_iterator_tag, T>
  {
    MapType<T>* m_mapPtr;
    int m_i;
  public:
    iterator(MapType<T>* m) : m_mapPtr(m), m_i(0) {}
    iterator(MapType<T>* m, int n) : m_mapPtr(m), m_i(n) {}
    iterator& operator++() { ++m_i; return *this; }
    iterator operator++(int) { iterator tmp(*this); operator++(); return tmp; }
    bool operator==(const iterator& rhs) { return m_mapPtr == rhs.m_mapPtr && m_i == rhs.m_i; }
    bool operator!=(const iterator& rhs) { return m_mapPtr != rhs.m_mapPtr || m_i != rhs.m_i; }
    T operator*() { return m_mapPtr->operator[](m_i); }
  };
  template<typename T>
  iterator<T> begin() { return iterator<T>(dynamic_cast<MapType<T>*>(m_mapPtr)); }
  template<typename T>
  iterator<T> end() { return iterator<T>(dynamic_cast<MapType<T>*>(m_mapPtr), m_mapPtr->size()); }
  

protected:
  FieldMapping m_fieldMapping;
  std::string m_arrayName;

  MapBaseType* m_mapPtr;
  
  MapBaseType* m_subsetmapPtr; //placeholder for a MapBuilder class
  MultiMat::RelationSet m_retSubset; // a temporary set for returning subset map.

  //SetType* m_set; //storing the set in order to index into the elemXmat map. 
                  //Eventually the map should be changed to allow access for this
  MultiMat* mm;

}; //end MultiMatArray class




template<class T>
class MultiMatTypedArray : public MultiMatArray
{
public:
  using SetType = MultiMat::SetType;
  using MapType = axom::slam::Map<T>;

  //MultiMatTypedArray(MultiMat* m, std::string name, SetType* s, FieldMapping f);
  MultiMatTypedArray(MultiMat* m, std::string name, SetType* s, FieldMapping f, T* arr);
  ~MultiMatTypedArray() {}

  T getValue(int c);
  T getValue(int c, int m);

  void setValue(int c, T val);
  void setValue(int c, int m, T val);

  const MapType& getMap() { return m_map; }
  const MapType& getSubsetMap(int c);

  RangeSetType getIndexingSetOfCell(int c);

  class iterator : public std::iterator<std::input_iterator_tag, T>
  {
    MapType* m_mapPtr;
    int m_i;

  public:
    iterator(MapType* m) :m_mapPtr(m), m_i(0) {}
    iterator(MapType* m, int n) : m_mapPtr(m), m_i(n) {}
    iterator& operator++() { ++m_i; return *this; }
    iterator operator++(int) { iterator tmp(*this); operator++(); return tmp; }
    bool operator==(const iterator& rhs) { return m_mapPtr == rhs.m_mapPtr && m_i == rhs.m_i; }
    bool operator!=(const iterator& rhs) { return m_mapPtr != rhs.m_mapPtr || m_i != rhs.m_i; }
    T operator*() { return m_mapPtr->operator[](m_i); }
  };
  iterator begin() { return iterator(&m_map); }
  iterator end() { return iterator(&m_map, m_map.size()); }

protected:
  MapType m_map;
  MapType m_subsetmap; //placeholder for a MapBuilder class
  //MultiMat::RelationSet m_retSubset; // a temporary set for returning subset map. //moved to base class'

  //MultiMat* mm; //moved to base class'

}; //end MultiMatTypedArray


//------------ MultiMatArray Template function definitions --------------//
//
template<typename T>
inline MultiMatArray::MultiMatArray(MultiMat * m, std::string name, FieldMapping f, T * data_arr):
  mm(m), m_arrayName(name), m_fieldMapping(f), m_mapPtr(nullptr), m_subsetmap(nullptr)
{
  MultiMat::SetType* s = m->get_mapped_set(f);
  auto mapPtr = new MapType<T>(s);
  m_mapPtr = mapPtr;

  //copy data
  for (int i = 0; i < mapPtr->size(); i++) {
    (*mapPtr)[i] = data_arr[i];
  }
}

template<typename T> 
T MultiMatArray::getValue(int c)
{
  assert(m_fieldMapping == PER_CELL || m_fieldMapping == PER_MAT);
  MapType<T>* map = dynamic_cast<MapType<T>*>(m_mapPtr);
  return (*map)[c];
}

template<typename T> 
T MultiMatArray::getValue(int c, int m)
{
  assert(m_fieldMapping == PER_CELL_MAT);
  MapType<T>* map = dynamic_cast<MapType<T>*>(m_mapPtr);
  auto set_ptr = dynamic_cast<const MultiMat::MappedRelationSetType*>(map->set());
  auto from_set_size = set_ptr->fromSetSize();
  auto to_set_size = set_ptr->toSetSize();
  assert(c >= 0 && c < from_set_size);
  assert(m >= 0 && m < to_set_size);

  auto i = set_ptr->elementIndex(c, m);
  if (i == -1) //this cell does not contain this material
    return 0;
  return (*map)[i];
}

template<typename T>
const MultiMatArray::MapType<T>& MultiMatArray::getSubsetMap(int c)
{
  //A better way to do this is to use a MapBuilder to build a subset-map of the original map
  // like how SetBuilder can build a subset of the map.
  //Since that's unavailable, the placeholder is having a copy of the subset map inside the array class.
  assert(m_fieldMapping == PER_CELL_MAT);

  m_retSubset = mm->getMatInCell(c); // m_cell2matRel[c];
  auto csetsize = m_retSubset.size();
  auto offset = m_retSubset.offset();

  auto mapPtr = dynamic_cast<MapType<T>*>(m_mapPtr);

  //delete(dynamic_cast<MapType<T>*>(m_subsetmap));
  //m_subsetmap = nullptr;
  
  MapType<T>& subsetmap = * static_cast<MapType<T>*>(m_subsetmapPtr);
  subsetmap = MapType<T>(&m_retSubset);
  
  //copy the map
  int j = 0;
  for (int i = offset; i < offset + csetsize; i++)
  {
    subsetmap[j++] = (*mapPtr)[i];
  }

  assert(j == csetsize);

  return subsetmap;
}


//------------ MultiMatTypedArray Template function definitions --------------//

//template<class T>
//inline MultiMatTypedArray<T>::MultiMatTypedArray(MultiMat* m, std::string name, SetType* s, FieldMapping f)
//  :MultiMatArray(name, f), m_map(s), mm(m)
//{
//}

template<class T>
inline MultiMatTypedArray<T>::MultiMatTypedArray(MultiMat* m, std::string name, SetType* s, FieldMapping f, T * data_arr)
  : MultiMatArray(m, name, f), m_map(s)//, mm(m)
{
  //base class variables
  m_mapPtr = &m_map;
  m_subsetmapPtr = &m_subsetmap;

  //copy the array
  for (int i = 0; i < m_map.size(); i++) {
    m_map[i] = data_arr[i];
  }
}

template<class T>
inline T MultiMatTypedArray<T>::getValue(int c)
{
  assert(m_fieldMapping == PER_CELL || m_fieldMapping == PER_MAT);
  return m_map[c];
}

template<class T>
inline T MultiMatTypedArray<T>::getValue(int c, int m)
{
  assert(m_fieldMapping == PER_CELL_MAT);
  auto set_ptr = dynamic_cast<const MultiMat::MappedRelationSetType*>(m_map.set());
  auto from_set_size = set_ptr->fromSetSize();
  auto to_set_size = set_ptr->toSetSize();
  assert(c >= 0 && c < from_set_size);
  assert(m >= 0 && m < to_set_size);
  
  auto i = set_ptr->elementIndex(c, m);
  if (i == -1) //this cell does not contain this material
    return 0;
  return m_map[i];
}

template<class T>
inline void MultiMatTypedArray<T>::setValue(int c, T val)
{
  assert(m_fieldMapping == PER_CELL || m_fieldMapping == PER_MAT);
  assert(c >= 0 && c < m_map.size()); 
  m_map[c] = val;
}

template<class T>
inline void MultiMatTypedArray<T>::setValue(int c, int m, T val)
{
  assert(m_fieldMapping == PER_CELL_MAT);
  auto set_ptr = dynamic_cast<MultiMat::MappedRelationSetType*>(m_map.set());
  auto i = set_ptr->elementIndex(c, m);
  m_map[i] = val;
}

template<typename T>
const typename MultiMatTypedArray<T>::MapType& MultiMatTypedArray<T>::getSubsetMap(int c)
{
  //A better way to do this is to use a MapBuilder to build a subset-map of the original map
  // like how SetBuilder can build a subset of the map.
  //Since that's unavailable, the placeholder is having a copy of the subset map inside the array class.
  assert(m_fieldMapping == PER_CELL_MAT);

  m_retSubset = mm->getMatInCell(c); // m_cell2matRel[c];
  auto csetsize = m_retSubset.size();
  auto offset = m_retSubset.offset();

  m_subsetmap = MapType(&m_retSubset);
  
  //copy the map
  int j = 0;
  for (int i = offset; i < offset+csetsize; i++)
    m_subsetmap[j++] = m_map[i];

  assert(j == csetsize);
  return m_subsetmap;
}

template<typename T>
typename MultiMatTypedArray<T>::RangeSetType MultiMatTypedArray<T>::getIndexingSetOfCell(int c) {
  assert(0 <= c && c < mm->m_ncells);
  int start_idx = mm->m_cell2matRel_beginsVec[c];
  int end_idx = mm->m_cell2matRel_beginsVec[c + 1];
  return RangeSetType::SetBuilder().range(start_idx, end_idx);
}


//--------------- MultiMat template function definitions -----------------//






//template<class T>
//MultiMatArray<T>* axom::multimat::MultiMat::newFieldArray(std::string arr_name, FieldMapping arr_mapping)
//{
//  SetType* map_set = get_mapped_set(arr_mapping);
//  auto new_arr = new MultiMatArray<T>(this, arr_name, map_set, arr_mapping);
//  new_arr->setName(arr_name);
//  m_fieldArrayVec.push_back(new_arr);
//
//  return new_arr;
//}

template<class T>
MultiMatTypedArray<T>* MultiMat::newFieldArray(std::string arr_name, FieldMapping arr_mapping, T* data_arr)
{
  //SetType* map_set = get_mapped_set(arr_mapping);
  auto new_arr_ptr = new MultiMatTypedArray<T>(this, arr_name, get_mapped_set(arr_mapping), arr_mapping, data_arr);
  //MultiMatArray* new_arr = new MultiMatArray(this, arr_name, arr_mapping, data_arr);
  m_fieldArrayVec.push_back(new_arr_ptr);
  
  return new_arr_ptr;
}


template<typename T>
MultiMatTypedArray<T>* MultiMat::getFieldArray(std::string arr_name)
{
  MultiMatArray* arr = getFieldArray(arr_name);
  MultiMatTypedArray<T>* typed_arr = static_cast<MultiMatTypedArray<T>*>(arr);
  return typed_arr;
}

template<class T>
T MultiMat::getFieldValue(std::string field_name, int cell_i, int mat_i)
{
  auto arr_ptr = MM_CAST_TO(T, getFieldArray(field_name));
  assert(arr_ptr != nullptr);
  return arr_ptr->getValue(cell_i, mat_i);
}


template<class T>
T MultiMat::getFieldValue(std::string field_name, int i)
{
  auto arr_ptr = MM_CAST_TO(T, getFieldArray(field_name));
  assert(arr_ptr != nullptr);
  return arr_ptr->getValue(i);
}


template<typename T>
const MultiMat::MapType<T>& MultiMat::getMap(std::string field_name)
{
  MultiMatArray* fieldArr = getFieldArray(field_name);
  const MultiMatArray::MapBaseType* mapBasePtr = fieldArr->getMap();
  const MultiMatArray::MapType<T>* mapPtr = dynamic_cast<const MultiMatArray::MapType<T>*>(mapBasePtr);
  return *mapPtr;
}

} //end namespace multimat
} //end namespace axom


#endif