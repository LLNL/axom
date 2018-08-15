/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/**
 * \file multimat.hpp
 *
 * \brief Contains the MultiMat library header and its template implementation
 *
 */
#ifndef MULTIMAT_H_
#define MULTIMAT_H_


#include "slam/RangeSet.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/StaticRelation.hpp"
#include "slam/DynamicVariableRelation.hpp"
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
enum class DataTypeSupported
{
  TypeUnknown, TypeInt, TypeDouble, TypeFloat, TypeUnsignChar
};

/**
 * \class MultiMat
 *
 * \brief A multimaterial data management class that provides storage in various
 * layouts (dense/sparse, and material-dominant/cell-dominant).
 *
 */
class MultiMat
{
private:
  // SLAM Set type definitions
  using SetType = slam::Set;
  using RangeSetType   = slam::RangeSet;
  using BivariateSetType = slam::BivariateSet;
  using SetPosType  = SetType::PositionType;
  using SetElemType = SetType::ElementType;
  // SLAM Relation typedef
  using STLIndirection = policies::STLVectorIndirection<SetPosType, SetPosType>;
  using VariableCardinality =
          policies::VariableCardinality<SetPosType, STLIndirection>;
  using StaticVariableRelationType = slam::StaticRelation<
          VariableCardinality, STLIndirection, RangeSetType, RangeSetType>;
  using DynamicVariableRelationType = slam::DynamicVariableRelation;
  using OrderedSetType = slam::OrderedSet<
          policies::RuntimeSize<SetType::PositionType>,
          policies::RuntimeOffset<SetType::PositionType>,
          policies::StrideOne<SetType::PositionType>,
          STLIndirection,
          policies::NoSubset>;
  using ProductSetType = slam::ProductSet;
  // SLAM RelationSet for the set of non-zero cell to mat variables
  using RelationSetType = slam::RelationSet<StaticVariableRelationType>;
  using RelationSetDynType = slam::RelationSet<DynamicVariableRelationType>;

  // SLAM Map type
  using MapStrideType = slam::policies::RuntimeStride<SetPosType>;
  using MapBaseType = slam::MapBase;
  template <typename T>
  using MapType = slam::Map<T, MapStrideType>;
  template <typename T>
  using BivariateMapType = slam::BivariateMap<T, MapStrideType>;
  template<typename T, typename M>
  using SubMap = slam::SubMap<T, M, MapStrideType>;

public:
  //public typedef
  template <typename T>
  using Field1D = MapType<T>;
  template <typename T>
  using Field2D = BivariateMapType<T>;
  template <typename T>
  using SubField = SubMap<T, Field2D<T> >;
  using IndexSet = RangeSetType; //For returning set of SparseIndex
  using IdSet = OrderedSetType;  //For returning set of DenseIndex

  //Constructors

  /**
   * \brief Constructor to create a new MultiMat object
   *
   * \param data_layout  Select material-dominant or cell-dominant layout to
   *                     store the data in. Default is cell-dominant.
   * \param sparcity_layout Select dense or sparse layout to store the data in.
   *                        Default is sparse.
   */
  MultiMat(DataLayout data_layout = DataLayout::CELL_CENTRIC,
           SparcityLayout sparcity_layout = SparcityLayout::SPARSE);
  /** Destructor **/
  ~MultiMat();
  /** Copy constructor (Deep copy). **/
  MultiMat(const MultiMat&);
  /** Assignment operator **/
  MultiMat& operator=(const MultiMat&);

  //Set-up functions
  /**
   * \brief Set the number of materials
   * \pre num_cells > 0
   */
  void setNumberOfMat(int num_cells);
  /**
   * \brief Set the number of cells
   * \pre num_mats > 0
   */
  void setNumberOfCell(int num_mats);

  /**
   * \brief Set the cell-material relation.
   *
   * \detail This function takes the input of boolean vectors and set-up
   * internal data structure. The relation_info vector is indexed by the cell
   * and material index: \n
   *   `idx = mat * num_cell + cell`  \n
   * or `idx = cell * num_mats + mat` \n
   * for material-centric and cell-centric layout, respectively.\n
   * Entry at `idx` should containing 'true' where materials are presented in
   * the cell.\n
   * The number of materials and cell must be set prior to calling this function
   * with setNumberOfMat(int) and setNumberOfCell(int)
   *
   * \param A boolean vector of size num_mats * num_cells containing information
   * on if a materials is present in a cell.
   *
   */
  void setCellMatRel(std::vector<bool>& relation_info);



  //functions related to field arrays

  /**
   * \brief Add a field to the MultiMat object
   *
   * \tparam T The data type (double, float...) of the field
   * \param field_name The name of the field, used to retrieve the field later
   * \param field_mapping
   * \param data_array The array containing data to the field. The length of the
   *            array should be `num_mats * num_cells * ncomp` if the current
   *            format is Dense, or `num_nonzero * ncomp` if the current format
   *            is Sprase.
   * \param (optional) ncomp The number of component of the field. Default is 1
   * \return int the index of the field, can be used to retrieve the field later
   */
  template<class T>
  int addField(const std::string& field_name, FieldMapping field_mapping,
               T* data_array, int ncomp = 1);

  /**
   * \brief Set the volume fraction field
   * \detail volume fraction field is assumed to be a double. Its field index is
   * always 0, and the name of the field is "Volfrac"
   *
   * \param data_array the array containing the volumn fraction information
   * \return int the volume fraction field index, which is always zero.
   */
  int setVolfracField(double* data_array);

  /**
   * \brief Search for and return the field index of the field given its name.
   *
   * \param field_name the name of the field
   * \return int the index of the field
   */
  int getFieldIdx(const std::string& field_name) const;

  /**
   * \brief Search for and return the field given the field name.
   * \detail the field is of type Field1D, containing an entry for each cell
   * or material. To retrieve a field of type Field2D, use get2dField()
   *
   * \tparam T The data type of the field
   * \param field_name the name of the field
   * \return Field1D<T>& the field reference
   */

  template<typename T>
  Field1D<T>& get1dField(const std::string& field_name);
  /**
   * \brief Search for and return the field given the field name.
   * \detail the field is of type Field2D, containing an entry for each cell and
   * each material. To retrieve a field of type Field1D, use get1dField()
   *
   * \tparam T The data type of the field
   * \param field_name the name of the field
   * \return Field1D<T>& the field reference
   */
  template<typename T>
  Field2D<T>& get2dField(const std::string& field_name);

  /**
   * \brief Get the volume fraction field
   */
  Field2D<double>& getVolfracField();

  /**
   * \brief Get a set of index for a Subfield.
   * \detail for a cell-dominant layout, this is equivalent to getting the set
   * of material presented in a cell. Vise versa, for a material-dominant
   * layout, this returns a set of cell that contains a material.\n
   * getMatInCell() and getCellContainingMat() are layout specific
   * calls for this function.
   *
   * \param idx the index of the subfield
   * \return IdSet A set of index in the subfield
   */
  IdSet getSubfieldIndex(int idx);

  /** Cell-dominant version of getSubfieldIndex() **/
  IdSet getMatInCell(int cell_id);
  /** Material-dominant version of getSubfieldIndex() **/
  IdSet getCellContainingMat(int mat_id);

  /**
   * \brief Get a set to index a Subfield in a Field2D.
   * \detail This returns a set of indices that can be used to directly access
   * a subfield in a Field2D without having to retrieve a Subfield first.\n
   * Accessing the Field2D using this indexing set should be done with the
   * bracket operator.
   *
   * e.g.
   * \code
   *   IndexSet indexing_set = getSubfieldIndexingSet( idx );
   *   for( int i : indexing_set );
   *       value = field2d[ i ];
   * \endcode
   *
   * Or to account for multiple components:
   * \code
   *   IndexSet indexing_set = getSubfieldIndexingSet( idx );
   *   for( int i : indexing_set );
   *       for( int c = 0; c < numComp; c++)
   *           value = field2d[ i * numComp + c ];
   * \endcode
   *
   * getIndexingSetOfCell() and getIndexingSetOfMat() are layout specific
   * calls for this function.
   *
   * \param idx The index of the subfield
   * \return IndexSet A Set of index for the subfield
   */
  IndexSet getSubfieldIndexingSet(int idx);

  /** Cell-dominant version of getSubfieldIndexingSet() **/
  IndexSet getIndexingSetOfCell(int cell_id);
  /** Material-dominant version of getSubfieldIndexingSet() **/
  IndexSet getIndexingSetOfMat(int mat_id);

  /** Return the number of material this object holds **/
  int getNumberOfMaterials() const { return m_nmats; };
  /** Return the number of cells this object holds **/
  int getNumberOfCells() const { return m_ncells; }




  //Layout modification functions

  /** Convert the data to be stored in cell-dominant layout. **/
  void convertLayoutToCellDominant();
  /** Convert the data to be stored in material-dominant layout. **/
  void convertLayoutToMaterialDominant();
  /** Convert the data stored from a material-dominant layout to cell-dominant
   *  layout, or vise versa. **/
  void transposeData();

  /** Convert the data to be stored in sparse/compact layout **/
  void convertLayoutToSparse();
  /** Convert the data to be stored in dense layout **/
  void convertLayoutToDense();

  /** Convert the data to be stored in the specified layout. **/
  void convertLayout(DataLayout, SparcityLayout);
  /** Return the data layout used currently **/
  DataLayout getDataLayout() const { return m_dataLayout; }
  /** Return the sparcity layout used currently **/
  SparcityLayout getSparcityLayout() const { return m_sparcityLayout; }
  /** Return the data layout used currently as string **/
  std::string getDataLayoutAsString() const;
  /** Return the sparcity layout used currently as string **/
  std::string getSparcityLayoutAsString() const;
  /** Return true if the current layout is sparse/compact, false otherwise. **/
  bool isSparse() const { return m_sparcityLayout == SparcityLayout::SPARSE; }
  /** Return true if the current layout is dense, false otherwise. **/
  bool isDense() const { return m_sparcityLayout == SparcityLayout::DENSE; }
  /** Return true if the current layout is cell-dominant, false otherwise. **/
  bool isCellDom() const { return m_dataLayout == DataLayout::CELL_CENTRIC; }
  /** Return true if the current layout is material-dominant, false otherwise.*/
  bool isMatDom() const { return m_dataLayout == DataLayout::MAT_CENTRIC; }
  /**
   * \brief Get the FieldMapping for a field.
   *
   * A FieldMapping of a field describes if there is an entry for each cell, for
   * each material, or for each cell x material.
   *
   * \param field_idx the index of the field
   */
  FieldMapping getFieldMapping(int field_idx) const {
    return m_fieldMappingVec[field_idx];
  }




  //Dynamic mode functions

  /** Convert the layout to use dynamic layout. In dynamic layout, users can
   *  call addEntry() or removeEntry() to change the materials in a cell.
   */
  void convertToDynamic();
  /** Convert the layout to use static layout. **/
  void convertToStatic();

  /**
   * \brief Add a material in a cell.
   *
   * In a cell-dominant layout, firstIdx is the cell index, secondIdx is the
   * material index. And it's reverse for material-dominant layout.
   *
   * \return true if the material is added in the. False if it already exists.
   */
  bool addEntry(int firstIdx, int secondIdx);
  /**
   * \brief Remove a material from a cell
   *
   * In a cell-dominant layout, firstIdx is the cell index, secondIdx is the
   * material index. And it's reverse for material-dominant layout.
   *
   * \return true if a material is removed. False if it was not in the cell.
   */
  bool removeEntry(int firstIdx, int secondIdx);


  /** Print the detail of this object */
  void print() const;
  /**
   * \brief  Return true if the object is valid, false otherwise.
   *
   * \param verboseOutput (Optional) if true, this function will print out its
   * valid/invalid detail.
   */
  bool isValid(bool verboseOutput = false) const;


private: //private functions
  //Return the Set pointer associalted with the given FieldMapping
  SetType* get_mapped_set(FieldMapping);
  //Return the BivariateSet used for the current layout.
  BivariateSetType* get_mapped_biSet();

  //helper functions
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
  MapBaseType* helper_copyField(const MultiMat&, int map_i);

private:
  unsigned int m_ncells, m_nmats;
  DataLayout m_dataLayout;
  SparcityLayout m_sparcityLayout;

  //slam set variables
  RangeSetType m_cellSet;
  RangeSetType m_matSet;
  //slam relation variables
  std::vector<SetPosType> m_cellMatRel_beginsVec;
  std::vector<SetPosType> m_cellMatRel_indicesVec;
  StaticVariableRelationType* m_cellMatRel;
  DynamicVariableRelationType* m_cellMatRelDyn;
  //slam bivariateSet variables
  RelationSetType* m_cellMatNZSet;
  ProductSetType* m_cellMatProdSet;

  //vector of information for each fields
  std::vector<std::string> m_arrNameVec;
  std::vector<FieldMapping> m_fieldMappingVec;
  std::vector<MapBaseType*> m_mapVec;
  std::vector<DataTypeSupported> m_dataTypeVec;

  //To store the static layout information when converting to dynamic.
  struct Layout
  {
    DataLayout data_layout;
    SparcityLayout sparcity_layout;
  };

  Layout m_static_layout; //Layout used during the static mode
  bool m_dynamic_mode; //true if in dynamic mode

}; //end MultiMat class



//--------------- MultiMat template function definitions -----------------//


/* Helper function to create a new slam::Map or slam::BivariateMap and add
 * it into the list of maps */
template<class T>
int MultiMat::addFieldArray_impl(const std::string& arr_name,
                                 FieldMapping arr_mapping, T* data_arr,
                                 int stride)
{
  int new_arr_idx = m_mapVec.size();

  if (arr_mapping == FieldMapping::PER_CELL_MAT)
  {
    BivariateSetType* s = get_mapped_biSet();
    SLIC_ASSERT(s != nullptr);

    Field2D<T>* new_map_ptr = new Field2D<T>(s, T(), stride);
    new_map_ptr->copy(data_arr);
    m_mapVec.push_back(new_map_ptr);
  }
  else
  {
    SLIC_ASSERT(arr_mapping == FieldMapping::PER_CELL ||
                arr_mapping == FieldMapping::PER_MAT);
    SetType* s = get_mapped_set(arr_mapping);
    Field1D<T>* new_map_ptr = new Field1D<T>(s, T(), stride);

    //copy data
    int i = 0;
    for (auto iter = new_map_ptr->begin() ; iter != new_map_ptr->end() ; iter++)
      for (auto s = 0 ; s < stride ; ++s)
        iter(s) = data_arr[i++];

    m_mapVec.push_back(new_map_ptr);
  }

  m_arrNameVec.push_back(arr_name);
  m_fieldMappingVec.push_back(arr_mapping);

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

  SLIC_ASSERT(m_mapVec.size() == m_arrNameVec.size());
  SLIC_ASSERT(m_mapVec.size() == m_fieldMappingVec.size());
  SLIC_ASSERT(m_mapVec.size() == m_dataTypeVec.size());
  return new_arr_idx;
}



template<class T>
int MultiMat::addField(const std::string& arr_name, FieldMapping arr_mapping,
                       T* data_arr, int stride)
{
  SLIC_ASSERT(stride > 0);

  //make sure the name does not conflict
  int fieldIdx = getFieldIdx(arr_name);
  if (fieldIdx == 0 && m_mapVec[0] == nullptr)
  { //this is the vol frac array. call setVolfrac instead
    SLIC_ASSERT(arr_mapping == FieldMapping::PER_CELL_MAT);
    SLIC_ASSERT(stride == 1);
    SLIC_ASSERT(data_arr != nullptr);
    setVolfracField(data_arr);
    return 0;
  }
  else if(fieldIdx > 0)
  {
    // There is already an array with the current name. And it's not Volfrac.
    // Don't add the new field.
    SLIC_ASSERT(false);

    return fieldIdx;
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

  if(fieldIdx < 0)
    throw std::invalid_argument("No field with this name is found");

  if (m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL ||
      m_fieldMappingVec[fieldIdx] == FieldMapping::PER_MAT)
  {
    return *dynamic_cast<Field1D<T>*>(m_mapVec[fieldIdx]);
  }
  else
  {
    SLIC_ASSERT(m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL_MAT);

    //Right now we're allowing Field2D (BivariateMap) to be returned as
    // a Field1D (Map) so it can be accessed like a 1d array, but the
    // indexing information would be lost.
    Field2D<T>* map_2d = dynamic_cast<Field2D<T>*>(m_mapVec[fieldIdx]);
    return *(map_2d->getMap());
  }
}


template<typename T>
MultiMat::Field2D<T>& MultiMat::get2dField(const std::string& field_name)
{
  int fieldIdx = getFieldIdx(field_name);

  if (fieldIdx < 0)
    throw std::invalid_argument("No field with this name is found");

  SLIC_ASSERT(m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL_MAT);
  return *dynamic_cast<Field2D<T>*>(m_mapVec[fieldIdx]);
}


template<typename DataType>
void MultiMat::convertToSparse_helper(int map_i)
{
  SLIC_ASSERT(m_sparcityLayout != SparcityLayout::SPARSE);

  MapBaseType* mapPtr = m_mapVec[map_i];

  //Skip if no volume fraction array is set-up
  if (map_i == 0 && mapPtr == nullptr)
    return;

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(m_cellMatRel->totalSize()*stride);
  int idx = 0;
  for (int i = 0 ; i < m_cellMatRel->fromSetSize() ; ++i)
  {
    auto relset = (*m_cellMatRel)[i];
    auto submap = old_map(i);
    for (int j = 0 ; j < relset.size() ; ++j)
    {
      for (int s = 0 ; s < stride ; ++s)
      {
        arr_data[idx++] = submap[relset[j]*stride + s];
      }
    }
  }
  SLIC_ASSERT(idx == m_cellMatRel->totalSize()*stride);
  Field2D<DataType>* new_field =
    new Field2D<DataType>(m_cellMatNZSet, DataType(), stride);
  new_field->copy(arr_data.data());

  delete m_mapVec[map_i];
  m_mapVec[map_i] = new_field;
}


template<typename DataType>
void MultiMat::convertToDense_helper(int map_i)
{
  SLIC_ASSERT(m_sparcityLayout != SparcityLayout::DENSE);

  MapBaseType* mapPtr = m_mapVec[map_i];

  //Skip if no volume fraction array is set-up
  if (map_i == 0 && mapPtr == nullptr)
  {
    return;
  }

  Field2D<DataType>& old_map = *dynamic_cast<Field2D<DataType>*>(mapPtr);
  int stride = old_map.stride();
  std::vector<DataType> arr_data(m_cellMatProdSet->size()*stride);
  for (int i = 0 ; i < old_map.firstSetSize() ; ++i)
  {
    for (auto iter = old_map.begin(i) ; iter != old_map.end(i) ; ++iter)
    {
      int elem_idx = i * old_map.secondSetSize() + iter.index();
      for (int c = 0 ; c < stride ; ++c)
      {
        arr_data[elem_idx*stride + c] = iter.value(c);
      }
    }
  }

  Field2D<DataType>* new_field =
    new Field2D<DataType>(m_cellMatProdSet, DataType(), stride);
  new_field->copy(&arr_data[0]);

  delete m_mapVec[map_i];
  m_mapVec[map_i] = new_field;
}



} //end namespace multimat
} //end namespace axom


#endif
