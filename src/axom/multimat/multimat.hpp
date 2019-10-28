// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file multimat.hpp
 *
 * \brief Contains the MultiMat library header and its template implementation
 *
 */
#ifndef MULTIMAT_H_
#define MULTIMAT_H_


#include "axom/slam.hpp"

#include <vector>
#include <cassert>
#include <stdexcept>

namespace axom
{
namespace multimat
{

enum class FieldMapping { PER_CELL, PER_MAT, PER_CELL_MAT };
enum class DataLayout { CELL_DOM, MAT_DOM };
enum class SparsityLayout { SPARSE, DENSE };
enum class DataTypeSupported
{
  TypeUnknown, TypeInt, TypeDouble, TypeFloat, TypeUnsignChar
};


//forward class declarations
template<typename T, typename BiSet>
class MMField2D;
template<typename Field2DType>
class MMSubField2D;
template <typename T, DataLayout D, typename B>
class MMField2DTemplated;


/**
 * \class MultiMat
 *
 * \brief A multimaterial data management class that provides storage in various
 * layouts (dense/sparse, and material-dominant/cell-dominant).
 *
 */
class MultiMat
{
public:
  using BivariateSetType = slam::BivariateSet<>;//RangeSetType, RangeSetType>;
  using ProductSetType = slam::ProductSet<>;//RangeSetType, RangeSetType>;

protected:
  // SLAM Set type definitions
  using SetPosType  = slam::DefaultPositionType;
  using SetElemType = slam::DefaultPositionType;
  using SetType = slam::Set<SetPosType,SetElemType>;
  using RangeSetType = slam::RangeSet<SetPosType,SetElemType>;
  // SLAM Relation typedef
  template<typename T>
  using IndPolicy = slam::policies::STLVectorIndirection<SetPosType, T>;
  using VariableCardinality =
        slam::policies::VariableCardinality<SetPosType, IndPolicy<SetElemType> >;
  using StaticVariableRelationType = slam::StaticRelation<
          SetPosType, SetElemType, VariableCardinality, 
          IndPolicy<SetElemType>, RangeSetType, RangeSetType>;

  using DynamicVariableRelationType =
          slam::DynamicVariableRelation<SetPosType,SetElemType>;
  using OrderedSetType = slam::OrderedSet<
          SetPosType,SetElemType,
          slam::policies::RuntimeSize<SetPosType>,
          slam::policies::RuntimeOffset<SetPosType>,
          slam::policies::StrideOne<SetPosType>,
          IndPolicy<SetElemType>,
          slam::policies::NoSubset>;

  // SLAM Map type
  using MapStrideType = slam::policies::RuntimeStride<SetPosType>;
  using MapBaseType = slam::MapBase;

  template <typename T>
  using MapType = slam::Map<T, SetType, IndPolicy<T>, MapStrideType>;

  template <typename T, typename BSet = BivariateSetType>
  using BivariateMapType =
          slam::BivariateMap<T,BSet,IndPolicy<T>,MapStrideType>;

public:
  using SparseRelationType = StaticVariableRelationType;

  // SLAM RelationSet for the set of non-zero cell to mat variables
  using RelationSetType = 
        slam::RelationSet<StaticVariableRelationType>;//, RangeSetType
  using RelationSetDynType = slam::RelationSet<DynamicVariableRelationType>;

public:
  //1D Field
  template <typename T>
  using Field1D = MapType<T>;

  //2D Field
  //old
  //template <typename T, typename BSet = BivariateSetType>
  //using Field2D = BivariateMapType<T,BSet>;
  template <typename T, typename BSet = BivariateSetType>
  using Field2D = MMField2D<T, BSet>;
  //special
  template <typename T>
  using SparseField2D = MMField2D<T, RelationSetType>;
  template <typename T>
  using DenseField2D = MMField2D<T, ProductSetType>;

  template <typename T, DataLayout D, typename B>
  using Field2DTemplated = MMField2DTemplated<T, D, B>;


  template<typename Field2DType>
  using SubField = MMSubField2D<Field2DType>;
  


  using IndexSet = RangeSetType; //For returning set of SparseIndex
  using IdSet = OrderedSetType;  //For returning set of DenseIndex

  //Constructors

  /**
   * \brief Constructor to create a new MultiMat object
   *
   * \param data_layout  Select material-dominant or cell-dominant layout to
   *                     store the data in. Default is cell-dominant.
   * \param sparsity_layout Select dense or sparse layout to store the data in.
   *                        Default is sparse.
   */
  MultiMat(DataLayout data_layout = DataLayout::CELL_DOM,
           SparsityLayout sparsity_layout = SparsityLayout::SPARSE);
  /** Destructor **/
  ~MultiMat();
  /** Copy constructor (Deep copy). **/
  MultiMat(const MultiMat&);
  /** Assignment operator **/
  MultiMat& operator=(const MultiMat&);

  //Set-up functions
  /**
   * \brief Set the number of materials
   * \pre num_mats > 0
   */
  void setNumberOfMaterials(int num_mats);
  /**
   * \brief Set the number of cells
   * \pre num_cells > 0
   */
  void setNumberOfCells(int num_cells);

  /// \brief Returns a pointer to the dense 2d field set
  const ProductSetType* getDense2dFieldSet(DataLayout layout) const
  {
    if (layout == DataLayout::CELL_DOM)
      return m_cellMatProdSet;
    else
      return m_matCellProdSet;
  }

  /// \brief Returns a pointer to the sparse 2d field set
  const RelationSetType* getSparse2dFieldSet(DataLayout layout) const
  {
    if (layout == DataLayout::CELL_DOM)
      return m_cellMatNZSet;
    else
      return m_matCellNZSet;
  }

  
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
   * with setNumberOfMaterials(int) and setNumberOfCells(int)
   *
   * \param A boolean vector of size num_mats * num_cells containing information
   * on if a materials is present in a cell.
   *
   */
  void setCellMatRel(std::vector<bool>& relation_info, DataLayout layout);



  //functions related to fields

  int getNumberOfFields() const { return m_mapVec.size(); }

  /**
   * \brief Add a field to the MultiMat object
   *
   * \tparam T The data type (double, float...) of the field
   * \param field_name The name of the field, used to retrieve the field later
   * \param field_mapping
   * \param data_layout
   * \param sparsity_layout
   * \param data_array The array containing data to the field. The length of the
   *            array should be `num_mats * num_cells * ncomp` if the current
   *            format is Dense, or `num_nonzero * ncomp` if the current format
   *            is Sprase.
   * \param (optional) ncomp The number of component of the field. Default is 1
   * \return int the index of the field, can be used to retrieve the field later
   */
  template<class T>
  int addField(const std::string& field_name, FieldMapping field_mapping,
    DataLayout data_layout, SparsityLayout sparsity_layout, 
    T* data_array, int ncomp = 1);

private:
  template<typename T>
  int addFieldArray_impl(const std::string&, FieldMapping, DataLayout, SparsityLayout, T*, int);

public:

  /**
   * \brief Set the volume fraction field
   * \detail volume fraction field is assumed to be a double. Its field index is
   * always 0, and the name of the field is "Volfrac"
   *
   * \param data_array the array containing the volumn fraction information
   * \return int the volume fraction field index, which is always zero.
   */
  int setVolfracField(double* data_array, DataLayout layout, SparsityLayout sparsity);

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

  template<typename T>
  DenseField2D<T> getDense2dField(const std::string& field_name);

  template<typename T>
  DenseField2D<T> getSparse2dField(const std::string& field_name);

  template<typename T, DataLayout D, typename B>
  Field2DTemplated<T, D, B> getTemplated2DField(const std::string& field_name);

  template<typename T, typename BSetType>
  slam::BivariateMap<T,BSetType> get2dFieldAsSlamBivarMap(const std::string& field_name);

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
  IndexSet getSubfieldIndexingSet(int idx, DataLayout layout, SparsityLayout sparsity);
  

  /** Cell-dominant version of getSubfieldIndexingSet() **/
  IndexSet getIndexingSetOfCell(int cell_id, SparsityLayout sparsity);
  /** Material-dominant version of getSubfieldIndexingSet() **/
  IndexSet getIndexingSetOfMat(int mat_id, SparsityLayout sparsity);

  /** Return the number of material this object holds **/
  int getNumberOfMaterials() const { return m_nmats; };
  /** Return the number of cells this object holds **/
  int getNumberOfCells() const { return m_ncells; }




  //Layout modification functions
  
  void convertFieldLayout(int field_idx, SparsityLayout, DataLayout);

  void convertFieldToSparse(int field_idx);
  void convertFieldToDense(int field_idx);
  SparsityLayout getFieldSparsityLayout(int field_idx);

  void transposeField(int field_idx);
  void convertFieldToMatDom(int field_idx);
  void convertFieldToCellDom(int field_idx);
  DataLayout getFieldDataLayout(int field_idx);

  std::string getFieldDataLayoutAsString(int field_idx) const;
  std::string getFieldSparsityLayoutAsString(int field_idx) const;

  /** Convert the data to be stored in the specified layout. **/
  void convertLayout(DataLayout, SparsityLayout);
  /** Convert the data to be stored in sparse/compact layout **/
  void convertLayoutToSparse();
  /** Convert the data to be stored in dense layout **/
  void convertLayoutToDense();
  /** Convert the data to be stored in cell-dominant layout. **/
  void convertLayoutToCellDominant();
  /** Convert the data to be stored in material-dominant layout. **/
  void convertLayoutToMaterialDominant();
  

  
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

protected:
  //Return the Set pointer associalted with the given FieldMapping or field idx
  SetType* get_mapped_set(FieldMapping);
  SetType* get_mapped_set(int field_idx);
  //Return the BivariateSet used for the specified layouts or the field
  BivariateSetType* get_mapped_biSet(DataLayout, SparsityLayout);
  BivariateSetType* get_mapped_biSet(int field_idx);
  //Return the relation for the specified field idx or the mapping
  StaticVariableRelationType* getRel(int field_idx);

private: //private functions
  //Given a relation (cell->mat or mat->cell), create the other relation
  void makeOtherRelation();

  
  //helper functions
  template<typename DataType>
  void convertToSparse_helper(int map_i);
  template<typename DataType>
  void convertToDense_helper(int map_i);

  template<typename DataType>
  void transposeField_helper(int field_idx);
  
  template<typename T>
  MapBaseType* helper_copyField(const MultiMat&, int map_i);

private:
  unsigned int m_ncells, m_nmats;
  
  //slam set variables
  RangeSetType m_cellSet;
  RangeSetType m_matSet;
  //slam relation variables (for sparse-layout fields)
  //Depending on the layout of the field, each field can be mapped to different Relations
  //Relation can be nullptr if no field is using said relation
  //cell to mat relations
  std::vector<SetPosType> m_cellMatRel_beginsVec;
  std::vector<SetPosType> m_cellMatRel_indicesVec;
  StaticVariableRelationType* m_cellMatRel;
  DynamicVariableRelationType* m_cellMatRelDyn;
  //mat to cell relations
  std::vector<SetPosType> m_matCellRel_beginsVec;
  std::vector<SetPosType> m_matCellRel_indicesVec;
  StaticVariableRelationType* m_matCellRel;
  DynamicVariableRelationType* m_matCellRelDyn;


  //slam bivariateSet variables
  //cell to mat versions
  RelationSetType* m_cellMatNZSet;  //for sparse layout
  ProductSetType* m_cellMatProdSet; //for dense layout
  //mat to cell vertions
  RelationSetType* m_matCellNZSet;  //for sparse layout
  ProductSetType* m_matCellProdSet; //for dense layout

  //std::vector of information for each fields
  std::vector<std::string> m_fieldNameVec;
  std::vector<FieldMapping> m_fieldMappingVec;
  std::vector<MapBaseType*> m_mapVec;
  std::vector<DataTypeSupported> m_dataTypeVec;
  std::vector<DataLayout> m_fieldDataLayoutVec;
  std::vector<SparsityLayout> m_fieldSparsityLayoutVec;

  //To store the static layout information when converting to dynamic.
  struct Layout
  {
    DataLayout data_layout;
    SparsityLayout sparsity_layout;
  };

  //Layout used during the static mode
  std::vector<Layout> m_layout_when_static;
  Layout m_static_layout; 

  bool m_dynamic_mode; //true if in dynamic mode

  template<typename T, typename S>
  friend class MMField2D; // every type of MMField2D is a friend
  template <typename Field2D>
  friend class MMSubField2D; // every type of MMSubField2D is a friend
  template <typename T, DataLayout D, typename S>
  friend class MMField2DTemplated;
}; //end MultiMat class





//--------------- MultiMat template function definitions -----------------//


template<class T>
int MultiMat::addField(const std::string& arr_name, FieldMapping arr_mapping,
  DataLayout data_layout, SparsityLayout sparsity_layout, T* data_arr, int stride)
{
  SLIC_ASSERT(stride > 0);

  //make sure the name does not conflict
  int fieldIdx = getFieldIdx(arr_name);
  if (fieldIdx == 0 && m_mapVec[0] == nullptr)
  { //this is the vol frac array. call setVolfrac instead
    SLIC_ASSERT(arr_mapping == FieldMapping::PER_CELL_MAT);
    SLIC_ASSERT(stride == 1);
    SLIC_ASSERT(data_arr != nullptr);
    setVolfracField(data_arr, data_layout, sparsity_layout);
    return 0;
  }
  else if (fieldIdx > 0)
  {
    // There is already an array with the current name. And it's not Volfrac.
    // Don't add the new field.
    SLIC_ASSERT(false);

    return fieldIdx;
  }
  else
  {
    //No field with this name. Proceed with adding the field
    return addFieldArray_impl<>(arr_name, arr_mapping, data_layout, sparsity_layout, data_arr, stride);
  }
}

template<typename T>
int MultiMat::addFieldArray_impl(const std::string& field_name, FieldMapping field_mapping,
  DataLayout data_layout, SparsityLayout sparsity_layout, T* data_arr, int stride)
{
  unsigned int new_arr_idx = m_mapVec.size();

  if (field_mapping == FieldMapping::PER_CELL_MAT)
  {
    BivariateSetType* s = get_mapped_biSet(data_layout, sparsity_layout);
    SLIC_ASSERT(s != nullptr);

    //old field2d
    //Field2D<T>* new_map_ptr = new Field2D<T>(s, T(), stride);
    //new_map_ptr->copy(data_arr);
    Field2D<T>* new_map_ptr = new Field2D<T>(*this, s,
      field_name, data_arr, stride);

    m_mapVec.push_back(new_map_ptr);

    
  }
  else
  {
    SLIC_ASSERT(field_mapping == FieldMapping::PER_CELL ||
      field_mapping == FieldMapping::PER_MAT);
    SetType* s = get_mapped_set(field_mapping);
    Field1D<T>* new_map_ptr = new Field1D<T>(s, T(), stride);

    //copy data
    int i = 0;
    for (auto iter = new_map_ptr->begin(); iter != new_map_ptr->end(); iter++)
      for (auto str = 0; str < stride; ++str)
        iter(str) = data_arr[i++];

    m_mapVec.push_back(new_map_ptr);
  }

  m_fieldNameVec.push_back(field_name);
  m_fieldMappingVec.push_back(field_mapping);
  m_fieldDataLayoutVec.push_back(data_layout);
  m_fieldSparsityLayoutVec.push_back(sparsity_layout);

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

  SLIC_ASSERT(m_mapVec.size() == m_fieldNameVec.size());
  SLIC_ASSERT(m_mapVec.size() == m_dataTypeVec.size());
  SLIC_ASSERT(m_mapVec.size() == m_fieldMappingVec.size());
  SLIC_ASSERT(m_mapVec.size() == m_fieldDataLayoutVec.size());
  SLIC_ASSERT(m_mapVec.size() == m_fieldSparsityLayoutVec.size());

  return new_arr_idx;
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
    auto* map_2d = dynamic_cast<Field2D<T>*>(m_mapVec[fieldIdx]);
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

template<typename T>
MultiMat::DenseField2D<T> MultiMat::getDense2dField(const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto& bmap = get2dField<T>(field_name);
  
  //create instance of that map
  ProductSetType* prod_set;
  if (m_fieldDataLayoutVec[getFieldIdx(field_name)] == DataLayout::CELL_DOM)
    prod_set = m_cellMatProdSet;
  else
    prod_set = m_matCellProdSet;

  DenseField2D<T> typedBMap(*this, prod_set, field_name, bmap.getMap()->data().data());

  return typedBMap;
}

template<typename T, DataLayout D, typename B>
MultiMat::Field2DTemplated<T, D, B> MultiMat::getTemplated2DField(const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto& bmap = get2dField<T>(field_name);
  Field2DTemplated<T, D, B> typedBMap(*this, field_name, bmap.getMap()->data().data());

  return typedBMap;
}


// Warning: The return type uses a compile time stride of one!
template<typename T, typename BSetType>
slam::BivariateMap<T,BSetType> MultiMat::get2dFieldAsSlamBivarMap(const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto& bmap = get2dField<T>(field_name);

  auto bset = bmap.template getBivariateSet<BSetType, RelationSetType>();

  // Create instance of templated BivariateMap
  slam::BivariateMap<T,BSetType> typedBMap(new BSetType(bset));
  typedBMap.setManagesBSetPtr(true);
  
  // Copy data from original map to templated map
  // WARNING: Map and BMap should take a pointer to the data
  // instead of copying it!
  typedBMap.copy(bmap.getMap()->data().data());

  return typedBMap;
}



} //end namespace multimat
} //end namespace axom


#include "axom/multimat/mmfield.hpp"

#endif
