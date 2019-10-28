
#ifndef MMFIELD_H_
#define MMFIELD_H_

#include "axom/multimat/multimat.hpp"
#include "axom/multimat/mmsubfield.hpp"


namespace axom {
namespace multimat {

/**
 * Class for MultiMat 2D Fields
 * Ideally this class would have template specialization for each type of 
 * layout (dense/sparse, mat/cell dom, and maybe more). 
 */
template<typename DataType, typename BiSet = MultiMat::BivariateSetType>
class MMField2D : public MultiMat::BivariateMapType<DataType, BiSet> 
{
public:
  using BiVarSetType = BiSet;
  using BiVarMapType = MultiMat::BivariateMapType<DataType, BiVarSetType>;
  using ProductSetType = MultiMat::ProductSetType;
  using RelationSetType = MultiMat::RelationSetType;

  using Field2DType = MMField2D<DataType, BiVarSetType>;
  using SubFieldType = MMSubField2D<Field2DType>;
  using ConstSubFieldType = const MMSubField2D<const Field2DType>;

  /**
  * \brief Constructor
  */
  MMField2D() = delete;
  // Constructor with specified layouts
  //MMField2D(MultiMat& mm, DataLayout data_layout, SparsityLayout sparsity_layout,
  //  const std::string& arr_name = "unnamed", 
  //  const DataType* data_arr = nullptr,
  //  int stride = 1);
  // Constructor given the bivariate set
  MMField2D(MultiMat& mm, const BiSet*, 
    const std::string& arr_name = "unnamed", 
    const DataType* data_arr = nullptr,
    int stride = 1);
  
  /** Destructor **/
  ~MMField2D() {};
  /** Copy constructor (Deep copy). **/
  MMField2D(const MMField2D&);
  /** Assignment operator **/
  //MMField2D& operator=(const MMField2D&);


  //subfield (instead of SubMap) 
  SubFieldType getSubfield(SetPosition firstIdx)
  {
    return operator()(firstIdx);
  }
  SubFieldType operator() (SetPosition firstIdx)
  {
    const bool hasInd = submapIndicesHaveIndirection();
    return SubFieldType(this, firstIdx, hasInd);
  }
  const ConstSubFieldType operator() (SetPosition firstIdx) const
  {
    const bool hasInd = submapIndicesHaveIndirection();
    return ConstSubFieldType(this, firstIdx, hasInd);
  }

  //Mimic BivariateMap operator(i) and return slam submap
  SubMapType getSlamSubMap(SetPosition firstIdx)
  {
    return BiVarMapType::operator(firstIdx);
  }

  
  std::string getName() { return m_field_name; };

  MultiMat::IndexSet getSubfieldIndexingSet(int idx) {
    return m_mm->getSubfieldIndexingSet
    (idx, m_data_layout, m_sparsity_layout);
  }

  bool isDense() const { return m_sparsity_layout == SparsityLayout::DENSE; }
  bool isSparse() const { return m_sparsity_layout == SparsityLayout::SPARSE; }
  bool isCellDom() const { return m_data_layout == DataLayout::CELL_DOM; }
  bool isMatDom() const { return m_data_layout == DataLayout::MAT_DOM; }


private:
  MultiMat* m_mm;

  DataLayout m_data_layout;
  SparsityLayout m_sparsity_layout;

  std::string m_field_name;

};

///////////////////////////////////////////////////////////////////////////

// Constructor given layouts
//template<typename DataType, typename BiSet>
//inline MMField2D<DataType, BiSet>::MMField2D(MultiMat& mm, 
//    DataLayout data_layout, SparsityLayout sparsity_layout,
//    const std::string& arr_name, const DataType* data_arr, int stride) :
//  MMField2D(mm, mm.get_mapped_biSet(data_layout, sparsity_layout),
//     arr_name, data_arr, stride)
//{ }


// Constructor given the biset
template<typename DataType, typename BiSet>
inline MMField2D<DataType, BiSet>::MMField2D(MultiMat& mm, const BiSet *biset,
    const std::string& arr_name, const DataType* data_arr, int stride) :
  //call Bivariate map constructor
  MultiMat::BivariateMapType<DataType, BiSet>(biset, DataType(), stride),
  m_mm(&mm),
  m_field_name(arr_name)
{
  SLIC_ASSERT(stride > 0);

  if (data_arr != nullptr)
    this->copy(data_arr);

  if (biset == mm.get_mapped_biSet(DataLayout::CELL_DOM, SparsityLayout::DENSE)) {
    m_data_layout = DataLayout::CELL_DOM;
    m_sparsity_layout = SparsityLayout::DENSE;
  }
  else if (biset == mm.get_mapped_biSet(DataLayout::CELL_DOM, SparsityLayout::SPARSE)) {
    m_data_layout = DataLayout::CELL_DOM;
    m_sparsity_layout = SparsityLayout::SPARSE;
  }
  else if (biset == mm.get_mapped_biSet(DataLayout::MAT_DOM, SparsityLayout::DENSE)) {
    m_data_layout = DataLayout::MAT_DOM;
    m_sparsity_layout = SparsityLayout::DENSE;
  }
  else if (biset == mm.get_mapped_biSet(DataLayout::MAT_DOM, SparsityLayout::SPARSE)) {
    m_data_layout = DataLayout::MAT_DOM;
    m_sparsity_layout = SparsityLayout::SPARSE;
  }
  else {
    SLIC_ASSERT(false);
  }

  m_mm = &mm;
}

//Copy constructor
template<typename DataType, typename BiSet>
inline MMField2D<DataType, BiSet>::MMField2D(const MMField2D & o):
  MMField2D(*o.m_mm, o.set(), o.m_field_name,
    o.getMap()->data().data(), o.stride())
  { }


////////////////////////////////////////////////////

// class for field2d in specific format
template<typename DataType, DataLayout DataLayoutT, typename BiSet = MultiMat::BivariateSetType>
class MMField2DTemplated : public MMField2D<DataType, BiSet>
{ };


//maybe there's a way for me to not repeat code for same DataLayout below, 
//but I didn't figure it out.


// CellDOM Dense specialization
template<typename DataType>
class MMField2DTemplated<DataType, DataLayout::CELL_DOM, MultiMat::ProductSetType> : 
  public MMField2D<DataType, MultiMat::ProductSetType>
{
  using Field2DType = MMField2D<DataType, MultiMat::ProductSetType>;
public:
  MMField2DTemplated(MultiMat& mm, const std::string& arr_name = "unnamed",
    const DataType* data_arr = nullptr, int stride = 1):
    Field2DType(mm, mm.getDense2dFieldSet(DataLayout::CELL_DOM), 
      arr_name, data_arr, stride)
  {}


};


// CellDOM Sparse specialization
template<typename DataType>
class MMField2DTemplated<DataType, DataLayout::CELL_DOM, MultiMat::RelationSetType> :
  public MMField2D<DataType, MultiMat::RelationSetType>
{
  using Field2DType = MMField2D<DataType, MultiMat::RelationSetType>;
public:
  MMField2DTemplated(MultiMat& mm, const std::string& arr_name = "unnamed",
    const DataType* data_arr = nullptr, int stride = 1) :
    Field2DType(mm, mm.getSparse2dFieldSet(DataLayout::CELL_DOM), 
      arr_name, data_arr, stride)
  {}
};

// MatDom Dense specialization
template<typename DataType>
class MMField2DTemplated<DataType, DataLayout::MAT_DOM, MultiMat::ProductSetType> :
  public MMField2D<DataType, MultiMat::ProductSetType>
{
  using Field2DType = MMField2D<DataType, MultiMat::ProductSetType>;
public:
  MMField2DTemplated(MultiMat& mm,
    const std::string& arr_name = "unnamed",
    const DataType* data_arr = nullptr,
    int stride = 1) :
    Field2DType(mm, mm.getDense2dFieldSet(DataLayout::MAT_DOM),
      arr_name, data_arr, stride)
  {}
};


// MatDom Sparse specialization
template<typename DataType>
class MMField2DTemplated<DataType, DataLayout::MAT_DOM, MultiMat::RelationSetType> :
  public MMField2D<DataType, MultiMat::RelationSetType>
{
  using Field2DType = MMField2D<DataType, MultiMat::RelationSetType>;
public:
  MMField2DTemplated(MultiMat& mm,
    const std::string& arr_name = "unnamed",
    const DataType* data_arr = nullptr,
    int stride = 1) :
    Field2DType(mm, mm.getSparse2dFieldSet(DataLayout::MAT_DOM),
      arr_name, data_arr, stride)
  {}
};



} //end namespace multimat
} //end namespace axom

#endif