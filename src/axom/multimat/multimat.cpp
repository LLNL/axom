#include "multimat/multimat.hpp" 

#include <iostream>
#include <iterator>
#include <algorithm>

#include "slic/slic.hpp"
#include <cassert>
#include "multimat.hpp"


using namespace std;
using namespace axom::multimat;


MultiMat::MultiMat():
  m_cellMatRel(nullptr),
  m_cellMatNZSet(nullptr),
  m_cellMatProdSet(nullptr)
{
  m_ncells = m_nmats = 0;
  m_dataLayout = DataLayout::CELL_CENTRIC;
  m_sparcityLayout = SparcityLayout::SPARSE;
}

MultiMat::MultiMat(DataLayout d, SparcityLayout s): MultiMat() {
  m_dataLayout = d;
  m_sparcityLayout = s;
}

MultiMat::~MultiMat()
{
  for (auto mapPtr : m_mapVec) {
    delete mapPtr;
  }
  delete m_cellMatProdSet;
  delete m_cellMatNZSet;
  delete m_cellMatRel;
}

template<typename T>
MultiMat::MapBaseType* MultiMat::helperfun_copyField(const MultiMat& mm, int map_i)
{
  MapBaseType* other_map_ptr = mm.m_mapVec[map_i];
  if (mm.getFieldMapping(map_i) == FieldMapping::PER_CELL_MAT)
  {
    BivariateSetType* biSetPtr = dynamic_cast<BivariateSetType*>(get_mapped_set(getFieldMapping(map_i)));

    MultiMat::Field2D<T>* typed_ptr = dynamic_cast<MultiMat::Field2D<T>*>(other_map_ptr);
    MultiMat::Field2D<T>* new_ptr = new MultiMat::Field2D<T>(biSetPtr, T(), typed_ptr->stride());
    new_ptr->copy(typed_ptr->getMap()->data().data());
    return new_ptr;
  }
  else
  {
    SetType* setPtr = get_mapped_set(getFieldMapping(map_i));
    MultiMat::Field1D<T>* typed_ptr = dynamic_cast<MultiMat::Field1D<T>*>(other_map_ptr);
    MultiMat::Field1D<T>* new_ptr = new MultiMat::Field1D<T>(setPtr, T(), typed_ptr->stride());
    new_ptr->copy(*typed_ptr);
    return new_ptr;
  }
}


MultiMat::MultiMat(const MultiMat& other) :
  m_ncells(other.m_ncells),
  m_nmats(other.m_nmats),
  m_dataLayout( other.m_dataLayout),
  m_sparcityLayout( other.m_sparcityLayout),
  m_matSet(0, other.m_nmats),
  m_cellSet(0,other.m_ncells),
  m_cellMatRel_beginsVec(other.m_cellMatRel_beginsVec),
  m_cellMatRel_indicesVec(other.m_cellMatRel_indicesVec),
  m_cellMatRel(nullptr),
  m_cellMatNZSet(nullptr),
  m_cellMatProdSet(nullptr),
  m_arrNameVec(other.m_arrNameVec),
  m_fieldMappingVec(other.m_fieldMappingVec),
  m_dataTypeVec(other.m_dataTypeVec)
{
  RangeSetType& set1 = (m_dataLayout == DataLayout::CELL_CENTRIC ? m_cellSet : m_matSet);
  RangeSetType& set2 = (m_dataLayout == DataLayout::CELL_CENTRIC ? m_matSet : m_cellSet);
  m_cellMatRel = new StaticVariableRelationType(&set1, &set2);
  m_cellMatRel->bindBeginOffsets(set1.size(), &m_cellMatRel_beginsVec);
  m_cellMatRel->bindIndices(m_cellMatRel_indicesVec.size(), &m_cellMatRel_indicesVec);
  m_cellMatNZSet = new RelationSetType(m_cellMatRel);
  m_cellMatProdSet = new ProductSetType(&set1, &set2);


  for (unsigned int map_i = 0; map_i < other.m_mapVec.size(); ++map_i)
  {
    MapBaseType* new_map_ptr = nullptr;
    FieldMapping fm = m_fieldMappingVec[map_i];
    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble) {
      new_map_ptr = helperfun_copyField<double>(other, map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat) {
      new_map_ptr = helperfun_copyField<float>(other, map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt) {
      new_map_ptr = helperfun_copyField<int>(other, map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar) {
      new_map_ptr = helperfun_copyField<unsigned char>(other, map_i);
    }
    else assert(false); //TODO
    
    SLIC_ASSERT(new_map_ptr != nullptr);
    m_mapVec.push_back(new_map_ptr);
  }
}

MultiMat& MultiMat::operator=(const MultiMat& other)
{
  if (this == &other) return *this;

  SLIC_ASSERT(false);

  return *this;
}


void MultiMat::setNumberOfMat(int n)
{
  assert(n > 0);
  m_nmats = n;

  m_matSet = RangeSetType(0, m_nmats);
  assert(m_matSet.isValid());
}

void MultiMat::setNumberOfCell(int c)
{
  assert(c > 0);
  m_ncells = c;

  m_cellSet = RangeSetType(0, m_ncells);
  assert(m_cellSet.isValid());
}

void MultiMat::setCellMatRel(vector<bool>& vecarr)
{
  //Setup the SLAM cell to mat relation
  //This step is necessary if the volfrac field is sparse

  SLIC_ASSERT(vecarr.size() == m_ncells * m_nmats); //This should be a dense matrix
  SLIC_ASSERT(m_cellMatRel == nullptr); //cellmatRel has not been set before

  RangeSetType& set1 = (m_dataLayout == DataLayout::CELL_CENTRIC ? m_cellSet : m_matSet);
  RangeSetType& set2 = (m_dataLayout == DataLayout::CELL_CENTRIC ? m_matSet : m_cellSet);

  //count the non-zeros
  int nz_count = 0;
  for (bool b : vecarr)
    nz_count += b;

  //Set-up the cell/mat relation
  m_cellMatRel_beginsVec.resize(set1.size() + 1, -1);
  m_cellMatRel_indicesVec.resize(nz_count);

  SetPosType curIdx = SetPosType();
  for (SetPosType i = 0; i < set1.size(); ++i)
  {
    m_cellMatRel_beginsVec[i] = curIdx;
    for (SetPosType j = 0; j < set2.size(); ++j)
    {
      if (vecarr[i*set2.size() + j]) {
        m_cellMatRel_indicesVec[curIdx] = j;
        ++curIdx;
      }
    }
  }
  m_cellMatRel_beginsVec[set1.size()] = curIdx;

  m_cellMatRel = new StaticVariableRelationType(&set1, &set2);
  m_cellMatRel->bindBeginOffsets(set1.size(), &m_cellMatRel_beginsVec);
  m_cellMatRel->bindIndices(m_cellMatRel_indicesVec.size(), &m_cellMatRel_indicesVec);

  assert(m_cellMatRel->isValid());
  
  cout << "indice total size: " << m_cellMatRel_indicesVec.size() << endl;
  cout << "cellmatrel total size: " << m_cellMatRel->totalSize() << endl;
  cout << "fromset size: " << set1.size() << endl;

  //Set-up both dense and sparse sets, since they don't take any extra memory...

  //a set of mapped relation
  m_cellMatNZSet = new RelationSetType(m_cellMatRel);
  
  // a cartesian set of cell x mat
  m_cellMatProdSet = new ProductSetType(&set1, &set2);
  
  
  //Create a field for VolFrac as the 0th field
  m_mapVec.push_back(nullptr);
  m_arrNameVec.push_back("Volfrac");
  m_fieldMappingVec.push_back(FieldMapping::PER_CELL_MAT);
  m_dataTypeVec.push_back(DataTypeSupported::TypeDouble);
  assert(m_mapVec.size() == 1);
  assert(m_arrNameVec.size() == 1);
  assert(m_fieldMappingVec.size() == 1);
  assert(m_dataTypeVec.size() == 1);

}


int MultiMat::setVolfracField(double* arr)
{
  //Assumes this is to CellxMat mapping, named "Volfrac", and is stride 1.
  int arr_i = addFieldArray_impl<double>("Volfrac", FieldMapping::PER_CELL_MAT, arr, 1);

  //checks volfrac values sums to 1
  auto& map = *dynamic_cast<Field2D<double>*>(m_mapVec[arr_i]);
  double tol = 10e-9;
  if (m_dataLayout == DataLayout::CELL_CENTRIC) 
  {
    for (int i = 0; i < map.firstSetSize(); ++i)
    {
      double sum = 0.0;
      for (auto iter = map.begin(i); iter != map.end(i); ++iter)
      {
        sum += iter.value();
      }

      SLIC_ASSERT(abs(sum - 1.0) < tol);
    }
  }
  else //material centric layout
  {
    std::vector<double> sum_vec(m_cellSet.size(), 0);
    for (int i = 0; i < map.firstSetSize(); ++i)
    {
      for (auto iter = map.begin(i); iter != map.end(i); ++iter)
      {
        sum_vec[iter.index()] += iter.value();
      }
    }

    for(int i=0; i<sum_vec.size(); ++i)
      SLIC_ASSERT(abs(sum_vec[i] - 1.0) < tol);
  }

  //move the data to the first one (index 0) in the list
  std::iter_swap(m_mapVec.begin(), m_mapVec.begin() + arr_i);
  std::iter_swap(m_dataTypeVec.begin(), m_dataTypeVec.begin() + arr_i);

  //remove the new entry...
  int nfield = m_mapVec.size() - 1;
  m_mapVec.resize(nfield); //TODO if delete is needed
  m_fieldMappingVec.resize(nfield);
  m_arrNameVec.resize(nfield);
  m_dataTypeVec.resize(nfield);

  return 0;
}


MultiMat::Field2D<double>& MultiMat::getVolfracField()
{
  return *dynamic_cast<Field2D<double>*>(m_mapVec[0]);
}


int MultiMat::getFieldIdx(std::string field_name)
{
  for (unsigned int i = 0; i < m_arrNameVec.size(); i++)
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

  if (m_sparcityLayout == SparcityLayout::SPARSE) {
    return (*m_cellMatRel)[c]; //returns a RelationSet / OrderedSet with STLindirection
  }
  else if (m_sparcityLayout == SparcityLayout::DENSE)
    return (*m_cellMatRel)[c]; //since the relation is currently only stored sparse, return the same thing.
  else assert(false);
}


MultiMat::IndexSet MultiMat::getIndexingSetOfCell(int c)
{
  assert(m_dataLayout == DataLayout::CELL_CENTRIC);

  if (m_sparcityLayout == SparcityLayout::SPARSE) {
    assert(0 <= c && c < m_ncells);
    int start_idx = m_cellMatRel_beginsVec[c];
    int end_idx = m_cellMatRel_beginsVec[c + 1];
    return RangeSetType::SetBuilder().range(start_idx, end_idx);
  }
  else if (m_sparcityLayout == SparcityLayout::DENSE) {
    //return m_cellMatProdSet.getRow(c);
    int size2 = m_cellMatProdSet->secondSetSize();
    return RangeSetType::SetBuilder().range(c*size2, (c + 1)*size2 - 1);
  }
  else assert(false);
}

void MultiMat::convertLayoutToCellDominant()
{
  if (m_dataLayout == DataLayout::CELL_CENTRIC) return;
  transposeData();
}

void MultiMat::convertLayoutToMaterialDominant() 
{ 
  if (m_dataLayout == DataLayout::MAT_CENTRIC) return;
  transposeData();
}

void MultiMat::transposeData()
{
  //create the new relation to pass into helper...
  //TODO
  
  for (unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    //no conversion needed unless the field is PER_CELL_MAT
    if (m_fieldMappingVec[map_i] != FieldMapping::PER_CELL_MAT)
      continue;

    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble) {
      transposeData_helper<double>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat) {
      transposeData_helper<float>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt) {
      transposeData_helper<int>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar) {
      transposeData_helper<unsigned char>(map_i);
    }
    else assert(false); //TODO
  }

  if(m_dataLayout == DataLayout::MAT_CENTRIC)
    m_dataLayout = DataLayout::CELL_CENTRIC;
  else
    m_dataLayout = DataLayout::MAT_CENTRIC;
}

void MultiMat::convertLayoutToSparse() 
{ 
  if (m_sparcityLayout == SparcityLayout::SPARSE) return;
  
  for (unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    //no conversion needed unless the field is PER_CELL_MAT
    if (m_fieldMappingVec[map_i] != FieldMapping::PER_CELL_MAT)
      continue;

    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble) {
      convertToSparse_helper<double>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat) {
      convertToSparse_helper<float>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt) {
      convertToSparse_helper<int>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar) {
      convertToSparse_helper<unsigned char>(map_i);
    }
    else assert(false); //TODO
  }
  m_sparcityLayout = SparcityLayout::SPARSE;
}


void MultiMat::convertLayoutToDense() 
{ 
  if(m_sparcityLayout == SparcityLayout::DENSE) return;

  for (unsigned int map_i = 0; map_i < m_fieldMappingVec.size(); map_i++)
  {
    //no conversion needed unless the field is PER_CELL_MAT
    if (m_fieldMappingVec[map_i] != FieldMapping::PER_CELL_MAT)
      continue;

    if (m_dataTypeVec[map_i] == DataTypeSupported::TypeDouble) {
      convertToDense_helper<double>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeFloat) {
      convertToDense_helper<float>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeInt) {
      convertToDense_helper<int>(map_i);
    }
    else if (m_dataTypeVec[map_i] == DataTypeSupported::TypeUnsignChar) {
      convertToDense_helper<unsigned char>(map_i);
    }
    else assert(false); //TODO
  }
  m_sparcityLayout = SparcityLayout::DENSE;
}


void MultiMat::convertLayout(DataLayout new_layout, SparcityLayout new_sparcity)
{
  if (new_layout == m_dataLayout && new_sparcity == m_sparcityLayout)
    return;

  //sparse/dense conversion
  if (m_sparcityLayout == SparcityLayout::DENSE && new_sparcity == SparcityLayout::SPARSE)
  {
    convertLayoutToSparse();
  }
  else if(m_sparcityLayout == SparcityLayout::SPARSE && new_sparcity == SparcityLayout::DENSE)
  {
    convertLayoutToDense();
  }

  //cell/mat centric conversion
  if (m_dataLayout == DataLayout::CELL_CENTRIC && new_layout == DataLayout::MAT_CENTRIC)
  {
    convertLayoutToMaterialDominant();
  }
  else if(m_dataLayout == DataLayout::MAT_CENTRIC && new_layout == DataLayout::CELL_CENTRIC)
  {
    convertLayoutToCellDominant();
  }
}


std::string MultiMat::getDataLayoutAsString() const
{
  switch (m_dataLayout) {
  case DataLayout::CELL_CENTRIC:
    return "Cell-Centric";
  case DataLayout::MAT_CENTRIC:
    return "Material-Centric";
  default:
    assert(false);
  }
}

std::string MultiMat::getSparcityLayoutAsString() const
{
  switch (m_sparcityLayout) {
  case SparcityLayout::DENSE:
    return "Dense";
  case SparcityLayout::SPARSE:
    return "Sparse";
  default:
    assert(false);
  }
}


void MultiMat::print() const
{
  printf("Multimat Object\n");
  printf("Number of materials: %d\n", m_nmats);
  printf("Number of cells:     %d\n", m_ncells);
  printf("Data layout: ");
  if (m_dataLayout == DataLayout::CELL_CENTRIC)
    printf("Cell-centric\n");
  else
    printf("Material-centric\n");

  printf("Sparcity layout: ");
  if (m_sparcityLayout == SparcityLayout::DENSE)
    printf("Dense\n");
  else
    printf("Sparse\n");

  printf("\n%d Fields:\n", m_mapVec.size());
  for (unsigned int i = 0; i < m_mapVec.size(); i++)
  {
    printf("Field %d - %s\n", i, m_arrNameVec[i].c_str());
    printf("  Mapping per ");
    switch (m_fieldMappingVec[i]) {
    case FieldMapping::PER_CELL: printf("cell"); break;
    case FieldMapping::PER_MAT: printf("material"); break;
    case FieldMapping::PER_CELL_MAT: printf("cellXmaterial"); break;
    }
  }
  printf("\n\n");
  
}


bool MultiMat::isValid(bool verboseOutput) const
{
  if (verboseOutput) print();

  //make sure there is a volfrac field filled out 
  //that matches the CellMatRelation
  //TODO

  return true;
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
    assert(false);
    return nullptr;
  }
  return set_ptr;
}

MultiMat::BivariateSetType* MultiMat::get_mapped_biSet()
{
  BivariateSetType* set_ptr = nullptr;
  if (m_sparcityLayout == SparcityLayout::SPARSE)
    set_ptr = m_cellMatNZSet;
  else if (m_sparcityLayout == SparcityLayout::DENSE)
    set_ptr = m_cellMatProdSet;
  
  SLIC_ASSERT(set_ptr != nullptr);
  return set_ptr;
}

