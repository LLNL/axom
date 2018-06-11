#include "multimat/multimat.hpp" 

#include <iostream>
#include <iterator>
#include <algorithm>

#include <cassert>


using namespace std;
using namespace axom::multimat;


MultiMat::MultiMat() {
  m_ncells = m_nmats = 0;
  m_dataLayout = LAYOUT_CELL_DOM;
  m_sparcityLayout = LAYOUT_SPARSE;
}

MultiMat::MultiMat(DataLayout d, SparcityLayout s): MultiMat() {
  m_dataLayout = d;
  m_sparcityLayout = s;
}

axom::multimat::MultiMat::~MultiMat()
{
  for (auto mapPtr : m_mapVec) {
    delete mapPtr;
  }
}

void MultiMat::setNumberOfMat(int n){
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

  assert(vecarr.size() == m_ncells * m_nmats); //This should be a dense matrix

  assert(m_dataLayout == LAYOUT_CELL_DOM); //for now assumes cell dominant

  //Set-up the cell/mat relation
  m_cell2matRel_beginsVec.resize(m_cellSet.size() + 1, -1);

  SetPosType curIdx = SetPosType();
  for (SetPosType i = 0; i < m_ncells; ++i)
  {
    m_cell2matRel_beginsVec[i] = curIdx;
    for (SetPosType j = 0; j < m_nmats; ++j)
    {
      if (vecarr[i*m_nmats + j]) {
        m_cell2matRel_indicesVec.push_back(j);
        ++curIdx;
      }
    }
  }
  m_cell2matRel_beginsVec[m_ncells] = curIdx;

  m_cell2matRel = StaticVariableRelationType(&m_cellSet, &m_matSet);
  m_cell2matRel.bindBeginOffsets(m_cellSet.size(), &m_cell2matRel_beginsVec);
  m_cell2matRel.bindIndices(m_cell2matRel_indicesVec.size(), &m_cell2matRel_indicesVec);

  assert(m_cell2matRel.isValid());

  cout << "indice total size: " << m_cell2matRel_indicesVec.size() << endl;
  cout << "cellmatrel total size: " << m_cell2matRel.totalSize() << endl;
  cout << "fromset size: " << m_cellSet.size() << endl;

  //Set-up both dense and sparse sets, since they don't take any memory...

  //a set of mapped relation
  m_cellMatNZSet = MappedRelationSetType(&m_cell2matRel);

  // a cartesian set of cell x mat
  if(m_dataLayout == LAYOUT_CELL_DOM)
    m_cellMatProdSet = ProductSetType(&m_cellSet, &m_matSet);
  else if(m_dataLayout==LAYOUT_MAT_DOM)
    m_cellMatProdSet = ProductSetType(&m_matSet, &m_cellSet);
  else assert(false);

}

MultiMat::IdSubSet MultiMat::getMatInCell(int c)
{
  if (m_dataLayout == LAYOUT_SPARSE) {
    return m_cell2matRel[c]; //returns a RelationSet / OrderedSet with STLindirection
  }
  //else if (m_dataLayout == LAYOUT_DENSE)
  //  return m_cell2matRel[c]; //since the relation is currently only stored sparse, return the same thing.
  else assert(false);
}

MultiMat::IndexSet axom::multimat::MultiMat::getIndexingSetOfCell(int c)
{
  if (m_sparcityLayout == LAYOUT_SPARSE) {
    assert(0 <= c && c < m_ncells);
    int start_idx = m_cell2matRel_beginsVec[c];
    int end_idx = m_cell2matRel_beginsVec[c + 1];
    return RangeSetType::SetBuilder().range(start_idx, end_idx);
  }
  else if (m_sparcityLayout == LAYOUT_DENSE) {
    //return m_cellMatProdSet.getRow(c);
    int size1 = m_cellMatProdSet.set1Size();
    int size2 = m_cellMatProdSet.set2Size();
    return RangeSetType::SetBuilder().range(c*size2, (c + 1)*size2 - 1);
  }
  else assert(false);
}

void axom::multimat::MultiMat::convertLayout(DataLayout new_layout, SparcityLayout new_sparcity)
{
  if (new_layout == m_dataLayout && new_sparcity == m_sparcityLayout)
    return;

  //assumes cell dom stays cell dom
  assert(new_layout == LAYOUT_CELL_DOM);

  if (m_sparcityLayout == LAYOUT_DENSE && new_sparcity == LAYOUT_SPARSE)
  { //convert from dense to sparse

    //go through each field, for every matXcell field, create a new map of sparse mat
    for (int i = 0; i < m_fieldMappingVec.size(); i++)
    {
      if (m_fieldMappingVec[i] != PER_CELL_MAT) 
        continue;

      //TODO
    }
  }
  else {
    assert(false);
    //TODO
  }
}

axom::multimat::DataLayout axom::multimat::MultiMat::getDataLayout()
{
  return m_dataLayout;
}

std::string axom::multimat::MultiMat::getLayoutAsString()
{
  switch (m_dataLayout) {
  case LAYOUT_CELL_DOM:
    return "LAYOUT_CELL_DOM";
  case LAYOUT_MAT_DOM:
    return "LAYOUT_MAT_DOM";
  default:
    assert(false);
    return "";
  }
}

SparcityLayout axom::multimat::MultiMat::getSparcityLayout()
{
  return m_sparcityLayout;
}

void axom::multimat::MultiMat::printSelf()
{
  printf("Multimat Object\n");
  printf("Number of materials: %d\n", m_nmats);
  printf("Number of cells:     %d\n", m_ncells);

  printf("\nFields:\n");
  for (int i = 0; i < m_mapVec.size(); i++)
  {
    printf("Field %d - %s\n", i, m_arrNameVec[i].c_str());
    printf("  Mapping type: %d\n", m_fieldMappingVec[i]);
  }
}

bool axom::multimat::MultiMat::isValid()
{
  //TODO
  return true;
}

axom::multimat::MultiMat::SetType* axom::multimat::MultiMat::get_mapped_set(FieldMapping fm)
{
  SetType* map_set;
  switch (fm)
  {
  case PER_CELL:
    map_set = &m_cellSet;
    break;
  case PER_MAT:
    map_set = &m_matSet;
    break;
  case PER_CELL_MAT:
    if (m_sparcityLayout == LAYOUT_SPARSE)
      map_set = &m_cellMatNZSet;
    else if (m_sparcityLayout == LAYOUT_DENSE)
      map_set = &m_cellMatProdSet;
    else assert(false);
    break;
  default:
    assert(false);
    return nullptr;
  }
  return map_set;
}


