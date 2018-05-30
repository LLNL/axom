#include "multimat/multimat.hpp" 

#include <iostream>
#include <iterator>
#include <algorithm>

#include <cassert>


using namespace std;
using namespace axom::multimat;


MultiMatArray::MultiMatArray(MultiMat* m, std::string name, FieldMapping f)
  :m_arrayName(name),  m_fieldMapping(f), mm(m)
{ }

//axom::multimat::MultiMatArray::MultiMatArray(MultiMat * m, std::string name, SetType * s, FieldMapping f)
//{
//  assert(false);
//
//}

MultiMat::MultiMat() {
  m_ncells = m_nmats = 0;
  m_dataLayout = LAYOUT_CELL_DOM;
  m_modifiedSinceLastBuild = false;
}

void MultiMat::setNumberOfMat(int n){
  assert(n > 0);
  m_nmats = n;
  m_modifiedSinceLastBuild = true;

  m_matSet = SetType(0, m_nmats);
  assert(m_matSet.isValid());
}

void MultiMat::setNumberOfCell(int c)
{
  assert(c > 0);
  m_ncells = c;
  m_modifiedSinceLastBuild = true;

  m_cellSet = SetType(0, m_ncells);
  assert(m_cellSet.isValid());
}

void MultiMat::setCellMatRel(vector<bool>& vecarr)
{
  
  m_modifiedSinceLastBuild = true;

  //Setup the SLAM cell to mat relation
  assert(vecarr.size() == m_ncells*m_nmats);
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

  //Setup the cartesian set (a hack right now) of cell x mat

  m_cellMatNZSet = MappedRelationSetType(&m_cell2matRel);
  
}

MultiMat::RelationSet MultiMat::getMatInCell(int c)
{
  return m_cell2matRel[c];
}

MultiMatArray * axom::multimat::MultiMat::getFieldArray(std::string arr_name)
{
  for (auto arr : m_fieldArrayVec)
    if (arr->getName() == arr_name)
      return arr;

  return nullptr;
}

MultiMatArray* MultiMat::getFieldArray(int arr_idx)
{
  assert(arr_idx >= 0 && arr_idx < m_fieldArrayVec.size());
  return m_fieldArrayVec[arr_idx];
}


void axom::multimat::MultiMat::convertLayout(DataLayout new_layout)
{
  if (new_layout == m_dataLayout)
    return;

  //TODO
}

axom::multimat::DataLayout axom::multimat::MultiMat::getLayout()
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

void axom::multimat::MultiMat::printSelf()
{
  printf("Multimat Object\n");
  printf("Number of materials: %d\n", m_nmats);
  printf("Number of cells:     %d\n", m_ncells);

  printf("\nFields:\n");
  for (int i = 0; i < m_fieldArrayVec.size(); i++)
  {
    printf("Field %d - %s\n", i, m_fieldArrayVec[i]->getName().c_str());
    printf("  Type: %d\n", m_fieldArrayVec[i]->getFieldMapping());
  }

}

bool axom::multimat::MultiMat::isValid()
{

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
    map_set = &m_cellMatNZSet;
    break;
  default:
    return nullptr;
  }
  return map_set;
}


