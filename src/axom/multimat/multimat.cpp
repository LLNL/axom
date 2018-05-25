#include "multimat/multimat.hpp" //Use this when committing
//#include "multimat.hpp"

#include <iostream>
#include <iterator>
#include <algorithm>

#include <cassert>


using namespace std;
using namespace axom::multimat;


MultiMatAbstractArray::MultiMatAbstractArray(std::string name, FieldMapping f)
	:m_arrayName(name),  m_fieldMapping(f)
{ }
MultiMatAbstractArray::~MultiMatAbstractArray() {}

MultiMat::MultiMat() {
	m_ncells = m_nmats = 0;
	m_layout = LAYOUT_CELL_DOM;
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
			if(vecarr[i*m_nmats+j])	{
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

	//Check the relation is set up correctly
	//for (int i = 0; i < m_ncells; i++) {
	//	auto a = m_cell2matRel[i];
	//	for (int j = 0; j < a.size(); j++){
	//		cout << i << " " << j << " " << a[j] << endl;;
	//	}
	//}

	cout << "indice total size: " << m_cell2matRel_indicesVec.size() << endl;
	cout << "cellmatrel total size: " << m_cell2matRel.totalSize() << endl;
	cout << "fromset size: " << m_cellSet.size() << endl;

	//Setup the cartesian set (a hack right now) of cell x mat

	m_cellMatNZSet = MappedRelationSetType(&m_cell2matRel);
	
}


MultiMatAbstractArray * axom::multimat::MultiMat::getFieldArray(std::string arr_name)
{
	for (auto arr : m_fieldArrayVec)
		if (arr->getName() == arr_name)
			return arr;

	return nullptr;
}

MultiMatAbstractArray* MultiMat::getFieldArray(int arr_idx)
{
	assert(arr_idx >= 0 && arr_idx < m_fieldArrayVec.size());
	return m_fieldArrayVec[arr_idx];
}


void axom::multimat::MultiMat::setLayout(DataLayout new_layout)
{
	if (new_layout == m_layout)
		return;

	//TODO
}

axom::multimat::DataLayout axom::multimat::MultiMat::getLayout()
{
	return m_layout;
}

std::string axom::multimat::MultiMat::getLayoutAsString()
{
	switch (m_layout) {
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
