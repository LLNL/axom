#ifndef MULTIMAT_H_
#define MULTIMAT_H_


#include "slam/RangeSet.hpp"
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

using SetType = slam::RangeSet;

enum FieldMapping { PER_CELL, PER_MAT, PER_CELL_MAT };
enum DataLayout { LAYOUT_CELL_DOM, LAYOUT_MAT_DOM };

class MultiMatAbstractArray //pure abstract class for returning array
{
public:
	MultiMatAbstractArray(std::string name, FieldMapping f);
	virtual ~MultiMatAbstractArray() = 0;
	void setName(std::string arr_name) { m_arrayName = arr_name; }
	std::string getName() { return m_arrayName; }
	FieldMapping getFieldMapping() { return m_fieldMapping; };

protected:
	FieldMapping m_fieldMapping;
	std::string m_arrayName; //name of the array
};


template<class T>
class MultiMatArray : public MultiMatAbstractArray
{
public:
	using MapType = axom::slam::Map<T>;

	MultiMatArray(std::string name, SetType* s, FieldMapping f);
	MultiMatArray(std::string name, SetType* s, FieldMapping f, T* arr);
	~MultiMatArray() {}

	T getValue(int c);
	T getValue(int c, int m);
	
	void setValue(int c, T val);
	void setValue(int c, int m, T val);

	const MapType getMap() { return m_map; }

	class iterator : public std::iterator<std::input_iterator_tag, T>
	{
	MapType* m_mapPtr;
	int m_i;
	public:
		iterator(MapType* m) : m_mapPtr(m), m_i(0) {}
		iterator(MapType* m, int n) : m_mapPtr(m), m_i(n) {}
		iterator& operator++() { ++m_i; return *this; }
		iterator operator++(int) { iterator tmp(*this); operator++(); return tmp; }
		bool operator==(const iterator& rhs) { return m_mapPtr == rhs.m_mapPtr && m_i == rhs.m_i; }
		bool operator!=(const iterator& rhs) { return m_mapPtr != rhs.m_mapPtr || m_i != rhs.m_i; }
		T operator*() { return m_mapPtr->operator[](m_i); }
	};
	iterator begin() { return iterator(&m_map); }
	iterator end() { return iterator(&m_map, m_map.size()); }

//private:
	MapType m_map;
	SetType* m_set; //storing the set in order to index into the elemXmat map. 
	//Eventually the map should be changed to allow access for this
};


class MultiMat
{
public:
	// SLAM Set type definitions
	//using SetType     = slam::RangeSet; //moved out
	using SetPosType  = SetType::PositionType;
	using SetElemType = SetType::ElementType;
	// SLAM Relation typedef (static variable relation)
	using STLIndirection             = policies::STLVectorIndirection<SetPosType, SetPosType>;
	using VariableCardinality        = policies::VariableCardinality<SetPosType, STLIndirection>;
	using StaticVariableRelationType = slam::StaticRelation<VariableCardinality, STLIndirection,
		slam::RangeSet, slam::RangeSet>;
	
	// SLAM MappedRelationSet for the set of non-zero cell to mat variables
	using MappedRelationSetType = slam::MappedRelationSet<StaticVariableRelationType>;

#define MM_CAST_TO( t, a) dynamic_cast<MultiMatArray<t>*>(a)

	MultiMat();

	//Set-up functions
	void setNumberOfMat(int); 
	void setNumberOfCell(int);
	void setNumberOfField(int);
	void setCellMatRel(std::vector<bool>&); //relation information

	void build(); //after all the set-up functions are called, build the internal structures

	//SLAM set-up functions
	template<class T>
	MultiMatArray<T>* newFieldArray(std::string arr_name, FieldMapping); //get the array
	template<class T>
	MultiMatArray<T>* newFieldArray(std::string arr_name, FieldMapping, T* arr); //get the array
	MultiMatAbstractArray* getFieldArray(std::string arr_name); //get the array
	MultiMatAbstractArray* getFieldArray(int arr_idx); //get the array
	
	//accessing functions
	template<class T>
	T getFieldValue(std::string field_name, int cell_i, int mat_i);
	template<class T>
	T getFieldValue(std::string field_name, int i);

	//Data modification functions
	//...

	//Layout modification functions
	void setLayoutToCellDominant() { setLayout(LAYOUT_CELL_DOM); }
	void setLayoutToMaterialDominant() { setLayout(LAYOUT_MAT_DOM); }
	void setLayout(DataLayout);
	DataLayout getLayout();
	std::string getLayoutAsString();

	bool isValid();
	void printSelf();
	
public: //private:
	int m_nmats, m_ncells;
	bool m_modifiedSinceLastBuild;
	DataLayout m_layout;

	//slam variables
	SetType m_cellSet;
	SetType m_matSet;
	StaticVariableRelationType m_cell2matRel;
	std::vector<SetPosType> m_cell2matRel_beginsVec; //to store the cell2mat relation
	std::vector<SetPosType> m_cell2matRel_indicesVec; //to store the cell2mat relation
	MappedRelationSetType m_cellMatNZSet; // set of non-zero entries in the cellXmat matrix

	std::vector<MultiMatAbstractArray*> m_fieldArrayVec; //list of MM-Arrays

};


//------------ MultiMatArray Template function definitions --------------//

template<class T>
inline MultiMatArray<T>::MultiMatArray(std::string name, SetType * s, FieldMapping f)
	:MultiMatAbstractArray(name, f), m_map(s), m_set(s)
{
}

template<class T>
inline MultiMatArray<T>::MultiMatArray(std::string name, SetType * s, FieldMapping f, T * data_arr)
	: MultiMatAbstractArray(name, f), m_map(s), m_set(s)
{
	//copy the array
	for (int i = 0; i < m_map.size(); i++) {
		setValue(i, data_arr[i]);
	}
}

template<class T>
inline T MultiMatArray<T>::getValue(int c)
{
	assert(m_fieldMapping == PER_CELL || m_fieldMapping == PER_MAT);
	return m_map[c];
}

template<class T>
inline T MultiMatArray<T>::getValue(int c, int m)
{
	assert(m_fieldMapping == PER_CELL_MAT);
	auto set_ptr = dynamic_cast<MultiMat::MappedRelationSetType*>(m_set);
	auto from_set_size = set_ptr->fromSetSize();
	auto to_set_size = set_ptr->toSetSize();
	assert(c >= 0 && c < from_set_size);
	assert(m >= 0 && m < to_set_size);
	
	auto i = set_ptr->at(c, m);
	if (i == -1) //this cell does not contain this material
		return 0;
	return m_map[i];
}

template<class T>
inline void MultiMatArray<T>::setValue(int c, T val)
{
	assert(m_fieldMapping == PER_CELL || m_fieldMapping == PER_MAT);
	assert(c >= 0 && c < m_map.size()); 
	m_map[c] = val;
}

template<class T>
inline void MultiMatArray<T>::setValue(int c, int m, T val)
{
	assert(m_fieldMapping == PER_CELL_MAT);
	auto set_ptr = dynamic_cast<MultiMat::MappedRelationSetType*>(m_set);
	auto i = set_ptr->at(c, m);
	m_map[i] = val;
}


//--------------- MultiMat template function definitions -----------------//

template<class T>
MultiMatArray<T>* axom::multimat::MultiMat::newFieldArray(std::string arr_name, FieldMapping arr_mapping)
{
	SetType* map_set;
	switch (arr_mapping)
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
	auto new_arr = new MultiMatArray<T>(arr_name, map_set, arr_mapping);
	new_arr->setName(arr_name);
	m_fieldArrayVec.push_back(new_arr);

	return new_arr;
}

template<class T>
MultiMatArray<T>* axom::multimat::MultiMat::newFieldArray(std::string arr_name, FieldMapping arr_mapping, T* data_arr)
{
	auto new_arr = newFieldArray<T>(arr_name, arr_mapping);

	//copy data
	for (int i = 0; i < new_arr->m_map.size(); i++) {
		new_arr->m_map[i] = data_arr[i];
	}

	return new_arr;
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


} //end namespace multimat
} //end namespace axom


#endif