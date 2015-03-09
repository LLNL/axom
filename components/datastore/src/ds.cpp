// Datastore
//  TODO: namespace
//
//   Groups do not know their own name.
//   DataObject does not have a single parent but a list of groups (not group names).
//   Datastore have pool of groups similar to objects with a gid.
//     helpful for Fortran, maybe helpful to have a constant id for restart (instead of address)

// C++ library headers
#include <string>
#include <stack>
#include <vector>
#include <map>
#include <set>

class DataStore;
class DataGroup;

// typedefs, in case we need to change things in the future
/*!
 * \brief Universal id of DataObjects.
 */
typedef size_t IDType; 


// --------------------------------------------------------------------------------
// DataObject.hpp ----------------------------------------

/**
 * \class DataObject
 *
 * \brief Class to access data.
 * Each DataObject has a universal id associated with it.
 * A DataObject may belong to any number of DataGroups.
 * A DataObject has no name except the name given to it by each DataGroup.
 */
class DataObject
{
public:
  //  typedef std::map< std::string, DataGroup* > GroupMap;
  /*!
   * \brief Set of DataGroup pointers.
   */
  typedef std::set< DataGroup* > GroupSet;

private:
  // universal identification - unique within a DataStore
  // XXX do not move Objects between datastores.
  IDType m_uid;
  //  GroupMap mDataGroups;  // parent groups: name->Group pointer
  GroupSet m_GroupSet;

public:
  //--- sample attribute, would need a more generic system
  bool m_dump;


  /*!
   * \brief Constructor
   * XXX User should not call directly (friend of DataStore?)
   */
  DataObject(IDType uid) :
    m_uid(uid),
    m_dump(false)
  {};

  /*!
   * \brief Return the univeral id for this DataObject.
   */
  IDType GetUID() { return m_uid; }

  /*!
   * @param grp Pointer to DataGroup to associate with this DataObject.
   * \brief Associate a DataGroup with this DataObject.
   */
  void AttachGroup( DataGroup *grp ) { m_GroupSet.insert(grp); }

  /*!
   * @param grp Pointer to DataGroup to disassociate with this DataObject.
   * \brief Disassociate a DataGroup with this DataObject.
   */
  void DetachGroup( DataGroup *grp ) { m_GroupSet.erase(grp); }

  /*!
   * @param grp Pointer to DataGroup to test for membership.
   * \brief return true is grp is attached.
   */
  bool IsAttachedGroup( DataGroup *grp )
  {
    GroupSet::iterator it = m_GroupSet.find(grp);
    if (it == m_GroupSet.end())
      return false;
    else
      return true;
  }

  /*!
   * \brief Return DataGroups attached to this DataObject.
   */
  GroupSet *GetDataGroups() { return &m_GroupSet; }
};

// DataGroup.hpp ----------------------------------------

/**
 * \class DataGroup
 *
 * \brief Class to access collections of DataObject.
 *  The DataGroup will name each DataObject as it is added.
 */
class DataGroup
{
public:
  /*!
   * \brief vector of DataObject pointers.
   */
  typedef std::vector< DataObject* > dataArrayType;

  /*!
   * \brief map of name to index of DataObject within this DataGroup.
   */
  typedef std::map< std::string, IDType > lookupType;

  /*!
   * \brief map of name to DataGroup pointer.
   */
  typedef std::map< std::string, DataGroup* > lookupGroup;

private:
  DataGroup *m_parent;
  DataStore *m_datastore;
  dataArrayType m_DataObjects;  // DataObjects by index
  lookupType m_DataObjectLookup;      // DataObjects name to Object pointer
  lookupGroup m_childGroups;  // child Groups: name->Group pointer
#if 0
  std::string m_name;
  DataShape m_dataShape;
#endif

public:

  /*!
   * @param parent name Pointer to DataGroup which contains this Group.
   * @param datastore Pointer to DataStore container.
   * \brief Constructor.
   */
  DataGroup(DataGroup *parent, DataStore *datastore) :
    m_parent(parent),
    m_datastore(datastore)
  {};

  /*!
   * @param name Name to check.
   * \brief Return true if the name exists in this DataGroup.
   */
  bool HasName(const std::string& name);

  /*!
   * @param name Name for created DataObject.
   * \brief Create a DataObject and add to this DataGroup.
   */
  DataObject *CreateDataObject(const std::string& name) { return AddDataObject(name, NULL); }

  /*!
   * @param name Name of DataObject to add.
   * @param obj  Pointer to an existing DataObject.
   * \brief Add existing DataObject to this DataGroup.
   */
  DataObject *AddDataObject(const std::string& name, DataObject *obj);

  /*!
   * @param name Name of DataObject to find.
   * \brief Return pointer to DataObject.
   */
  DataObject *GetDataObject(const std::string& name)
  {
    const IDType indx = m_DataObjectLookup.at(name);
    return m_DataObjects[indx];
  }

  /*!
   * @param indx Index of DataObject within this DataGroup.
   * \brief Return pointer to DataObject.
   */
  DataObject *GetDataObject(const IDType indx)
  {
    DataObject *obj = m_DataObjects[indx];
    if (obj == NULL) {
      // Object has been deleted and index is a hole in the table.
      throw std::exception();
    }
    return obj;
  }

  /*!
   * @param name Name of DataObject to find.
   * \brief Return index of DataObject in this DataGroup.
   */
  IDType IndexDataObject(const std::string& name)
  {
    return m_DataObjectLookup.at(name);
  }

  /*!
   * @param obj Name of DataObject to find.
   * \brief Return name of DataObject in this DataGroup.
   */
  std::string const & NameDataObject(DataObject *obj);

  /*!
   * @param name Name of DataObject to remove.
   * \brief Remove named DataObject from the index.
   *   The DataObject still exists in the DataStore.
   */
  DataObject *RemoveDataObject(const std::string& name);

  /*!
   * @param obj Pointer to DataObject to remove.
   * \brief Remove obj from the index.
   *   The DataObject still exists in the DataStore.
   */
  DataObject *RemoveDataObject(DataObject *obj);

  /*!
   * \brief Remove all DataObjects from this DataGroup.
   */
  void ClearDataObjects();

  /*!
   * @param name Name of DataGroup to create.
   * \brief Create a new DataGroup within this DataGroup.
   */
  DataGroup *CreateDataGroup(const std::string& name);

  /*!
   * @param name Name of DataGroup to find.
   * \brief Return pointer to DataGroup.
   */
  DataGroup *GetDataGroup(const std::string& name)
  {
    return m_childGroups.at(name);
  }

  /*!
   * \brief Return number of DataObjects contained in this DataGroup.
   */
  size_t CountObjects() { return m_DataObjects.size(); }

  /*!
   * \brief Return number of DataGroups contained in this DataGroup.
   */
  size_t CountGroups() { return m_childGroups.size(); }

  /*!
   * \brief Return DataObjects contained in this DataGroup.
   */
  lookupType GetDataObjects() { return m_DataObjectLookup; }

  /*!
   * \brief Return DataGroups contained in this DataGroup.
   */
  lookupGroup GetDataGroups() { return m_childGroups; }
  

};

// DataStore.hpp ----------------------------------------

/**
 * \class DataStore
 *
 * \brief Class to own DataObject and a root DataGroup.
 */
class DataStore
{
public:
  /*!
   * \brief vector of DataObject pointers.
   */
  typedef std::vector< DataObject* > dataObjectContainerType;

  // constructor
  // creates empty root data group and names it "/"
  // calls reserve on data object vector to minimize later resizes
private:

  // use a counter to generate a new unique id
  // may also be useful to know the number of active data objects
  IDType m_IDCounter;

  // Root data group, created automatically with datastore.
  DataGroup m_RootGroup;

  // container of data object pointers
  // as long as we recycle ids, this should not have many vacancies with NULL pointers
  // if it's an issue, change this to a std::map (or boost/std::unordered_map)
  dataObjectContainerType m_DataObjects;

  // stack of unique ids that can be recycled
  std::stack< DataObject* > m_AvailableDataObjects;

public:
  /*!
   * \brief Constructor.
   */
  DataStore() :
    m_IDCounter(0),
    m_RootGroup(NULL, this) {};

  /*!
   * \brief Destructor.
   */
  ~DataStore();

  // copy constructor
  //DataStore( const DataStore* store );

  /*!
   * \brief Create an empty DataObject.
   *    It is assigned a universal id and owned by the DataStore
   */
  DataObject *CreateDataObject();

  /*!
   * @param obj DataObject to delete.
   * \brief Remove a DataObject from the DataStore.
   *   It is disassociated with all groups and returned to the free pool.
   */
  void DeleteDataObject( DataObject *obj);

  /*!
   * @param id  Universal id of the DataObject.
   * \brief Remove a DataObject from the DataStore.
   *   It is disassociated with all groups and returned to the free pool.
   */
  void DeleteDataObject(IDType id);

  /*!
   * \brief Return pointer to the root DataGroup.
   */
  DataGroup* GetRootDataGroup() { return &m_RootGroup; };

};

// --------------------------------------------------------------------------------
// DataObject.cpp ----------------------------------------


// DataGroup.cpp ----------------------------------------

bool DataGroup::HasName(const std::string& name) {
  // XXX must perform two lookups: objects and groups
  DataGroup::lookupType::iterator it = m_DataObjectLookup.find(name);
  if (it != m_DataObjectLookup.end())
    return true;
  DataGroup::lookupGroup::iterator itg = m_childGroups.find(name);
  if (itg != m_childGroups.end())
    return true;
  return false;
}

DataObject *DataGroup::AddDataObject(const std::string& name, DataObject *obj)
{
  if (HasName(name)) {
      throw std::exception();
  }
  if (obj == NULL) {
    obj = m_datastore->CreateDataObject();
  }
  m_DataObjectLookup[name] = m_DataObjects.size();  // map name to index
  m_DataObjects.push_back(obj);
  obj->AttachGroup(this);
  // XXX how does user get index? Search IndexDataObject(name) or return from here?
  return obj;
}

DataObject *DataGroup::RemoveDataObject(const std::string& name)
{
  DataGroup::lookupType::iterator it = m_DataObjectLookup.find(name);
  if (it != m_DataObjectLookup.end()) {
    IDType indx = it->second;
    DataObject *obj = m_DataObjects[indx];
    obj->DetachGroup(this);
    m_DataObjectLookup.erase(it);
    m_DataObjects[indx] = NULL;    // XXX remove from m_DataObjects
    // XXX this makes a hole in the table, but perserved existing indexes.
  } else {
    throw std::exception();
  }
}

DataObject *DataGroup::RemoveDataObject(DataObject *obj)
{
  std::string const & name = NameDataObject(obj);
  RemoveDataObject(name);
}

DataGroup *DataGroup::CreateDataGroup(const std::string& name)
{
  if (HasName(name)) {
      throw std::exception();
  }
  DataGroup *grp = new DataGroup(this, this->m_datastore);
  m_childGroups[name] = grp;
  return grp;
}

void DataGroup::ClearDataObjects()
{
  m_DataObjects.clear();
  m_DataObjectLookup.clear();
}

std::string const & DataGroup::NameDataObject(DataObject *obj)
{
  DataGroup::lookupType::iterator it;

  // XXX brute force search for obj
  for (it=m_DataObjectLookup.begin(); it != m_DataObjectLookup.end(); ++it) {
    IDType indx = it->second;
    DataObject *obj1 = m_DataObjects[indx];
    if (obj1 == obj) {
      return it->first;
    }
  }
  throw std::exception();
}


// DataStore.cpp ----------------------------------------

DataStore::~DataStore()
{
  for (size_t i=0; i < m_DataObjects.size(); i++) {
    delete m_DataObjects[i];
  }
}

DataObject *DataStore::CreateDataObject()
{
  // TODO: implement pool, look for free nodes.  Allocate in blocks.
  DataObject *obj;
  if (m_AvailableDataObjects.empty()) {
    obj = new DataObject(m_IDCounter);
    ++m_IDCounter;
    m_DataObjects.push_back(obj);
  } else {
    obj = m_AvailableDataObjects.top();
    m_AvailableDataObjects.pop();
  }
  
  return obj;
}

void DataStore::DeleteDataObject( DataObject *obj )
{
  DataObject::GroupSet *gset = obj->GetDataGroups();
  for (DataObject::GroupSet::iterator it = gset->begin(); it != gset->end(); ++it) {
    DataGroup *grp = *it;
    grp->RemoveDataObject(obj);
  }

  m_AvailableDataObjects.push(obj);
  return;
}

// --------------------------------------------------------------------------------
// test.cpp  ----------------------------------------

#include <cstdio>
#include <iostream>
#include <assert.h>
#include <string.h>

char blanks[] = "                 ";
// print nodes which have attribute
void print_tree(const std::string& name, DataGroup *grp, int indent, DataGroup *attr)
{
  int nobjs = grp->CountObjects();
  int ngrps = grp->CountGroups();

  printf("%.*sgroup %s   nobjs=%d ngrps=%d\n", indent, blanks, name.c_str(), nobjs, ngrps);

  indent += 2;

  {
    DataGroup::lookupType objs = grp->GetDataObjects();
    DataGroup::lookupType::iterator it;

    for (it=objs.begin(); it != objs.end(); ++it) {
      IDType indx = it->second;
      DataObject *obj = grp->GetDataObject(indx);
      const char *name = it->first.c_str();
      if (attr == NULL) {
	printf("%.*sobj %s  gid=%d uid=%ld\n", indent, blanks, name, indx, (long) obj->GetUID());
      } else if (obj->IsAttachedGroup(attr)) {
	printf("%.*sobj %s  uid=%ld\n", indent, blanks, name, (long) obj->GetUID());
      }
    }
  }

  {
    DataGroup::lookupGroup objs = grp->GetDataGroups();
    DataGroup::lookupGroup::iterator it;

    for (it=objs.begin(); it != objs.end(); ++it) {
      print_tree(it->first.c_str(), it->second, indent, attr);
    }
  }

  return;
}

void dump_search(DataGroup *grp, DataGroup *attrgrp)
{
  {
    DataGroup::lookupType objs = grp->GetDataObjects();
    DataGroup::lookupType::iterator it;

    for (it=objs.begin(); it != objs.end(); ++it) {
      IDType indx = it->second;
      DataObject *obj = grp->GetDataObject(indx);
      const char *name = it->first.c_str();
      if (obj->m_dump) {
	attrgrp->AddDataObject(it->first, obj);
      }
    }
  }

  {
    DataGroup::lookupGroup objs = grp->GetDataGroups();
    DataGroup::lookupGroup::iterator it;

    for (it=objs.begin(); it != objs.end(); ++it) {
      dump_search(it->second, attrgrp);
    }
  }

  return;
}


// Test routine with a simple tree
void test1()
{
  std::cout << "Test1" << std::endl;

  DataStore ds;

  DataGroup *root = ds.GetRootDataGroup();

  // XXX Id would be possible to skip calling IndexDataObject if CreatDataObject set id1
  DataObject *obj1 = root->CreateDataObject("one");
  IDType id1 = root->IndexDataObject("one");
  assert(id1 == 0);
  assert(obj1 == root->GetDataObject(id1));
  assert(obj1->GetUID() == 0);

  DataObject *obj2 = root->CreateDataObject("two");
  IDType id2 = root->IndexDataObject("two");
  assert(id2 == 1);
  assert(obj2 == root->GetDataObject(id2));
  assert(obj2->GetUID() == 1);

  DataGroup *grp3 = root->CreateDataGroup("three");
  //  DataGroup *grp3 = root->CreateDataGroup("three", &index); // return index relative to group
  DataObject *obj4 = grp3->CreateDataObject("four");

  assert(root->HasName("one") == true);    // object
  assert(root->HasName("three") == true);  // group
  assert(root->HasName("nosuchname") == false);

  assert(strcmp("one", root->NameDataObject(obj1).c_str()) == 0);

#if 0
  // raise exception because of duplicate names
  obj1 = root->CreateDataObject("one");
  grp3 = root->CreateDataGroup("three");
#endif

  print_tree("root", root, 0, NULL);
}

// Test moving a DataObject around in the tree.
void test2()
{
  std::cout << "Test2" << std::endl;

  // Create a tree
  DataStore ds;
  DataObject *obj = ds.CreateDataObject();  // unattached node
  DataGroup *root = ds.GetRootDataGroup();
  DataGroup *grp1 = root->CreateDataGroup("a");
  DataGroup *grp2 = root->CreateDataGroup("b");

  // Add object to two different groups
  grp1->AddDataObject("a", obj);
  grp2->AddDataObject("b", obj);
  //-----

  print_tree("root", root, 0, NULL);

  grp1->RemoveDataObject("a");
  print_tree("root", root, 0, NULL);
  grp2->RemoveDataObject("b");
  print_tree("root", root, 0, NULL);

  grp2->ClearDataObjects();

#if 0
  grp2->RemoveDataObject("nosuchname");
  // raise exception
#endif

  ds.DeleteDataObject(obj);

  // Should reuse DataObject just deleted.
  DataObject *obj2 = ds.CreateDataObject();
  assert(obj->GetUID() == 0);
}

// Test deleting object while still attached to groups.
void test3()
{
  std::cout << "Test3" << std::endl;

  // Create a tree
  DataStore ds;
  DataObject *obj = ds.CreateDataObject();  // unattached node
  DataGroup *root = ds.GetRootDataGroup();
  DataGroup *grp1 = root->CreateDataGroup("a");
  DataGroup *grp2 = root->CreateDataGroup("b");

  // Add to two different groups
  grp1->AddDataObject("a", obj);
  grp2->AddDataObject("b", obj);
  //-----

  print_tree("root-before", root, 0, NULL);

  // Delete while still associated with groups
  // This willl disassociate obj from groups
  ds.DeleteDataObject(obj);

  print_tree("root-after", root, 0, NULL);

}

void sample_tree(DataGroup *root)
{
  root->CreateDataObject("a");
  root->CreateDataObject("b");
  root->CreateDataObject("c");

  DataGroup *grpc = root->CreateDataGroup("d");
  grpc->CreateDataObject("e");
  grpc->CreateDataObject("f");

  DataGroup *grpg = root->CreateDataGroup("g");
  grpg->CreateDataObject("h");
  grpg->CreateDataObject("i");

  DataGroup *grpj = grpg->CreateDataGroup("j");
  grpj->CreateDataObject("k");
  grpj->CreateDataObject("l");

}

// Sample dump two ways.
// Attach object to a 'dump' group.
// Traverse tree looking for dump 'attribute'.
int test_dump()
{
  std::cout << "Test Dump" << std::endl;

  DataStore ds;
  DataGroup *root = ds.GetRootDataGroup();

  DataGroup *code = root->CreateDataGroup("code");
  DataGroup *attrs = root->CreateDataGroup("attrs");
  DataGroup *dump = attrs->CreateDataGroup("dump");
  DataGroup *dump2 = attrs->CreateDataGroup("dump2");

  sample_tree(code);
  print_tree("root", root, 0, NULL);

  DataGroup *grp;
  DataObject *obj;

  // Add some dump 'attributes'
  // When adding to to dump group, must create a name for each node since they may not be unique
  obj = code->GetDataObject("a");
  dump->AddDataObject("aa", obj);
  obj->m_dump = true;

  grp = code->GetDataGroup("g");
  obj = grp->GetDataObject("h");
  dump->AddDataObject("hh", obj);
  obj->m_dump = true;

  grp = grp->GetDataGroup("j");
  obj = grp->GetDataObject("k");
  dump->AddDataObject("kk", obj);
  obj->m_dump = true;

  // fill in group dump2 by doing a 'query'
  // Duplicate names are a problem because groups require a name.
  // A 'nameless' group would be useful.  No get by name, on get by index and iterate.
  dump_search(code, dump2);
  print_tree("root", root, 0, NULL);


  // Only print (or dump) the objects which associate the dump group (attribute)
  print_tree("code w/dump", code, 0, dump);

}

int main(int argc, char *argv[])
{


  test1();
  test2();
  test3();
  test_dump();


#if 0
  DataObject *x = root["one"];       // root.GetDataObject("one")
  DataGroup  *y = root["three"];
  //  DataObject *z = root["three/four"];  // want to have a version which will not split name
    // vs root.FindPath("three/four")  to avoid searching for / with default case
    //    ds.FindPath("/three/four")   only use paths with datastore - assume from /
#endif

  return 0;
}
