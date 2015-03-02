/**
 *
 * test.cpp
 *
 *  Created on: Nov 17, 2014
 */

#include<iostream>
#include "DatastoreInterface.hpp"
#include<vector>


class AnyOldClass
{

};
/**
 *   Group1
 *       Group1a
 *   Group2

  have some data.

  Show how to construct that by,
  	  1) standard method. Create Groups, add data members.
      2) add data members with attributes.
 *
 * @return
 */
int approach1()
{
  //********************************************************************************************************************
  //***** CREATE DATASTORE *****
  //********************************************************************************************************************

  // Create DataStore
  DataStoreNS::CreateDataStore("myDS1");

  // get reference to a data store object.
  DataStoreNS::DataStore* const myDS1 = DataStoreNS::GetDataStore("myDS1");


  //********************************************************************************************************************
  // ***** CREATE GROUP *****
  //********************************************************************************************************************


  // use top level datastore group to create group
  myDS1->CreateDataGroup("group1");

  // get reference to the group
  DataStoreNS::DataGroup* const group1 = myDS1->GetDataGroup("group1");

  // use top level datastore group to create group, and assign to local reference
  DataStoreNS::DataGroup* const group2 = myDS1->CreateDataGroup("group2");

  // use group to create group, and assign to local reference
  DataStoreNS::DataGroup* const group1a = group1->CreateDataGroup("group1a");


  //********************************************************************************************************************
  // ***** CREATE OBJECTS *****
  //********************************************************************************************************************

  // create blank DataObject
  group1->CreateDataObject("data1");

  // get reference to blank DataObject
  DataStoreNS::DataObject* const dataObj1 = group1->GetDataObject("data1");

  // specify size of dataObject
  dataObj1->SetDataShape( DataStoreNS::DataShape() );

  // allocate data
  dataObj1->Allocate();

  // do all that stuff in one line
  DataStoreNS::DataObject* const dataObj2 = group1->CreateDataObject("data2")->SetDataShape(DataStoreNS::DataShape())->Allocate();

  // again
  DataStoreNS::DataObject* const dataObj3 = group1a->CreateDataObject("data3")->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj4 = group1a->CreateDataObject("data4")->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj5 = group2->CreateDataObject("data5")->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj6 = group2->CreateDataObject("data6")->SetDataShape(DataStoreNS::DataShape())->Allocate();


  // create regular c++ object and insert into datastore
  AnyOldClass myClass;
  group1->CreateDataObject("AnyOldClass")->SetDataPointer(&myClass);

  //********************************************************************************************************************
  // ***** ACCESS DATA *****
  //********************************************************************************************************************

  // use DataObject to get data
  auto * const data1 = dataObj1->GetData<DataStoreNS::int32*>();
  auto * const data2 = dataObj2->GetData<DataStoreNS::real64*>();

  // assuming we didn't want to have a DataObject* const, and didn't want to have to get a DataObject* const at all.
  auto * const data3 = group1a->GetDataObject("data3")->GetData<DataStoreNS::int32*>();
  // - or, really what you might want is:
  auto * const data4 = group1a->GetData<DataStoreNS::real64*>("data4");

  // a pointer is great, but what about the shape?
  DataStoreNS::DataShape dataDesc5 = group2->GetDataObject("data5")->GetDataShape();
  auto * const data5 = DataStoreNS::rtTypes::CastPtr<DataStoreNS::int32*>(dataDesc5.m_dataPtr);
  const int ndims5 = dataDesc5.m_numDimensions;
  std::size_t* const dims5 = dataDesc5.m_dimensions;

  // -again, you might want to cut out the middle man
  DataStoreNS::DataShape dataDesc6 = group2->GetDataShape("data6");
  auto * const data6 = DataStoreNS::rtTypes::CastPtr<DataStoreNS::int32*>(dataDesc6.m_dataPtr);
  const int ndims6 = dataDesc6.m_numDimensions;
  std::size_t* const dims6 = dataDesc6.m_dimensions;

  // get the AnyOldClass object out of storage
  AnyOldClass* const p_myClass = group1->GetData<AnyOldClass*>("AnyOldClass");
  // DO SOMETHING WITH THE DATA

  return 0;
}



/*
// All DataObjects are owned at the top level, and attributes are used to create groupings.
int approach2()
{
  //********************************************************************************************************************
  //***** CREATE DATASTORE *****
  //********************************************************************************************************************

  // Create DataStore
  DataStoreNS::CreateDataStore("myDS1");

  // get reference to a data store object.
  DataStoreNS::DataStore* const myDS1 = DataStoreNS::GetDataStore("myDS1");


  //********************************************************************************************************************
  // ***** DO NOT CREATE GROUPS *****
  //********************************************************************************************************************


  //********************************************************************************************************************
  // ***** CREATE OBJECTS *****
  //********************************************************************************************************************

  // create blank DataObject
  myDS1->CreateDataObject("data1");

  // get reference to blank DataObject
  DataStoreNS::DataObject* const dataObj1 = myDS1->GetDataObject("data1");

  // specify attribute "group1"
  dataObj1->SetAttribute( DataStoreNS::Attribute("group1") );

  // specify size of dataObject
  dataObj1->SetDataShape( DataStoreNS::DataShape() );

  // allocate data
  dataObj1->Allocate();

  // do all that stuff in one line
  DataStoreNS::DataObject* const dataObj2 = myDS1->CreateDataObject("data2")->SetAttribute( DataStoreNS::Attribute("group1") )->SetDataShape(DataStoreNS::DataShape())->Allocate();

  // again
  DataStoreNS::DataObject* const dataObj3 = myDS1->CreateDataObject("data3")->SetAttribute( DataStoreNS::Attribute("group1a") )->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj4 = myDS1->CreateDataObject("data4")->SetAttribute( DataStoreNS::Attribute("group1a") )->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj5 = myDS1->CreateDataObject("data5")->SetAttribute( DataStoreNS::Attribute("group2") )->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj6 = myDS1->CreateDataObject("data6")->SetAttribute( DataStoreNS::Attribute("group2") )->SetDataShape(DataStoreNS::DataShape())->Allocate();

  // create regular c++ object and insert into datastore
  AnyOldClass myClass;
  myDS1->CreateDataObject("AnyOldClass")->SetDataPointer<AnyOldClass*>(&myClass);

  //********************************************************************************************************************
  // ***** ACCESS DATA *****
  //********************************************************************************************************************

  // Get groups based on attributes
  DataStoreNS::DataGroup* const group1 = myDS1->GetDataGroup( DataStoreNS::Attribute("group1") );
  DataStoreNS::DataGroup* const group1a = myDS1->GetDataGroup( DataStoreNS::Attribute("group1a") );
  DataStoreNS::DataGroup* const group2 = myDS1->GetDataGroup( DataStoreNS::Attribute("group2") );


  // use DataObject to get data
  auto * const data1 = dataObj1->GetData<DataStoreNS::int32*>();
  auto * const data2 = dataObj2->GetData<DataStoreNS::real64*>();

  // assuming we didn't want to have a DataObject* const, and didn't want to have to get a DataObject* const at all.
  auto * const data3 = group1a->GetDataObject("data3")->GetData<DataStoreNS::int32*>();
  // - or, really what you might want is:
  auto * const data4 = group1a->GetData<DataStoreNS::real64*>("data4");

  // a pointer is great, but what about the shape?
  DataStoreNS::DataShape dataDesc5 = group2->GetDataObject("data5")->GetDataShape();
  auto * const data5 = DataStoreNS::rtTypes::CastPtr<DataStoreNS::int32*>(dataDesc5.m_dataPtr);
  const int ndims5 = dataDesc5.m_numDimensions;
  std::size_t* const dims5 = dataDesc5.m_dimensions;

  // -again, you might want to cut out the middle man
  DataStoreNS::DataShape dataDesc6 = group2->GetDataShape("data6");
  auto * const data6 = DataStoreNS::rtTypes::CastPtr<DataStoreNS::int32*>(dataDesc5.m_dataPtr);
  const int ndims6 = dataDesc6.m_numDimensions;
  std::size_t* const dims6 = dataDesc6.m_dimensions;



  // DO SOMETHING WITH THE DATA

  return 0;
}

// same as approach2, but with detached interface.
int approach3()
{
  //********************************************************************************************************************
  //***** CREATE DATASTORE *****
  //********************************************************************************************************************

  // Create DataStore
  DataStoreNS::DataStore* myDS1 = DataStoreNS::CreateDataStore("myDS1");


  //********************************************************************************************************************
  // ***** DO NOT CREATE GROUPS *****
  //********************************************************************************************************************


  //********************************************************************************************************************
  // ***** CREATE OBJECTS *****
  //********************************************************************************************************************

  // create blank DataObject
  DataStoreNS::CreateDataObject(myDS1, "data1");

  // get reference to blank DataObject
  DataStoreNS::DataObject* const dataObj1 = DataStoreNS::GetDataObject(myDS1,"data1");

  // specify attribute "group1"
  DataStoreNS::SetAttribute( dataObj1, DataStoreNS::Attribute("group1") );

  // specify size of dataObject
  DataStoreNS::SetDataShape( dataObj1, DataStoreNS::DataShape() );

  // allocate data
  DataStoreNS::Allocate(dataObj1);

  // do all that stuff in one line
  DataStoreNS::DataObject* const dataObj2 = DataStoreNS::CreateDataObject(myDS1,"data2");
  DataStoreNS::SetAttribute( dataObj2, DataStoreNS::Attribute("group1") );
  DataStoreNS::SetDataShape( dataObj2, DataStoreNS::DataShape());
  DataStoreNS::Allocate(dataObj2);

  // combo detached and member
  DataStoreNS::DataObject* const dataObj3 = DataStoreNS::CreateDataObject(myDS1,"data3")->SetAttribute( DataStoreNS::Attribute("group1a") )->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj4 = DataStoreNS::CreateDataObject(myDS1,"data4")->SetAttribute( DataStoreNS::Attribute("group1a") )->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj5 = DataStoreNS::CreateDataObject(myDS1,"data5")->SetAttribute( DataStoreNS::Attribute("group2") )->SetDataShape(DataStoreNS::DataShape())->Allocate();
  DataStoreNS::DataObject* const dataObj6 = DataStoreNS::CreateDataObject(myDS1,"data6")->SetAttribute( DataStoreNS::Attribute("group2") )->SetDataShape(DataStoreNS::DataShape())->Allocate();

  //********************************************************************************************************************
  // ***** ACCESS DATA *****
  //********************************************************************************************************************

  // Get groups based on attributes
  DataStoreNS::DataGroup* const group1 =  DataStoreNS::GetDataGroup( myDS1, DataStoreNS::Attribute("group1") );
  DataStoreNS::DataGroup* const group1a = DataStoreNS::GetDataGroup( myDS1, DataStoreNS::Attribute("group1a") );
  DataStoreNS::DataGroup* const group2 =  DataStoreNS::GetDataGroup( myDS1, DataStoreNS::Attribute("group2") );


  // use DataObject to get data
  auto * const data1 = DataStoreNS::GetData<DataStoreNS::int32*>(dataObj1);
  auto * const data2 = DataStoreNS::GetData<DataStoreNS::real64*>(dataObj2);

  // assuming we didn't want to have a DataObject* const, and didn't want to have to get a DataObject* const at all.
  auto * const data3 = DataStoreNS::GetData<DataStoreNS::int32*>(dataObj3);
  // - or, really what you might want is:
  auto * const data4 = DataStoreNS::GetData<DataStoreNS::real64*>(dataObj4);

  // a pointer is great, but what about the shape?
  DataStoreNS::DataShape dataDesc5 = *(DataStoreNS::GetDataShape(group2));
  auto * const data5 = DataStoreNS::rtTypes::CastPtr<DataStoreNS::int32*>(dataDesc5.m_dataPtr);
  const int ndims5 = dataDesc5.m_numDimensions;
  std::size_t* const dims5 = dataDesc5.m_dimensions;

  // -again, you might want to cut out the middle man
  DataStoreNS::DataShape dataDesc6 = group2->GetDataShape("data6");
  auto * const data6 = DataStoreNS::rtTypes::CastPtr<DataStoreNS::int32*>(dataDesc5.m_dataPtr);
  const int ndims6 = dataDesc6.m_numDimensions;
  std::size_t* const dims6 = dataDesc6.m_dimensions;


  // DO SOMETHING WITH THE DATA

  return 0;
}
*/


#if 0


int usecases1()
{
  using DataStoreNS::DataGroup;
  using DataStoreNS::DataObject;
  using DataStoreNS::DataShape;
  using DataStoreNS::real64;
  using DataStoreNS::Attribute;
//  using namespace DataStore;

#define ALT_DEF 0
  /********************************************************************************************************************
   ****** Definitions of terms *****
   * - "namespace DataStore" contains all classes and interface function intrinsic to the Datastore.
   *
   * - "class DataObject" refers to a fundamental type that contains data, a description of the data, and information
   * about where the object resides in the overall structure of the datastore.
   *  XXX - name as Node instead of DataObject?
   *      contains connectivity and attribute information and reference to DataNode/ArrayNode
   *
   * - "class DataGroup : public DataObject" is a type that is used to interface with a collection of DataObjects. It
   * may or may not actually have the ability to own the collection of DataObjects that it provides access to.
   *  XXX - name as Directory instead of Group?
   *
   */
  //********************************************************************************************************************
  //***** DATASTORE CREATION and EXTRACTION *****
  //********************************************************************************************************************






  // use the interface to create a data store object
  DataStoreNS::CreateDataStore("/myDS1");

  // get reference to a data store object.
  DataStoreNS::DataGroup* const myDS1 = DataStoreNS::GetDataStore("/myDS1");



  //********************************************************************************************************************
  //***** DATA OBJECT CREATION/DECLARATION *****
  //********************************************************************************************************************

  // document the use cases and behavior for error cases.


  // using DataStore interface
  DataStoreNS::CreateDataObject(myDS1,"dataObject1");


  // using DataStore Object
  myDS1.CreateDataObject("dataObject1").SetDataShape( DataStoreNS::DataShape() ).Allocate();
  // maybe creation without indication of what is stored is not useful. Discuss later.
  myDS1.CreateDataObject("dataObject1");


  // add object with it
//  myDS1.CreateDataObject("dataObject1");

  // and apply data description as well
  myDS1.CreateDataObject("dataObject1", DataStoreNS::DataShape() );




  //********************************************************************************************************************
  // USE CASE TO CREATE AND WORK WITH GROUP
  //********************************************************************************************************************

  // use detached interface function to create group
  DataStoreNS::CreateDataGroup("/myDS1","myGroup");
//  DataStore::CreateDataGroup(myDS1,"myGroup");

  // use top level datastore group to create group
  myDS1.CreateDataGroup("myGroup");

  // get access to the group
  DataStoreNS::DataGroup* const myGroup = myDS1.GetDataGroup("myGroup");

  // don't to this sort of thing
  //DataStore::DataGroup* const myGroup1 = myDS1.GetDataObject<DataGroup>("myGroup");




  //********************************************************************************************************************
  //***** GET DATA OBJECT *****
  //********************************************************************************************************************

  // through data store object
  DataStoreNS::DataObject* const dataObj1 = myDS1.GetDataObject("dataObject1");

  DataStoreNS::DataObject* const dataObj1 = myDS1.GetDataGroup("myGroup").GetDataObject("dataObject1");

  // through interface
  DataStoreNS::DataObject* const dataObj1a = DataStoreNS::DataQuery1("/myDS1/dataObject1");
  DataStoreNS::DataObject* const dataObj1a = DataStoreNS::DataQuery1(myDS1, "/Group/dataObject1");



  //********************************************************************************************************************
  // USE CASE ATTACH attribute
  //********************************************************************************************************************
  // tablebase.cxx in vista has attribute bool expression parsing. look at keybase.cxx.

  // set attribute

  DataStoreNS::AttachAttribute("myDS1", Attribute());

  DataStoreNS::AttachAttribute("myGroup", Attribute());

  myDS1.SetAttribute("dump", true, 1);

  // XXX every node with a dump attribute needs to use the same type
  //     boolean in this case.
  //     Should attributes be defined first?
  root.DefineBoolAttribute("dump");

  // or the first call to set it may implicitly define it and any other type
  // would be an error.
  myDS1.SetAttribute("dump", 1);   // error


  msDS1.HasAttribute("dump");

  //********************************************************************************************************************
  //***** DATA DESCRIPTION/ALLOCATION/INSERTION *****
  //********************************************************************************************************************



  // apply data descriptor through object
  dataObj1.SetDataShape(DataStoreNS::DataShape());

  // apply data descriptor to group, and propagate to data objects...and sub groups???
  myGroup.SetDataShape(DataStoreNS::DataShape());



  // ***** allocate data on DS managed object *****
  dataObj1.Allocate();
  dataObj1.Allocate( DataStoreNS::DataShape() );

  // should we allow for chaining commands...thus you have to return a reference to the object for every set function.
  myDS1.CreateDataObject("dataObject2").SetDataShape(DataStoreNS::DataShape()).Allocate();




  // if the caller wants to own the data
  std::vector<DataStoreNS::real64> myData1(50);
  DataStoreNS::real64* myData2 = new DataStoreNS::real64[50];


  // this should create the DataObject implicitly

  // option to transfer ownership??

  // We can insert a std::vector data container via the DataObject
  myDS1.insert( "mydata1", myData1 );

  // Or we can insert through the Interface
  DataStoreNS::insert( "/myDS1/mydata1", myData1 );

  // insert some local managed data from raw pointer
  DataStoreNS::insert( "/myDS1/mydata2", myData2, DataStoreNS::DataShape() );







  //********************************************************************************************************************
  // USE CASE TO GET DATA
  //********************************************************************************************************************
  // get data objects
  DataShape* const dataObj1Desc = dataObj1.GetDataDescriptor();


  real64* const data1  = dataObj1.GetData<real64>();
  real64* const data1a = myDS1.GetData<real64>("mydata1");

//  DataDescriptor* const data1a = myDS1.GetDataDesc<real8>("mydata1");

  // scalars are treated as arrays with length1. must pass address.






  // how do we check to make sure that the data we may be trying to insert is
  // not duplicated in the DS. The keys are unique to the DS, so one key can only
  // point to one DataObject.
  myGroup.insert("mydata1a", mydata1a);
  // or we can just enforce that data cannot be inserted through the group construct

  // insert into a specific slot of the group to match some enum in user code
  myGroup.insert("mydata1a", mydata1a, index);

  myGroup.delete("mydata1a");

  // Remove from a group and keep a reference to the node.
  DataStoreNS::Node* const ref = myGroup.remove("mydata1a");
  // The datastore no longer owns the node, it now belongs to the caller.
  unique_ptr<DataStoreNS::Node> ptr = myGroup.remove("mydata1a");

  int i = node.GetIndex("mydata1a");

  real64* const data1 = node[i].getData<real64>("mydata1");


  // duplicates node, duplicating a group with recursively duplicate members
  parent.CopyNode("myGroup", "myGroup2");

  // renames myGroup as myGroup2
  parent.MoveNode("myGroup", "myGroup2");


  // inerate over a group
  for (auto myGroup) {
  }





}

int usecases2()
{

#if ALT_DEF == 1
  // explicitly initialize the library
  // creates the root directory, /, which is identical to any other Group
  DataStoreNS::InitLibrary();

  // create directory directly below root
  DataStoreNS::CreateDirectory("/myDS1");

  // get reference to a data store object.
  DataStoreNS::DataDirectory* const myDS1 = DataStoreNS::GetDirectory("/myDS1");

  DataStoreNS::DataDirectory* const root = DataStoreNS::GetDirectory("/");
  DataStoreNS::DataDirectory* const myDS1 = root.GetDirectory("myDS1");

  

  //********************************************************************************************************************

  // Foreign Nodes
  // A foreign directory looks up names in its own data structures and returns
  // Datastore objects which wrap its own data.
  // This would apply to things like Conduit, Python, PDBlib, HDFlib

  // Looking up nodes ends up calling __getattr__ in Python

  PythonStore::Directory *dir =  new PythonStore::Directory;
  DataStoreNS::Group* const py = root.CreateNode("/python").SetObject(dir);

  root['sys'];     // references Python's sys module (assuming it was imported)

  root['sys']['os']['path'];


  //********************************************************************************************************************
  // explicitly finalize the library
  DataStoreNS::FinLibrary();


  myDS1.CreateNode("dataObject2").SetObject(DataStoreNS::DataShape());
  myDS1.CreateNode("fcnObject1").SetObject(DataStoreNS::FunctionDescriptor());
  myDS1.CreateNode("type:typeObject1").SetObject(DataStoreNS::TypeDescriptor());


  // Allow users to extend functionality by creating their own scope
  myDS1.CreateDataObject("user:userObject1").SetObject(user)


#endif

#if ALT_DEF == 1

  // The attribute server deals with "index" objects.
  // Similar to Group object.
  // They contain references to other Nodes, can be iterated over
  // but cannot be accessed by name.

  // Return an index object
  DataStoreNS::DataIndex* index = DataStoreNS::AttrServer(root, "dump");

  // Append to existing index object
  DataStoreNS::DataIndex *index = new DataStoreNS::DataIndex;
  DataStoreNS::AttrServer(dir1, "dump", index);
  DataStoreNS::AttrServer(dir2, "dump", index);


  // call actor instead of create index
  DataStoreNS::AttrServer("/dir1", "dump", actor);

  // call selector to decided which notes to use, return index
  DataStoreNS::AttrServer(dir1, selector);

  // call selector to decided which notes to use, call actor for selected nodes
  DataStoreNS::AttrServer(dir1, selector, actor);

  // query its location in the tree
  DataStoreNS::Node* const root = node.GetRoot();
  DataStoreNS::Node* const parent = node.GetParent();

  // Return "/dir1/dir2/node"
  string path = node.GetPath()

#endif


}
#endif

int main2()
{




	std::cout<<"test"<<std::endl;

}
