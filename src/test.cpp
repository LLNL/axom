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
  DataStore::CreateDataStore("myDS1");

  // get reference to a data store object.
  DataStore::DataGroup* const myDS1 = DataStore::GetDataStore("myDS1");


  //********************************************************************************************************************
  // ***** CREATE GROUP *****
  //********************************************************************************************************************


  // use top level datastore group to create group
  myDS1->CreateDataGroup("group1");

  // get reference to the group
  DataStore::DataGroup* const group1 = myDS1->GetDataGroup("group1");

  // use top level datastore group to create group, and assign to local reference
  DataStore::DataGroup* const group2 = myDS1->CreateDataGroup("group2");

  // use group to create group, and assign to local reference
  DataStore::DataGroup* const group1a = group1->CreateDataGroup("group1a");


  //********************************************************************************************************************
  // ***** CREATE OBJECTS *****
  //********************************************************************************************************************

  // create blank DataObject
  group1->CreateDataObject("data1");

  // get reference to blank DataObject
  DataStore::DataObject* const dataObj1 = group1->GetDataObject("data1");

  // specify size of dataObject
  dataObj1->SetDataShape( DataStore::DataShape() );

  // allocate data
  dataObj1->Allocate();

  // do all that stuff in one line
  DataStore::DataObject* const dataObj2 = group1->CreateDataObject("data2")->SetDataShape(DataStore::DataShape())->Allocate();

  // again
  DataStore::DataObject* const dataObj3 = group1a->CreateDataObject("data3")->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj4 = group1a->CreateDataObject("data4")->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj5 = group2->CreateDataObject("data5")->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj6 = group2->CreateDataObject("data6")->SetDataShape(DataStore::DataShape())->Allocate();


  // create regular c++ object and insert into datastore
  AnyOldClass myClass;
  group1->CreateDataObject("AnyOldClass")->SetDataPointer(&myClass);

  //********************************************************************************************************************
  // ***** ACCESS DATA *****
  //********************************************************************************************************************

  // use DataObject to get data
  auto * const data1 = dataObj1->GetData<DataStore::int32*>();
  auto * const data2 = dataObj2->GetData<DataStore::real64*>();

  // assuming we didn't want to have a DataObject* const, and didn't want to have to get a DataObject* const at all.
  auto * const data3 = group1a->GetDataObject("data3")->GetData<DataStore::int32*>();
  // - or, really what you might want is:
  auto * const data4 = group1a->GetData<DataStore::real64*>("data4");

  // a pointer is great, but what about the shape?
  DataStore::DataShape dataDesc5 = group2->GetDataObject("data5")->GetDataShape();
  auto * const data5 = DataStore::rtTypes::CastPtr<DataStore::int32*>(dataDesc5.m_dataPtr);
  const int ndims5 = dataDesc5.m_numDimensions;
  std::size_t* const dims5 = dataDesc5.m_dimensions;

  // -again, you might want to cut out the middle man
  DataStore::DataShape dataDesc6 = group2->GetDataShape("data6");
  auto * const data6 = DataStore::rtTypes::CastPtr<DataStore::int32*>(dataDesc6.m_dataPtr);
  const int ndims6 = dataDesc6.m_numDimensions;
  std::size_t* const dims6 = dataDesc6.m_dimensions;

  // get the AnyOldClass object out of storage
  AnyOldClass* const p_myClass = group1->GetData<AnyOldClass*>("AnyOldClass");
  // DO SOMETHING WITH THE DATA

  return 0;
}




// All DataObjects are owned at the top level, and attributes are used to create groupings.
int approach2()
{
  //********************************************************************************************************************
  //***** CREATE DATASTORE *****
  //********************************************************************************************************************

  // Create DataStore
  DataStore::CreateDataStore("myDS1");

  // get reference to a data store object.
  DataStore::DataGroup* const myDS1 = DataStore::GetDataStore("myDS1");


  //********************************************************************************************************************
  // ***** DO NOT CREATE GROUPS *****
  //********************************************************************************************************************


  //********************************************************************************************************************
  // ***** CREATE OBJECTS *****
  //********************************************************************************************************************

  // create blank DataObject
  myDS1->CreateDataObject("data1");

  // get reference to blank DataObject
  DataStore::DataObject* const dataObj1 = myDS1->GetDataObject("data1");

  // specify attribute "group1"
  dataObj1->SetAttribute( DataStore::Attribute("group1") );

  // specify size of dataObject
  dataObj1->SetDataShape( DataStore::DataShape() );

  // allocate data
  dataObj1->Allocate();

  // do all that stuff in one line
  DataStore::DataObject* const dataObj2 = myDS1->CreateDataObject("data2")->SetAttribute( DataStore::Attribute("group1") )->SetDataShape(DataStore::DataShape())->Allocate();

  // again
  DataStore::DataObject* const dataObj3 = myDS1->CreateDataObject("data3")->SetAttribute( DataStore::Attribute("group1a") )->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj4 = myDS1->CreateDataObject("data4")->SetAttribute( DataStore::Attribute("group1a") )->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj5 = myDS1->CreateDataObject("data5")->SetAttribute( DataStore::Attribute("group2") )->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj6 = myDS1->CreateDataObject("data6")->SetAttribute( DataStore::Attribute("group2") )->SetDataShape(DataStore::DataShape())->Allocate();

  // create regular c++ object and insert into datastore
  AnyOldClass myClass;
  myDS1->CreateDataObject("AnyOldClass")->SetDataPointer<AnyOldClass*>(&myClass);

  //********************************************************************************************************************
  // ***** ACCESS DATA *****
  //********************************************************************************************************************

  // Get groups based on attributes
  DataStore::DataGroup* const group1 = myDS1->GetDataGroup( DataStore::Attribute("group1") );
  DataStore::DataGroup* const group1a = myDS1->GetDataGroup( DataStore::Attribute("group1a") );
  DataStore::DataGroup* const group2 = myDS1->GetDataGroup( DataStore::Attribute("group2") );


  // use DataObject to get data
  auto * const data1 = dataObj1->GetData<DataStore::int32*>();
  auto * const data2 = dataObj2->GetData<DataStore::real64*>();

  // assuming we didn't want to have a DataObject* const, and didn't want to have to get a DataObject* const at all.
  auto * const data3 = group1a->GetDataObject("data3")->GetData<DataStore::int32*>();
  // - or, really what you might want is:
  auto * const data4 = group1a->GetData<DataStore::real64*>("data4");

  // a pointer is great, but what about the shape?
  DataStore::DataShape dataDesc5 = group2->GetDataObject("data5")->GetDataShape();
  auto * const data5 = DataStore::rtTypes::CastPtr<DataStore::int32*>(dataDesc5.m_dataPtr);
  const int ndims5 = dataDesc5.m_numDimensions;
  std::size_t* const dims5 = dataDesc5.m_dimensions;

  // -again, you might want to cut out the middle man
  DataStore::DataShape dataDesc6 = group2->GetDataShape("data6");
  auto * const data6 = DataStore::rtTypes::CastPtr<DataStore::int32*>(dataDesc5.m_dataPtr);
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
  DataStore::DataGroup* myDS1 = DataStore::CreateDataStore("myDS1");


  //********************************************************************************************************************
  // ***** DO NOT CREATE GROUPS *****
  //********************************************************************************************************************


  //********************************************************************************************************************
  // ***** CREATE OBJECTS *****
  //********************************************************************************************************************

  // create blank DataObject
  DataStore::CreateDataObject(myDS1, "data1");

  // get reference to blank DataObject
  DataStore::DataObject* const dataObj1 = DataStore::GetDataObject(myDS1,"data1");

  // specify attribute "group1"
  DataStore::SetAttribute( dataObj1, DataStore::Attribute("group1") );

  // specify size of dataObject
  DataStore::SetDataShape( dataObj1, DataStore::DataShape() );

  // allocate data
  DataStore::Allocate(dataObj1);

  // do all that stuff in one line
  DataStore::DataObject* const dataObj2 = DataStore::CreateDataObject(myDS1,"data2");
  DataStore::SetAttribute( dataObj2, DataStore::Attribute("group1") );
  DataStore::SetDataShape( dataObj2, DataStore::DataShape());
  DataStore::Allocate(dataObj2);

  // combo detached and member
  DataStore::DataObject* const dataObj3 = DataStore::CreateDataObject(myDS1,"data3")->SetAttribute( DataStore::Attribute("group1a") )->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj4 = DataStore::CreateDataObject(myDS1,"data4")->SetAttribute( DataStore::Attribute("group1a") )->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj5 = DataStore::CreateDataObject(myDS1,"data5")->SetAttribute( DataStore::Attribute("group2") )->SetDataShape(DataStore::DataShape())->Allocate();
  DataStore::DataObject* const dataObj6 = DataStore::CreateDataObject(myDS1,"data6")->SetAttribute( DataStore::Attribute("group2") )->SetDataShape(DataStore::DataShape())->Allocate();

  //********************************************************************************************************************
  // ***** ACCESS DATA *****
  //********************************************************************************************************************

  // Get groups based on attributes
  DataStore::DataGroup* const group1 =  DataStore::GetDataGroup( myDS1, DataStore::Attribute("group1") );
  DataStore::DataGroup* const group1a = DataStore::GetDataGroup( myDS1, DataStore::Attribute("group1a") );
  DataStore::DataGroup* const group2 =  DataStore::GetDataGroup( myDS1, DataStore::Attribute("group2") );


  // use DataObject to get data
  auto * const data1 = DataStore::GetData<DataStore::int32*>(dataObj1);
  auto * const data2 = DataStore::GetData<DataStore::real64*>(dataObj2);

  // assuming we didn't want to have a DataObject* const, and didn't want to have to get a DataObject* const at all.
  auto * const data3 = DataStore::GetData<DataStore::int32*>(dataObj3);
  // - or, really what you might want is:
  auto * const data4 = DataStore::GetData<DataStore::real64*>(dataObj4);

  // a pointer is great, but what about the shape?
  DataStore::DataShape dataDesc5 = *(DataStore::GetDataShape(group2));
  auto * const data5 = DataStore::rtTypes::CastPtr<DataStore::int32*>(dataDesc5.m_dataPtr);
  const int ndims5 = dataDesc5.m_numDimensions;
  std::size_t* const dims5 = dataDesc5.m_dimensions;

  // -again, you might want to cut out the middle man
  DataStore::DataShape dataDesc6 = group2->GetDataShape("data6");
  auto * const data6 = DataStore::rtTypes::CastPtr<DataStore::int32*>(dataDesc5.m_dataPtr);
  const int ndims6 = dataDesc6.m_numDimensions;
  std::size_t* const dims6 = dataDesc6.m_dimensions;


  // DO SOMETHING WITH THE DATA

  return 0;
}


#if 0


int usecases1()
{
  using DataStore::DataGroup;
  using DataStore::DataObject;
  using DataStore::DataShape;
  using DataStore::real64;
  using DataStore::Attribute;
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
  DataStore::CreateDataStore("/myDS1");

  // get reference to a data store object.
  DataStore::DataGroup* const myDS1 = DataStore::GetDataStore("/myDS1");



  //********************************************************************************************************************
  //***** DATA OBJECT CREATION/DECLARATION *****
  //********************************************************************************************************************

  // document the use cases and behavior for error cases.


  // using DataStore interface
  DataStore::CreateDataObject(myDS1,"dataObject1");


  // using DataStore Object
  myDS1.CreateDataObject("dataObject1").SetDataShape( DataStore::DataShape() ).Allocate();
  // maybe creation without indication of what is stored is not useful. Discuss later.
  myDS1.CreateDataObject("dataObject1");


  // add object with it
//  myDS1.CreateDataObject("dataObject1");

  // and apply data description as well
  myDS1.CreateDataObject("dataObject1", DataStore::DataShape() );




  //********************************************************************************************************************
  // USE CASE TO CREATE AND WORK WITH GROUP
  //********************************************************************************************************************

  // use detached interface function to create group
  DataStore::CreateDataGroup("/myDS1","myGroup");
//  DataStore::CreateDataGroup(myDS1,"myGroup");

  // use top level datastore group to create group
  myDS1.CreateDataGroup("myGroup");

  // get access to the group
  DataStore::DataGroup* const myGroup = myDS1.GetDataGroup("myGroup");

  // don't to this sort of thing
  //DataStore::DataGroup* const myGroup1 = myDS1.GetDataObject<DataGroup>("myGroup");




  //********************************************************************************************************************
  //***** GET DATA OBJECT *****
  //********************************************************************************************************************

  // through data store object
  DataStore::DataObject* const dataObj1 = myDS1.GetDataObject("dataObject1");

  DataStore::DataObject* const dataObj1 = myDS1.GetDataGroup("myGroup").GetDataObject("dataObject1");

  // through interface
  DataStore::DataObject* const dataObj1a = DataStore::DataQuery1("/myDS1/dataObject1");
  DataStore::DataObject* const dataObj1a = DataStore::DataQuery1(myDS1, "/Group/dataObject1");



  //********************************************************************************************************************
  // USE CASE ATTACH attribute
  //********************************************************************************************************************
  // tablebase.cxx in vista has attribute bool expression parsing. look at keybase.cxx.

  // set attribute

  DataStore::AttachAttribute("myDS1", Attribute());

  DataStore::AttachAttribute("myGroup", Attribute());

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
  dataObj1.SetDataShape(DataStore::DataShape());

  // apply data descriptor to group, and propagate to data objects...and sub groups???
  myGroup.SetDataShape(DataStore::DataShape());



  // ***** allocate data on DS managed object *****
  dataObj1.Allocate();
  dataObj1.Allocate( DataStore::DataShape() );

  // should we allow for chaining commands...thus you have to return a reference to the object for every set function.
  myDS1.CreateDataObject("dataObject2").SetDataShape(DataStore::DataShape()).Allocate();




  // if the caller wants to own the data
  std::vector<DataStore::real64> myData1(50);
  DataStore::real64* myData2 = new DataStore::real64[50];


  // this should create the DataObject implicitly

  // option to transfer ownership??

  // We can insert a std::vector data container via the DataObject
  myDS1.insert( "mydata1", myData1 );

  // Or we can insert through the Interface
  DataStore::insert( "/myDS1/mydata1", myData1 );

  // insert some local managed data from raw pointer
  DataStore::insert( "/myDS1/mydata2", myData2, DataStore::DataShape() );







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
  DataStore::Node* const ref = myGroup.remove("mydata1a");
  // The datastore no longer owns the node, it now belongs to the caller.
  unique_ptr<DataStore::Node> ptr = myGroup.remove("mydata1a");

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
  DataStore::InitLibrary();

  // create directory directly below root
  DataStore::CreateDirectory("/myDS1");

  // get reference to a data store object.
  DataStore::DataDirectory* const myDS1 = DataStore::GetDirectory("/myDS1");

  DataStore::DataDirectory* const root = DataStore::GetDirectory("/");
  DataStore::DataDirectory* const myDS1 = root.GetDirectory("myDS1");

  

  //********************************************************************************************************************

  // Foreign Nodes
  // A foreign directory looks up names in its own data structures and returns
  // Datastore objects which wrap its own data.
  // This would apply to things like Conduit, Python, PDBlib, HDFlib

  // Looking up nodes ends up calling __getattr__ in Python

  PythonStore::Directory *dir =  new PythonStore::Directory;
  DataStore::Group* const py = root.CreateNode("/python").SetObject(dir);

  root['sys'];     // references Python's sys module (assuming it was imported)

  root['sys']['os']['path'];


  //********************************************************************************************************************
  // explicitly finalize the library
  DataStore::FinLibrary();


  myDS1.CreateNode("dataObject2").SetObject(DataStore::DataShape());
  myDS1.CreateNode("fcnObject1").SetObject(DataStore::FunctionDescriptor());
  myDS1.CreateNode("type:typeObject1").SetObject(DataStore::TypeDescriptor());


  // Allow users to extend functionality by creating their own scope
  myDS1.CreateDataObject("user:userObject1").SetObject(user)


#endif

#if ALT_DEF == 1

  // The attribute server deals with "index" objects.
  // Similar to Group object.
  // They contain references to other Nodes, can be iterated over
  // but cannot be accessed by name.

  // Return an index object
  DataStore::DataIndex* index = DataStore::AttrServer(root, "dump");

  // Append to existing index object
  DataStore::DataIndex *index = new DataStore::DataIndex;
  DataStore::AttrServer(dir1, "dump", index);
  DataStore::AttrServer(dir2, "dump", index);


  // call actor instead of create index
  DataStore::AttrServer("/dir1", "dump", actor);

  // call selector to decided which notes to use, return index
  DataStore::AttrServer(dir1, selector);

  // call selector to decided which notes to use, call actor for selected nodes
  DataStore::AttrServer(dir1, selector, actor);

  // query its location in the tree
  DataStore::Node* const root = node.GetRoot();
  DataStore::Node* const parent = node.GetParent();

  // Return "/dir1/dir2/node"
  string path = node.GetPath()

#endif


}
#endif

int main2()
{




	std::cout<<"test"<<std::endl;

}
