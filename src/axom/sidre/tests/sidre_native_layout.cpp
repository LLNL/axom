// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre/core/sidre.hpp"

#include "conduit_blueprint.hpp"

#include <vector>
#include <string>
#include <iostream>


using axom::sidre::Buffer;
using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
using axom::sidre::indexIsValid;
using axom::sidre::InvalidName;
using axom::sidre::nameIsValid;

//------------------------------------------------------------------------------

namespace
{

const axom::sidre::TypeID DOUBLE_ID = axom::sidre::DOUBLE_ID;
const axom::sidre::TypeID INT32_ID = axom::sidre::INT32_ID;

/**
 * \brief Simple utility function to initialize an array
 * \param data A pointer to the beginning of the array
 * \param size The number of elements in the array
 * \param initVal The starting value of the element
 * \param intDiv An integer divisor for the data's index
 * \param scaleFac A scaling factor for the data's index
 *
 * The value of element idx will be initVal + ( idx / intDiv) * scaleFac
 */
template<typename T>
void setData(T* data, int size, T initVal=T(0), int intDiv=1, T scaleFac=T(1))
{
  for(int i=0 ; i< size ; ++i)
    data[i] = initVal + ((i+1) / intDiv) * scaleFac;
}

/**
 * \brief Utility function to compare the corresponding node and view
 * \param path The path from the datastore root to the view
 * \param rootNode The root node of the layout
 * \param rootGroup The root group of the datastore
 */
template<typename T>
void checkPointersAndData(const std::string& path
                          , axom::sidre::Node& rootNode
                          , axom::sidre::Group* rootGroup)
{
  axom::sidre::Node& node = rootNode[path];
  T* nD = static_cast<T*>(node.element_ptr(0));

  axom::sidre::View* view = rootGroup->getView(path);
  EXPECT_TRUE(nullptr != view);

  T* vD = view->getData<T*>();
  EXPECT_TRUE(nullptr != vD);

  EXPECT_EQ(nD, vD)
    << "Error when comparing pointer address between "
    << "conduit native layout and Sidre View "
    << "for path " << path
    <<".\n\t Conduit address: " << nD
    <<"\n\t Datastore view address: " << vD;

  EXPECT_EQ(*nD, *vD)
    << "Error when comparing values between "
    << "conduit native layout and Sidre View "
    << "for path " << path
    <<".\n\t Conduit value: " << *nD
    <<"\n\t View value: " << *vD;
}

/**
 * \brief Template specialization of checkPointersAndData for scalar views
 *        containing strings
 */
template<>
void checkPointersAndData<std::string>(const std::string& path
                                       , axom::sidre::Node& rootNode
                                       ,
                                       axom::sidre::Group* rootGroup)
{
  axom::sidre::Node& node = rootNode[path];
  std::string nD = node.as_string();

  axom::sidre::View* view = rootGroup->getView(path);
  EXPECT_TRUE(nullptr != view);
  EXPECT_TRUE(view->isString());

  std::string vD(view->getString());

  EXPECT_EQ(nD, vD)
    << "Error when comparing strings between "
    << "conduit native layout and Sidre View "
    << "for path " << path
    <<".\n\t Conduit string: " << nD
    <<"\n\t Datastore view string: " << vD;

}
}


TEST(sidre_native_layout,empty_layout)
{
  DataStore* ds   = new DataStore();

  SLIC_INFO("Testing sidre_layout function on empty datastore");

  SLIC_INFO("****** Visual diagnostics ******");
  SLIC_INFO("Printing datastore");
  ds->print();

  SLIC_INFO("Printing datastore native layout:");
  axom::sidre::Node node;
  ds->getRoot()->createNativeLayout(node);
  node.to_json_stream(std::cout);

  SLIC_INFO("****** done ******");


  if(ds != nullptr)
  {
    delete ds;
    ds = nullptr;
  }
}

TEST(sidre_native_layout,generate_native_layout)
{
  /// Allocate and initialize two external arrays
  // The REAL_BUF has 10 copies of each integer
  const int EXT_REAL_BUF_SIZE = 20;
  double* extRealPtr = new double[EXT_REAL_BUF_SIZE];
  setData(extRealPtr, EXT_REAL_BUF_SIZE, 0., 1, 1/10.);

  // The INT_BUF is a decreasing sequence
  const int EXT_INT_BUF_SIZE = 5;
  int* extIntPtr = new int[EXT_INT_BUF_SIZE];
  setData(extIntPtr, EXT_INT_BUF_SIZE, EXT_INT_BUF_SIZE, 1, -1);


  DataStore* ds   = new DataStore();
  Group* root = ds->getRoot();

  // Setup a buffer that will have two views
  const int REAL_BUF_SIZE = 100;
  Buffer* realBuf = ds->createBuffer(DOUBLE_ID, REAL_BUF_SIZE)->allocate();
  setData<double>(realBuf->getData(), REAL_BUF_SIZE, 0., 10, 1.);

  // create the views using the path syntax
  root->createView("Ga/Va")->attachBuffer(realBuf)->apply(50,0);
  root->createView("Ga/Vb")->attachBuffer(realBuf)->apply(20, 50);

  root->createView("Ga/Gc/Ve",DOUBLE_ID, 10)->allocate();
  root->createView("Ga/Gc/Vf")->setExternalDataPtr(DOUBLE_ID, EXT_REAL_BUF_SIZE,
                                                   extRealPtr);

  root->createView("Gb/Vc",INT32_ID, 10)->allocate();
  root->createView("Gb/Vd")->setExternalDataPtr(INT32_ID, EXT_INT_BUF_SIZE,
                                                extIntPtr);

  // Set some data in the attached buffers so we can distinguish in our tests
  setData<double>(root->getView("Ga/Gc/Ve")->getData(), 10, 0., 1, 1.1);
  setData<int>(root->getView("Gb/Vc")->getData(), 10, 0, 1, 2);

  SLIC_INFO("****** Visual diagnostics ******");
  SLIC_INFO("Printing datastore");
  ds->print();

  SLIC_INFO("Printing datastore native layout:");
  axom::sidre::Node node;
  ds->getRoot()->createNativeLayout(node);
  node.to_json_stream(std::cout);
  //node.save("nativeLayoutTest.conduit");

  SLIC_INFO("done.");

  SLIC_INFO("*** Checking that addresses and first elements match.");
  checkPointersAndData<double>("Ga/Va", node, root);
  checkPointersAndData<double>("Ga/Va", node, root);
  checkPointersAndData<double>("Ga/Gc/Ve", node, root);
  checkPointersAndData<double>("Ga/Gc/Vf", node, root);
  checkPointersAndData<int>("Gb/Vc", node, root);
  checkPointersAndData<int>("Gb/Vd", node, root);
  SLIC_INFO("done.");



  /// Clean up memory
  if(extRealPtr != nullptr)
  {
    delete[] extRealPtr;
    extRealPtr = nullptr;
  }
  if(extIntPtr != nullptr)
  {
    delete[] extIntPtr;
    extIntPtr = nullptr;
  }

  if(ds != nullptr)
  {
    delete ds;
    ds = nullptr;
  }
}


TEST(sidre_native_layout,native_layout_with_scalars)
{

  DataStore* ds   = new DataStore();
  Group* root = ds->getRoot();

  root->createView("Garray/Vdbl20",DOUBLE_ID, 20)->allocate();
  root->createView("Garray/Vint10",INT32_ID, 10)->allocate();

  // Set some data in the attached buffers so we can distinguish in our tests
  setData<double>(root->getView("Garray/Vdbl20")->getData(), 20, 0., 1, 1.1);
  setData<int>(root->getView("Garray/Vint10")->getData(), 10, 0, 1, 2);

  root->createView("Gscalar/Vint")->setScalar(1);
  root->createView("Gscalar/Vdbl")->setScalar(2.);
  root->createView("Gscalar/Vstr")->setString("Hello");


  SLIC_INFO("****** Visual diagnostics ******");
  SLIC_INFO("Printing datastore::info():");
  ds->print();

  SLIC_INFO("Printing datastore::layout():");
  axom::sidre::Node node;
  ds->getRoot()->createNativeLayout(node);
  node.to_json_stream(std::cout);
  //node.save("nativeLayoutTest.conduit");

  SLIC_INFO("done.");
  SLIC_CHECK_MSG(false,"Some message");


  SLIC_INFO("*** Checking that addresses and first elements match.");
  checkPointersAndData<double>("Garray/Vdbl20", node, root);
  checkPointersAndData<int>("Garray/Vint10", node, root);
  checkPointersAndData<int>("Gscalar/Vint", node, root);
  checkPointersAndData<double>("Gscalar/Vdbl", node, root);
  checkPointersAndData<std::string>("Gscalar/Vstr", node, root);

  if(ds != nullptr)
  {
    delete ds;
    ds = nullptr;
  }
}


//------------------------------------------------------------------------------
TEST(sidre_native_layout,export_import_conduit)
{
  // Checks that we can import an exported conduit tree into a separate
  // datastore using Group::createNativeLayout() and Group::importConduitTree()

  // Tests with non-zero- and with zero-sized arrays
  for( int SZ : {10, 0})
  {
    SLIC_INFO("Checking export/import from conduit with array of size " << SZ);

    axom::sidre::DataStore ds1, ds2;
    conduit::Node node1, node2;

    int* conduit_data = nullptr;

    std::string grpName = "some_group";
    std::string viewName = "array_view";

    // Simple lambda to return the expected value at an index
    auto expValue = [=](int i) {
                      return SZ-i;
                    };

    // Create datastore with an array view and export to conduit node
    {
      // Create views and groups
      auto* root = ds1.getRoot();
      auto* grp  = root->createGroup(grpName);
      auto* view = grp->createViewAndAllocate(viewName, INT32_ID, SZ);

      // initialize the data
      int* sidre_data = view->getData();
      EXPECT_NE(nullptr, sidre_data);
      for(int i=0 ; i<SZ ; ++i)
      {
        sidre_data[i] = expValue(i);
      }

      // Dump to console
      std::cout << "Datastore1 after initialization: " << std::endl;
      root->print();
      std::cout << std::endl;

      // Convert datastore to conduit layout and compare data
      bool success = root->createNativeLayout(node1);
      EXPECT_TRUE(success);

      conduit_data = node1[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, conduit_data);
      for(int i=0 ; i<SZ ; ++i)
      {
        EXPECT_EQ(conduit_data[i], sidre_data[i]);
      }

      // Dump to console
      std::cout << "Node1 after createNativeLayout: " << std::endl;
      node1.print();
      std::cout << std::endl;
    }

    // import conduit node into new datastore and check validity
    {
      auto* root = ds2.getRoot();
      bool success = root->importConduitTree(node1);
      EXPECT_TRUE(success);

      // Dump to console
      std::cout << "Datastore2 after importing conduit tree: " << std::endl;
      root->print();
      std::cout << std::endl;

      // Check for expected views and groups
      EXPECT_TRUE(root->hasGroup(grpName));
      auto* grp = root->getGroup(grpName);
      EXPECT_NE(nullptr, grp);

      EXPECT_TRUE(grp->hasView(viewName));
      auto* view = grp->getView(viewName);
      EXPECT_NE(nullptr, view);

      EXPECT_TRUE(view->isAllocated());
      EXPECT_EQ(SZ, view->getNumElements());

      // compare data
      int* sidre_data = view->getData();
      EXPECT_NE(nullptr, sidre_data);
      for(int i=0 ; i<SZ ; ++i)
      {
        EXPECT_EQ(conduit_data[i], sidre_data[i]);
      }

      // confirm that this is a distinct copy of the first data store's data
      auto* otherView = ds1.getRoot()->getGroup(grpName)->getView(viewName);
      EXPECT_NE(sidre_data, otherView->getData<int*>() );

      // Convert to conduit layout and compare data
      success = root->createNativeLayout(node2);
      EXPECT_TRUE(success);

      conduit_data = node2[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, conduit_data);
      for(int i=0 ; i<SZ ; ++i)
      {
        EXPECT_EQ(conduit_data[i], sidre_data[i]);
      }

      // Dump to console
      std::cout << "Node2 after createNativeLayout: " << std::endl;
      node2.print();
      std::cout << std::endl;
    }

    // compare the two conduit nodes
    // first, check compatiblity
    bool compat = node1.compatible(node2);
    EXPECT_TRUE(compat);

    // Finally, check the actual array data
    {
      int* arr1 = node1[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, arr1);
      int* arr2 = node2[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, arr2);

      for(int i=0 ; i<SZ ; ++i)
      {
        EXPECT_EQ(arr1[i], arr2[i]);
        EXPECT_EQ(expValue(i), arr2[i]);
      }
    }
  }
}

TEST(sidre_native_layout,import_conduit_and_verify_protocol)
{
  // This tests is similar to the 'export_import_conduit' test above,
  // but uses conduit's blueprint_verify on the result
  // which is expected to be a blueprint multicomponent array (mcarray)

  std::string grpName = "multicomponent_array";
  const int NCOMP = 5;  // number of components in the array

  // Tests with non-zero- and with zero-sized arrays
  for(int SZ : {10, 0})
  {
    SLIC_INFO("Checking export of blueprint protocol for 'mcarray' with "
              << NCOMP << " components, each of size " << SZ);

    axom::sidre::DataStore ds1,ds2;
    conduit::Node node1,node2;

    // Simple lambda to get/check expected value at index (i,j)
    auto expValue = [=](int i, int j) {
                      return j*SZ + (SZ-i);
                    };
    // Simple lambda to get the name of the i^th component group
    auto viewName = [=](int j) {
                      return std::string(1,'a' + j); // i.e. 'a', 'b', 'c', ...
                    };

    // Initialize datastore, export to conduit and verify 'mcarray
    {
      auto* root = ds1.getRoot();
      auto* grp  = root->createGroup(grpName);
      for(int j=0 ; j<NCOMP ; ++j)
      {
        auto* view = grp->createViewAndAllocate(viewName(j), INT32_ID, SZ);

        int* sidre_data = view->getData();
        EXPECT_NE(nullptr, sidre_data);
        for(int i=0 ; i<SZ ; ++i)
        {
          sidre_data[i] = expValue(i,j);
        }
      }

      // Convert to conduit
      bool success = root->createNativeLayout(node1);
      EXPECT_TRUE(success);

      std::cout<<"Conduit native layout of datastore after exporting: \n";
      node1.print();
      std::cout << std::endl;

      // Call blueprint::verify to check that it complies with 'mcarray' schema
      conduit::Node info;
      bool verified =
        conduit::blueprint::verify("mcarray", node1[grpName], info);
      EXPECT_TRUE(verified);

      std::cout<<"Blueprint verify info for 'mcarray' protocol on node 1: \n";
      info.print();
      std::cout << std::endl;
    }

    // Import into sidre, export again to conduit, then verify blueprint
    {
      // Import into sidre
      auto* root = ds1.getRoot();
      bool success = root->importConduitTree(node1);
      EXPECT_TRUE(success);

      // export back to conduit
      success = root->createNativeLayout(node2);
      EXPECT_TRUE(success);

      std::cout<<"Conduit layout after export-import-export from Sidre: \n";
      node2.print();
      std::cout << std::endl;

      // verify the mcarray protocol
      conduit::Node info;
      bool verified =
        conduit::blueprint::verify("mcarray", node2[grpName], info);
      EXPECT_TRUE(verified);

      std::cout<<"Blueprint verify info for 'mcarray' protocol on node2: \n";
      info.print();
      std::cout << std::endl;
    }

    // Check the results
    {
      // Verify the two blueprints are compatible
      bool compat = node1.compatible(node2);
      EXPECT_TRUE(compat);

      // Check the data in the final conduit node
      conduit::Node& mcarr = node2[grpName];
      for(int j=0 ; j< NCOMP ; ++j)
      {
        std::string name = viewName(j);

        EXPECT_TRUE(mcarr.has_child(name));
        EXPECT_EQ(SZ, mcarr[name].dtype().number_of_elements());

        int* arr = mcarr[name].as_int_ptr();
        EXPECT_NE(nullptr, arr);

        for(int i=0 ; i<SZ ; ++i)
        {
          EXPECT_EQ( expValue(i,j), arr[i]);
        }
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger, finalized when
                          // exiting main scope
  axom::slic::setLoggingMsgLevel( axom::slic::message::Debug);

  result = RUN_ALL_TESTS();

  return result;
}
