/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "sidre/sidre.hpp"

using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;
using asctoolkit::sidre::INT_ID;

using namespace conduit;

#if 0
//------------------------------------------------------------------------------

TEST(sidre_view,create_views)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();

  DataView * dv_0 = root->createViewAndAllocate("field0",
                                                asctoolkit::sidre::INT_ID, 1);
  DataView * dv_1 = root->createViewAndAllocate("field1",
                                                asctoolkit::sidre::INT_ID, 1);


  DataBuffer * db_0 = dv_0->getBuffer();
  DataBuffer * db_1 = dv_1->getBuffer();

  EXPECT_EQ(db_0->getIndex(), 0);
  EXPECT_EQ(db_1->getIndex(), 1);
  delete ds;
}
#endif

//------------------------------------------------------------------------------

TEST(sidre_view,create_view_from_path)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();

  // Verify create works when groups must be created on demand.
  DataView * baz = root->createView("foo/bar/baz");
  // Groups should have been created.
  EXPECT_TRUE( root->hasGroup("foo") );
  EXPECT_TRUE( root->getGroup("foo")->hasGroup("bar") );

  DataGroup * bar = root->getGroup("foo")->getGroup("bar");
  EXPECT_TRUE( bar->hasView("baz") );
  EXPECT_EQ( bar->getView("baz"), baz );

  (void) baz;

  delete ds;

#if 0
  ds = new DataStore();
  root = ds->getRoot();

  // Verify create works when groups already exist.
  baz = root->createView("foo/bar/baz");
  EXPECT_TRUE( baz != ATK_NULLPTR );

  EXPECT_TRUE( root->hasGroup("foo") );

  std::cerr << "HERE2" << std::endl;

  foo = root->getGroup("foo");
  EXPECT_TRUE( foo->hasGroup("bar") );

  std::cerr << "HERE3" << std::endl;

  bar = foo->getGroup("bar");
  EXPECT_TRUE( bar->hasView("baz"));
  EXPECT_EQ( baz, bar->getView("baz"));

  std::cerr << "HERE4" << std::endl;

  delete ds;
#endif
}

//------------------------------------------------------------------------------

TEST(sidre_view,scalar_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  asctoolkit::sidre::SidreLength dims[2];

  // integer scalar
  DataView * i0view = root->createView("i0")->setScalar(1);
  EXPECT_EQ(i0view->getNumElements(), 1u);
  EXPECT_EQ(i0view->getNumDimensions(), 1);
  EXPECT_TRUE(i0view->getShape(1, dims) == 1 && dims[0] == 1);

  // Should report an error
  i0view->apply(INT_ID, 1);

  // string
  DataView * s0view = root->createView("s0")->setString("I am a string");
  EXPECT_EQ(s0view->getNumElements(), 14u);
  EXPECT_EQ(s0view->getNumDimensions(), 1);
  EXPECT_TRUE(s0view->getShape(1, dims) == 1 && dims[0] == 14);

  // using convenience functions
  // integer scalar
  DataView * i1view = root->createViewScalar("i1", 1);
  EXPECT_EQ(i1view->getNumElements(), 1u);
  EXPECT_EQ(i1view->getNumDimensions(), 1);
  EXPECT_TRUE(i1view->getShape(1, dims) == 1 && dims[0] == 1);

  // string
  DataView * s1view = root->createViewString("s1", "I am a string");
  EXPECT_EQ(s1view->getNumElements(), 14u);
  EXPECT_EQ(s1view->getNumDimensions(), 1);
  EXPECT_TRUE(s1view->getShape(1, dims) == 1 && dims[0] == 14);

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view,int_buffer_from_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataView * dv = root->createViewAndAllocate("u0", DataType::c_int(10));

  EXPECT_EQ(dv->getTypeID(), asctoolkit::sidre::INT_ID);
  int * data_ptr = dv->getData();

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  dv->print();

  EXPECT_EQ(dv->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>( sizeof(int) * 10) );
  delete ds;

}

//------------------------------------------------------------------------------

TEST(sidre_view,int_buffer_from_view_conduit_value)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataView * dv = root->createViewAndAllocate("u0", DataType::c_int(10));
  int * data_ptr = dv->getData();

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  dv->print();

  EXPECT_EQ(dv->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>( sizeof(int) * 10) );
  delete ds;

}

//------------------------------------------------------------------------------

TEST(sidre_view,int_array_strided_views)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataBuffer * dbuff = ds->createBuffer(asctoolkit::sidre::INT_ID, 10);

  dbuff->allocate();
  int * data_ptr = static_cast<int *>(dbuff->getData());

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i;
  }

  dbuff->print();

  EXPECT_EQ(dbuff->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int) * 10));

  DataView * dv_e = root->createView("even",dbuff);
  DataView * dv_o = root->createView("odd",dbuff);
  EXPECT_EQ(dbuff->getNumViews(), 2);

  // c_int(num_elems, offset [in bytes], stride [in bytes])
  dv_e->apply(DataType::c_int(5,0,8));

  // c_int(num_elems, offset [in bytes], stride [in bytes])
  dv_o->apply(DataType::c_int(5,4,8));

  dv_e->print();
  dv_o->print();

  // Check base pointer case:
  int * v_e_ptr = dv_e->getData();
  int * v_o_ptr = dv_o->getData();
  for(int i=0 ; i<10 ; i += 2)
  {
    std::cout << "idx:" <<  i
              << " e:" << v_e_ptr[i]
              << " o:" << v_o_ptr[i]
              << " em:" << v_e_ptr[i]  % 2
              << " om:" << v_o_ptr[i]  % 2
              << std::endl;

    EXPECT_EQ(v_e_ptr[i] % 2, 0);
    EXPECT_EQ(v_o_ptr[i] % 2, 1);
  }

  // Check Conduit mem-map struct case:
  int_array dv_e_ptr = dv_e->getData();
  int_array dv_o_ptr = dv_o->getData();
  for(int i=0 ; i<5 ; ++i)
  {
    std::cout << "idx:" <<  i
              << " e:" << dv_e_ptr[i]
              << " o:" << dv_o_ptr[i]
              << " em:" << dv_e_ptr[i]  % 2
              << " om:" << dv_o_ptr[i]  % 2
              << std::endl;

    EXPECT_EQ(dv_e_ptr[i] % 2, 0);
    EXPECT_EQ(dv_o_ptr[i] % 2, 1);
  }

  // Run similar test to above with different view apply method
  DataView * dv_e1 = root->createView("even1",dbuff);
  DataView * dv_o1 = root->createView("odd1",dbuff);
  EXPECT_EQ(dbuff->getNumViews(), 4);

  // (num_elems, offset [in # elems], stride [in # elems])
  dv_e1->apply(asctoolkit::sidre::INT_ID, 5,0,2);

  // (num_elems, offset [in # elems], stride [in # elems])
  dv_o1->apply(asctoolkit::sidre::INT_ID, 5,1,2);

  dv_e1->print();
  dv_o1->print();

  // Check base pointer case:
  int * v_e1_ptr = dv_e1->getData();
  int * v_o1_ptr = dv_o1->getData();
  for(int i=0 ; i<10 ; i += 2)
  {
    std::cout << "idx:" <<  i
              << " e1:" << v_e1_ptr[i]
              << " oj:" << v_o1_ptr[i]
              << " em1:" << v_e1_ptr[i]  % 2
              << " om1:" << v_o1_ptr[i]  % 2
              << std::endl;

    EXPECT_EQ(v_e1_ptr[i], v_e_ptr[i]);
    EXPECT_EQ(v_o1_ptr[i], v_o_ptr[i]);
  }

  // Check Conduit mem-map struct case:
  int_array dv_e1_ptr = dv_e1->getData();
  int_array dv_o1_ptr = dv_o1->getData();
  for(int i=0 ; i<5 ; i++)
  {
    std::cout << "idx:" <<  i
              << " e1:" << dv_e1_ptr[i]
              << " o1:" << dv_o1_ptr[i]
              << " em1:" << dv_e1_ptr[i]  % 2
              << " om1:" << dv_o1_ptr[i]  % 2
              << std::endl;

    EXPECT_EQ(dv_e1_ptr[i], dv_e_ptr[i]);
    EXPECT_EQ(dv_o1_ptr[i], dv_o_ptr[i]);
  }


  ds->print();
  delete ds;

}

//------------------------------------------------------------------------------

TEST(sidre_view,int_array_depth_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  const asctoolkit::sidre::SidreLength depth_nelems = 10;
  DataBuffer * dbuff = ds->createBuffer(asctoolkit::sidre::INT_ID,
                                        4 * depth_nelems);

  // Allocate buffer to hold data for 4 "depth" views
  dbuff->allocate();
  int * data_ptr = static_cast<int *>(dbuff->getData());

  for(size_t i = 0 ; i < 4 * depth_nelems ; ++i)
  {
    data_ptr[i] = i / depth_nelems;
  }

  dbuff->print();

  EXPECT_EQ(dbuff->getNumElements(), 4 * depth_nelems);

  // create 4 "depth" views and apply offsets into buffer
  DataView * views[4];
  std::string view_names[4] = { "depth_0", "depth_1", "depth_2", "depth_3" };

  for (int id = 0 ; id < 2 ; ++id)
  {
    views[id] = root->createView(view_names[id], dbuff)->apply(depth_nelems,
                                                               id*depth_nelems);
  }
  //
  // call path including type
  for (int id = 2 ; id < 4 ; ++id)
  {
    views[id] = root->createView(view_names[id], dbuff)->apply(
      asctoolkit::sidre::INT_ID,
      depth_nelems,
      id*depth_nelems);
  }
  EXPECT_EQ(dbuff->getNumViews(), 4);

  // print depth views...
  for (int id = 0 ; id < 4 ; ++id)
  {
    views[id]->print();
  }

  // check values in depth views...
  for (int id = 0 ; id < 4 ; ++id)
  {
    int * dv_ptr = views[id]->getData();
    for (size_t i = 0 ; i < depth_nelems ; ++i)
    {
      EXPECT_EQ(dv_ptr[i], id);
    }
  }

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view,int_array_view_attach_buffer)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  const asctoolkit::sidre::SidreLength field_nelems = 10;

  // create 2 "field" views with type and # elems
  asctoolkit::sidre::SidreLength elem_count = 0;
  DataView * field0 = root->createView("field0",
                                       asctoolkit::sidre::INT_ID, field_nelems);
  elem_count += field0->getNumElements();
  DataView * field1 = root->createView("field1",
                                       asctoolkit::sidre::INT_ID, field_nelems);
  elem_count += field1->getNumElements();
  EXPECT_EQ(elem_count, 2 * field_nelems);

  // create buffer to hold data for fields and allocate
  DataBuffer * dbuff = ds->createBuffer()->allocate(asctoolkit::sidre::INT_ID,
                                                    elem_count);
  EXPECT_EQ(dbuff->getNumElements(), elem_count);

  // Initilize buffer data for testing below.
  int * b_ptr = dbuff->getData();
  for(asctoolkit::sidre::SidreLength i = 0 ; i < elem_count ; ++i)
  {
    b_ptr[i] = i / field_nelems;
  }

  dbuff->print();

  // attach field views to buffer and apply offsets into buffer
  field0->attachBuffer(dbuff)->apply(field_nelems, 0 * field_nelems);
  field1->attachBuffer(dbuff)->apply(field_nelems, 1 * field_nelems);
  EXPECT_EQ(dbuff->getNumViews(), 2);

  // print field views...
  field0->print();
  field1->print();

  // check values in field views...
  int * f0_ptr = field0->getData();
  for (size_t i = 0 ; i < field_nelems ; ++i)
  {
    EXPECT_EQ(f0_ptr[i], 0);
  }
  int * f1_ptr = field1->getData();
  for (size_t i = 0 ; i < field_nelems ; ++i)
  {
    EXPECT_EQ(f1_ptr[i], 1);
  }

  ds->print();
  delete ds;
}


//------------------------------------------------------------------------------

TEST(sidre_view,int_array_multi_view_resize)
{
  ///
  /// This example creates a 4 * 10 buffer of ints,
  /// and 4 views that point the 4 sections of 10 ints
  ///
  /// We then create a new buffer to support 4*12 ints
  /// and 4 views that point into them
  ///
  /// after this we use the old buffers to copy the values
  /// into the new views
  ///

  // create our main data store
  DataStore * ds = new DataStore();
  // get access to our root data Group
  DataGroup * root = ds->getRoot();

  // create a group to hold the "old" or data we want to copy
  DataGroup * r_old = root->createGroup("r_old");
  // create a view to hold the base buffer and allocate
  DataView * base_old = r_old->createViewAndAllocate("base_data",
                                                     DataType::c_int(40));

  // we will create 4 sub views of this array
  int * data_ptr = base_old->getData();


  // init the buff with values that align with the
  // 4 subsections.
  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = 1;
  }
  for(int i=10 ; i<20 ; i++)
  {
    data_ptr[i] = 2;
  }
  for(int i=20 ; i<30 ; i++)
  {
    data_ptr[i] = 3;
  }
  for(int i=30 ; i<40 ; i++)
  {
    data_ptr[i] = 4;
  }


  /// setup our 4 views
  DataBuffer * buff_old = base_old->getBuffer();
  buff_old->print();
  DataView * r0_old = r_old->createView("r0",buff_old);
  DataView * r1_old = r_old->createView("r1",buff_old);
  DataView * r2_old = r_old->createView("r2",buff_old);
  DataView * r3_old = r_old->createView("r3",buff_old);

  // each view is offset by 10 * the # of bytes in a int
  // c_int(num_elems, offset)
  index_t offset =0;
  r0_old->apply(DataType::c_int(10,offset));

  offset += sizeof(int) * 10;
  r1_old->apply(DataType::c_int(10,offset));

  offset += sizeof(int) * 10;
  r2_old->apply(DataType::c_int(10,offset));

  offset += sizeof(int) * 10;
  r3_old->apply(DataType::c_int(10,offset));

  /// check that our views actually point to the expected data
  //
  int * r0_ptr = r0_old->getData();
  for(int i=0 ; i<10 ; i++)
  {
    EXPECT_EQ(r0_ptr[i], 1);
    // check pointer relation
    EXPECT_EQ(&r0_ptr[i], &data_ptr[i]);
  }

  int * r3_ptr = r3_old->getData();
  for(int i=0 ; i<10 ; i++)
  {
    EXPECT_EQ(r3_ptr[i], 4);
    // check pointer relation
    EXPECT_EQ(&r3_ptr[i], &data_ptr[i+30]);
  }

  // create a group to hold the "old" or data we want to copy into
  DataGroup * r_new = root->createGroup("r_new");
  // create a view to hold the base buffer and allocate
  DataView * base_new = r_new->createViewAndAllocate("base_data",
                                                     DataType::c_int(4 * 12));

  int * base_new_data = base_new->getData();
  for (int i = 0 ; i < 4 * 12 ; ++i)
  {
    base_new_data[i] = 0;
  }

  DataBuffer * buff_new = base_new->getBuffer();
  buff_new->print();

  // create the 4 sub views of this array
  DataView * r0_new = r_new->createView("r0",buff_new);
  DataView * r1_new = r_new->createView("r1",buff_new);
  DataView * r2_new = r_new->createView("r2",buff_new);
  DataView * r3_new = r_new->createView("r3",buff_new);

  // apply views to r0,r1,r2,r3
  // each view is offset by 12 * the # of bytes in a int

  // c_int(num_elems, offset)
  offset =0;
  r0_new->apply(DataType::c_int(12,offset));

  offset += sizeof(int) * 12;
  r1_new->apply(DataType::c_int(12,offset));

  offset += sizeof(int) * 12;
  r2_new->apply(DataType::c_int(12,offset));

  offset += sizeof(int) * 12;
  r3_new->apply(DataType::c_int(12,offset));

  /// update r2 as an example first
  buff_new->print();
  r2_new->print();

  /// copy the subset of value
  r2_new->getNode().update(r2_old->getNode());
  r2_new->getNode().print();
  buff_new->print();


  /// check pointer values
  int * r2_new_ptr = r2_new->getData();

  for(int i=0 ; i<10 ; i++)
  {
    EXPECT_EQ(r2_new_ptr[i], 3);
  }

  for(int i=10 ; i<12 ; i++)
  {
    EXPECT_EQ(r2_new_ptr[i], 0);     // assumes zero-ed alloc
  }


  /// update the other views
  r0_new->getNode().update(r0_old->getNode());
  r1_new->getNode().update(r1_old->getNode());
  r3_new->getNode().update(r3_old->getNode());

  buff_new->print();


  ds->print();
  delete ds;

}


//------------------------------------------------------------------------------

TEST(sidre_view,int_array_realloc)
{
  ///
  /// info
  ///

  // create our main data store
  DataStore * ds = new DataStore();
  // get access to our root data Group
  DataGroup * root = ds->getRoot();

  // create a view to hold the base buffer
  DataView * a1 = root->createViewAndAllocate("a1",DataType::c_float(5));
  DataView * a2 = root->createViewAndAllocate("a2",DataType::c_int(5));

  float * a1_ptr = a1->getData();
  int * a2_ptr = a2->getData();

  for(int i=0 ; i<5 ; i++)
  {
    a1_ptr[i] =  5.0;
    a2_ptr[i] = -5;
  }

  EXPECT_EQ(a1->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(float)*5));
  EXPECT_EQ(a2->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int)*5));


  a1->reallocate(DataType::c_float(10));
  a2->reallocate(DataType::c_int(15));

  a1_ptr = a1->getData();
  a2_ptr = a2->getData();

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(a1_ptr[i],5.0);
    EXPECT_EQ(a2_ptr[i],-5);
  }

  for(int i=5 ; i<10 ; i++)
  {
    a1_ptr[i] = 10.0;
    a2_ptr[i] = -10;
  }

  for(int i=10 ; i<15 ; i++)
  {
    a2_ptr[i] = -15;
  }

  EXPECT_EQ(a1->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(float)*10));
  EXPECT_EQ(a2->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int)*15));

  // Try some errors
  // XXX  a1->reallocate(DataType::c_int(20));

  ds->print();
  delete ds;

}

//------------------------------------------------------------------------------

TEST(sidre_view,simple_opaque)
{
  // create our main data store
  DataStore * ds = new DataStore();
  // get access to our root data Group
  DataGroup * root = ds->getRoot();
  int * src_data = new int[1];

  src_data[0] = 42;

  void * src_ptr = (void *)src_data;

  DataView * opq_view = root->createView("my_opaque", src_ptr);

  // External pointers are held in the view, should not have a buffer.
  EXPECT_EQ(ds->getNumBuffers(), 0u);

  EXPECT_TRUE(opq_view->isExternal());
  EXPECT_TRUE(!opq_view->isApplied());
  EXPECT_TRUE(opq_view->isOpaque());

  void * opq_ptr = opq_view->getVoidPtr();

  int * out_data = (int *)opq_ptr;
  EXPECT_EQ(opq_ptr,src_ptr);
  EXPECT_EQ(out_data[0],42);

  ds->print();
  delete ds;
  delete [] src_data;
}
