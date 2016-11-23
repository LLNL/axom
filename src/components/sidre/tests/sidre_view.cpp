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
using asctoolkit::sidre::SidreLength;
using asctoolkit::sidre::IndexType;
using asctoolkit::sidre::TypeID;
using asctoolkit::sidre::NO_TYPE_ID;
using asctoolkit::sidre::INT_ID;
using asctoolkit::sidre::DOUBLE_ID;
using asctoolkit::sidre::CHAR8_STR_ID;

using namespace conduit;

#define BLEN 10

// Match ViewBuffer states
enum State
{
  EMPTY,
  BUFFER,
  EXTERNAL,
  SCALAR,
  STRING,
  NOTYPE
};


// Check state of View based on booleans since m_state is private.

static State getState(DataView * view)
{
  if (view->isEmpty())
  {
    return EMPTY;
  }
  else if (view->hasBuffer())
  {
    return BUFFER;
  }
  else if (view->isExternal())
  {
    return EXTERNAL;
  }
  else if (view->isScalar())
  {
    return SCALAR;
  }
  else if (view->isString())
  {
    return STRING;
  }
  else
  {
    return NOTYPE;
  }
}

//
// Test various paths to creating a buffer.
// Assume all are int[len]
//
// Since this function is called in several places, any test failure
// will be ambiguous. To use, wrap in the EXPECT_TRUE macro.
//    EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
// Any failures will have two reports, one from this routine and one
// from the call site.
//
static bool checkViewValues(DataView * view,
                            State state,
                            bool isDescribed, bool isAllocated,
                            bool isApplied,
                            SidreLength len)
{
  bool rv = true;
  SidreLength dims[2];

  if (getState(view) != state)
  {
    EXPECT_EQ(getState(view), state);
    rv = false;
  }

  if (view->isDescribed() != isDescribed)
  {
    EXPECT_EQ(isDescribed, view->isDescribed());
    rv = false;
  }
  if (view->isAllocated() != isAllocated)
  {
    EXPECT_EQ(isAllocated, view->isAllocated());
    rv = false;
  }
  if (view->isApplied() != isApplied)
  {
    EXPECT_EQ(isApplied, view->isApplied());
    rv = false;
  }

  if (view->getNumElements() != len)
  {
    EXPECT_EQ(len, view->getNumElements());
    rv = false;
  }

  if (view->isDescribed())
  {
    if (view->getTypeID() != INT_ID)
    {
      EXPECT_EQ(INT_ID, view->getTypeID());
      rv = false;
    }
    if (view->getNumDimensions() != 1)
    {
      EXPECT_EQ(1, view->getNumDimensions());
      rv = false;
    }
    if (view->getShape(1, dims) != 1 || dims[0] != len)
    {
      EXPECT_TRUE(view->getShape(1, dims) == 1 && dims[0] == len);
      rv = false;
    }
    if (view->getTotalBytes() != static_cast<SidreLength>( sizeof(int) * len) )
    {
      EXPECT_EQ(view->getTotalBytes(),
                static_cast<SidreLength>( sizeof(int) * len) );
      rv = false;
    }
  }
  else
  {
    TypeID id = view->getTypeID();
    if (id != NO_TYPE_ID)
    {
      EXPECT_EQ(NO_TYPE_ID, id);
      rv = false;
    }
  }

  if (view->isApplied())
  {
    // Fill with data to help the print function
    int * data_ptr = view->getData();

    for(int i=0 ; i<len ; i++)
    {
      data_ptr[i] = i*i;
    }
  }

  //  view->print();

  return rv;
}

#if 0
//------------------------------------------------------------------------------

TEST(sidre_view,create_views)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();

  DataView * dv_0 = root->createViewAndAllocate("field0", INT_ID, 1);
  DataView * dv_1 = root->createViewAndAllocate("field1", INT_ID, 1);


  DataBuffer * db_0 = dv_0->getBuffer();
  DataBuffer * db_1 = dv_1->getBuffer();

  EXPECT_EQ(db_0->getIndex(), 0);
  EXPECT_EQ(db_1->getIndex(), 1);
  delete ds;
}
#endif

//------------------------------------------------------------------------------

TEST(sidre_view,get_path_name)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataView * v1 = root->createView("test/a/b/v1");
  DataView * v2 = root->createView("test/v2");
  DataView * v3 = root->createView("v3");

  EXPECT_EQ(std::string("v1"), v1->getName());
  EXPECT_EQ(std::string("test/a/b"), v1->getPath());
  EXPECT_EQ(std::string("test/a/b/v1"), v1->getPathName());

  EXPECT_EQ(std::string("v2") , v2->getName());
  EXPECT_EQ(std::string("test") , v2->getPath());
  EXPECT_EQ(std::string("test/v2") , v2->getPathName());

  EXPECT_EQ(std::string("v3") , v3->getName());
  EXPECT_EQ(std::string("") , v3->getPath());
  EXPECT_EQ(std::string("v3") , v3->getPathName());
}

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

static void checkScalarValues(DataView * view,
                              State state,
                              bool isDescribed, bool isAllocated,
                              bool isApplied,
                              TypeID type, SidreLength len)
{
  SidreLength dims[2];

  SCOPED_TRACE(view->getName());

  EXPECT_EQ(getState(view), state);

  EXPECT_EQ(view->isDescribed(), isDescribed);
  EXPECT_EQ(view->isAllocated(), isAllocated);
  EXPECT_EQ(view->isApplied(), isApplied);

  EXPECT_EQ(view->getTypeID(), type);
  EXPECT_EQ(view->getNumElements(), len);
  EXPECT_EQ(view->getNumDimensions(), 1);
  EXPECT_TRUE(view->getShape(1, dims) == 1 && dims[0] == len);
}

TEST(sidre_view,scalar_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  int i;
  const char * s;

  DataView * i0view = root->createView("i0")->setScalar(1);
  checkScalarValues(i0view, SCALAR, true, true, true, INT_ID, 1);
  i = i0view->getScalar();
  EXPECT_EQ( 1, i);

  DataView * i1view = root->createViewScalar("i1", 2);
  checkScalarValues(i1view, SCALAR, true, true, true, INT_ID, 1);
  i = i1view->getScalar();
  EXPECT_EQ( 2, i);

  DataView * s0view = root->createView("s0")->setString("I am a string");
  checkScalarValues(s0view, STRING, true, true, true, CHAR8_STR_ID, 14);
  s = s0view->getString();
  EXPECT_TRUE( strcmp(s, "I am a string") == 0);

  DataView * s1view = root->createViewString("s1", "I too am a string");
  checkScalarValues(s1view, STRING, true, true, true, CHAR8_STR_ID, 18);
  s = s1view->getString();
  EXPECT_TRUE( strcmp(s, "I too am a string") == 0);

  // Check illegal operations
  i0view->apply(INT_ID, 1);
  i0view->allocate();
  i0view->deallocate();

  s0view->apply(INT_ID, 1);
  s0view->allocate();
  s0view->deallocate();

  DataView * empty = root->createView("empty");
#if 0
  try
  {
    int * j = empty->getScalar();
    //int j = empty->getScalar();
    EXPECT_EQ(0, *j);
  }
  catch ( conduit::Error e)
  {}
#endif
  const char * svalue = empty->getString();
  EXPECT_EQ(NULL, svalue);

  delete ds;
}

//------------------------------------------------------------------------------

// Most tests deallocate via the DataStore destructor
// This is an explicit deallocate test

TEST(sidre_view,dealloc)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataBuffer * dbuff;
  DataView * dv;

  //----------  EMPTY(F,F,F)
  dv = root->createView("e1");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  //----------  EMPTY(T,F,F)
  dv = root->createView("e2", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  //----------  BUFFER(F,F,F)
  dv = root->createView("b1");
  dbuff = ds->createBuffer();
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));
  EXPECT_FALSE(dv->getBuffer()->isAllocated());

  //---------- BUFFER(T,F,F)
  dv = root->createView("b2", INT_ID, BLEN);
  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));
  EXPECT_FALSE(dv->getBuffer()->isAllocated());

  //---------- BUFFER(T,T,T)
  dv = root->createView("b3", INT_ID, BLEN);
  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));
  EXPECT_FALSE(dv->getBuffer()->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------

// allocate/reallocate with zero items results in allocated (yet zero).

TEST(sidre_view,alloc_zero_items)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataView * dv;

  // Allocate zero items
  dv = root->createView("z0");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(INT_ID, 0);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, false, 0));
  EXPECT_TRUE(dv->getBuffer()->isAllocated());

  // Reallocate zero items
  dv = root->createView("z1");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));
  EXPECT_TRUE(dv->getBuffer()->isAllocated());
  dv->reallocate(0);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, false, 0));
  EXPECT_TRUE(dv->getBuffer()->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------

// Test allocate, reallocate, and deallocate when there are multiple views
// attached to a buffer.  All operations should be no-ops.

TEST(sidre_view,alloc_and_dealloc_multiview)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataBuffer * dbuff;
  DataView * dv1, * dv2;
  void * baddr;

  //---------- allocate
  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  EXPECT_FALSE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  baddr = dbuff->getVoidPtr();

  dv1 = root->createView("dv1alloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  dv2 = root->createView("dv2alloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 2);

  dv1->allocate(INT_ID, BLEN+10);
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));

  // buffer is unchanged
  EXPECT_FALSE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  EXPECT_EQ(dbuff->getVoidPtr(), baddr);

  //---------- reallocate
  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  baddr = dbuff->getVoidPtr();

  dv1 = root->createView("dv1realloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  dv2 = root->createView("dv2realloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 2);

  dv1->reallocate(BLEN+10);
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));

  // buffer is unchanged
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  EXPECT_EQ(dbuff->getVoidPtr(), baddr);

  //---------- deallocate
  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  baddr = dbuff->getVoidPtr();

  dv1 = root->createView("dv1dealloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  dv2 = root->createView("dv2dealloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 2);

  dv1->deallocate();
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));

  // buffer is unchanged
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  EXPECT_EQ(dbuff->getVoidPtr(), baddr);

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view,int_alloc_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataView * dv;
  long shape[] = { BLEN };

  dv = root->createView("u0");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));
#if 0
  dv = root->createView("u1");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(INT_ID, 1, shape);  // XXX - missing overload
  EXPECT_TRUE(checkViewValues(dv, BUFFEER, true, true, true, BLEN));
#endif
  dv = root->createView("u2");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(DataType::c_int(BLEN));
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));


  dv = root->createView("v0", INT_ID, 10);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createView("v1", INT_ID, 1, shape);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createView("v2", DataType::c_int(BLEN));
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));


  dv = root->createViewAndAllocate("a0", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createViewAndAllocate("a1", INT_ID, 1, shape);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createViewAndAllocate("a2", DataType::c_int(BLEN));
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  delete ds;
}

//------------------------------------------------------------------------------

// Test some association/state transitions.
// All of these tests have only one view per buffer.

TEST(sidre_view,int_buffer_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataBuffer * dbuff, * otherbuffer;
  DataView * dv;
  IndexType bindex;
  //  long shape[] = { BLEN };

  //     view                        buffer
  //                  undescribed  described  allocated
  //     undescribed      1           2           3
  //     described        4           5           6
  //
  //  buffer description incomptable with view description
  //                               described  allocated
  //     described                    7           8


  //---------- 1
  // Attach undescribed buffer to undescribed view
  dv = root->createView("u1");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  // no-op, attach NULL buffer to EMPTY view
  dv->attachBuffer(NULL);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dbuff = ds->createBuffer();
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));
  EXPECT_EQ(dv->getBuffer(), dbuff);  // sanity check

  dv->allocate();  // no-op, no description
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));

  dv->allocate(INT_ID, 10);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv->reallocate(BLEN+5);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN+5));

  // After deallocate, description is intact.
  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN+5));

  //---------- 2
  // Attach described buffer to undescribed view
  dv = root->createView("u2");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));

  //---------- 3
  // Attach allocated buffer to undescribed view
  dv = root->createView("u3");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  bindex = dbuff->getIndex();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 1);

  // no-op, attaching a buffer to a view which already has a buffer
  otherbuffer = ds->createBuffer();
  dv->attachBuffer(otherbuffer);
  EXPECT_EQ(dv->getBuffer(), dbuff); // same buffer as before
  EXPECT_EQ(otherbuffer->getNumViews(), 0);

  // The current attached buffer will be destroyed.
  dv->attachBuffer(NULL);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  EXPECT_TRUE(ds->getBuffer(bindex) == NULL);

  //---------- 4
  // Attach undescribed buffer to described view
  dv = root->createView("u4", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dbuff = ds->createBuffer();
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  //---------- 5
  // Attach described buffer to described view
  dv = root->createView("u5", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  //---------- 6
  // Attach allocated buffer to described view
  dv = root->createView("u6", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  // Deallocate the buffer which will update the view.
  dbuff->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  // Allocate the buffer via the view.
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));
  EXPECT_TRUE(dbuff->isAllocated());

  //---------- 7
  // Attach incompatable described buffer to described view
  dv = root->createView("u7", INT_ID, BLEN+5);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN+5));

  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN+5));

  //---------- 8
  // Attach incompatable allocated buffer to described view
  dv = root->createView("u8", INT_ID, BLEN+5);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN+5));

  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, false, BLEN+5));  // XXX - how is isAllocated useful

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view,int_array_strided_views)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataBuffer * dbuff = ds->createBuffer(INT_ID, 10);

  dbuff->allocate();
  int * data_ptr = dbuff->getData();

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i;
  }

  dbuff->print();

  EXPECT_EQ(dbuff->getTotalBytes(), static_cast<SidreLength>(sizeof(int) * 10));

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
  dv_e1->apply(INT_ID, 5,0,2);

  // (num_elems, offset [in # elems], stride [in # elems])
  dv_o1->apply(INT_ID, 5,1,2);

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

  // Delete buffer and make sure views are no longer allocated.
  dbuff->deallocate();
  EXPECT_FALSE(dv_e->isAllocated());
  EXPECT_FALSE(dv_o->isAllocated());
  EXPECT_FALSE(dv_e1->isAllocated());
  EXPECT_FALSE(dv_o1->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view,int_array_depth_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  const SidreLength depth_nelems = 10;
  DataBuffer * dbuff = ds->createBuffer(INT_ID, 4 * depth_nelems);

  // Allocate buffer to hold data for 4 "depth" views
  dbuff->allocate();
  int * data_ptr = dbuff->getData();

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
    views[id] = root->createView(view_names[id], dbuff)
                ->apply(INT_ID, depth_nelems, id*depth_nelems);
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

TEST(sidre_view,view_offset_and_stride)
{
    DataStore * ds = new DataStore();
    DataGroup * root = ds->getRoot();

    const SidreLength nelems = 15;
    DataBuffer * dbuff = ds->createBuffer(DOUBLE_ID, nelems)
                           ->allocate();

    double * data_ptr = dbuff->getData();
    for(int i=0; i< nelems; ++i)
    {
        data_ptr[i] = 1.01 * i;
    }

    // Array layout:
    //   [-,a,b,c,a,b,c,a,b,c,d,d,d,d,d]
    const int NUM_GROUPS   = 4;
    int size[NUM_GROUPS]   = {3, 3, 3, 5};
    int stride[NUM_GROUPS] = {3, 3, 3, 1};
    int offset[NUM_GROUPS] = {1, 2, 3, 10};
    std::string names[NUM_GROUPS] = {"a","b","c","d"};


    // -- Test offsets and strides on buffer-based views

    DataGroup* bufferGroup = root->createGroup("buffer");

    // Create the views off of the buffer
    for(int i=0; i < NUM_GROUPS; ++i)
    {
        bufferGroup->createView(names[i], dbuff)
                   ->apply(size[i], offset[i], stride[i]);
    }

    // Test the sizes, offsets and strides of the views
    for(int i=0; i < NUM_GROUPS; ++i)
    {
        DataView* view = bufferGroup->getView(names[i]);
        EXPECT_EQ(size[i],   view->getNumElements());
        EXPECT_EQ(offset[i], view->getOffset());
        EXPECT_EQ(stride[i], view->getStride());
    }

    // test the offsets and strides applied to the data
    // Note: DataView::getData() already incorporates the offset
    //       of the view, so we are comparing against
    //       the offsets applied to the raw buffer
    for(int i=0; i < NUM_GROUPS; ++i)
    {
        DataView* view = bufferGroup->getView(names[i]);
        EXPECT_EQ(view->getData<double*>(), data_ptr + view->getOffset());
        EXPECT_EQ(stride[i], view->getStride());

    }


    // -- Test offsets and strudes on external pointer based views
    DataGroup* extGroup = root->createGroup("ext");

    // Create the views off of external data pointer
    for(int i=0; i < NUM_GROUPS; ++i)
    {
        extGroup->createView(names[i], data_ptr)
                ->apply(DOUBLE_ID, size[i], offset[i], stride[i]);
    }

    // Test the sizes, offsets and strides of the views
    for(int i=0; i < NUM_GROUPS; ++i)
    {
        DataView* view = extGroup->getView(names[i]);
        EXPECT_EQ(size[i],   view->getNumElements());
        EXPECT_EQ(offset[i], view->getOffset());
        EXPECT_EQ(stride[i], view->getStride());
    }

    // test the offsets and strides applied to the data
    for(int i=0; i < NUM_GROUPS; ++i)
    {
        DataView* view = extGroup->getView(names[i]);
        EXPECT_EQ(view->getData<double*>(), data_ptr + view->getOffset());
        EXPECT_EQ(stride[i], view->getStride());

    }

    // -- Test offset and stride on the other view types:
    //          string, scalar, empty, opaque
    DataGroup* othersGroup = root->createGroup("others");

    typedef std::vector<DataView*> ViewVec;
    ViewVec views;
    asctoolkit::sidre::detail::sidre_uint8  ui8  = 3;
    asctoolkit::sidre::detail::sidre_uint16 ui16 = 4;
    asctoolkit::sidre::detail::sidre_uint32 ui32 = 5;
    asctoolkit::sidre::detail::sidre_uint64 ui64 = 6;
    asctoolkit::sidre::detail::sidre_int8   i8   = -3;
    asctoolkit::sidre::detail::sidre_int16 i16   = -4;
    asctoolkit::sidre::detail::sidre_int32 i32   = -5;
    asctoolkit::sidre::detail::sidre_int64 i64   = -6;
    asctoolkit::sidre::detail::sidre_float32 f32 = 7.7f;
    asctoolkit::sidre::detail::sidre_float64 f64 = 8.8;


    views.push_back( othersGroup->createView("key_empty"));
    views.push_back( othersGroup->createView("key_opaque", data_ptr)); // not described
    views.push_back( othersGroup->createViewString("key_string", "string_value"));

    views.push_back( othersGroup->createViewScalar("key_uint8",   ui8));
    views.push_back( othersGroup->createViewScalar("key_uint16",  ui16));
    views.push_back( othersGroup->createViewScalar("key_uint32",  ui32));
    views.push_back( othersGroup->createViewScalar("key_uint64",  ui64));

    views.push_back( othersGroup->createViewScalar("key_int8",    i8));
    views.push_back( othersGroup->createViewScalar("key_int16",   i16));
    views.push_back( othersGroup->createViewScalar("key_int32",   i32));
    views.push_back( othersGroup->createViewScalar("key_int64",   i64));

    views.push_back( othersGroup->createViewScalar("key_float32", f32));
    views.push_back( othersGroup->createViewScalar("key_float64", f64));


    for(ViewVec::iterator it=views.begin(); it != views.end(); ++it)
    {
        DataView* view = *it;
        EXPECT_EQ(0, view->getOffset());
        EXPECT_EQ(1, view->getStride());
    }


    // -- cleanup

    ds->print();
    delete ds;
}



//------------------------------------------------------------------------------

TEST(sidre_view,int_array_view_attach_buffer)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  const SidreLength field_nelems = 10;

  // create 2 "field" views with type and # elems
  SidreLength elem_count = 0;
  DataView * field0 = root->createView("field0", INT_ID, field_nelems);
  elem_count += field0->getNumElements();
  DataView * field1 = root->createView("field1",INT_ID, field_nelems);
  elem_count += field1->getNumElements();
  EXPECT_EQ(elem_count, 2 * field_nelems);

  // create buffer to hold data for fields and allocate
  DataBuffer * dbuff = ds->createBuffer()->allocate(INT_ID, elem_count);
  EXPECT_EQ(dbuff->getNumElements(), elem_count);

  // Initialize buffer data for testing below.
  int * b_ptr = dbuff->getData();
  for(SidreLength i = 0 ; i < elem_count ; ++i)
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

  EXPECT_EQ(a1->getTotalBytes(), static_cast<SidreLength>(sizeof(float)*5));
  EXPECT_EQ(a2->getTotalBytes(), static_cast<SidreLength>(sizeof(int)*5));


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

  EXPECT_EQ(a1->getTotalBytes(), static_cast<SidreLength>(sizeof(float)*10));
  EXPECT_EQ(a2->getTotalBytes(), static_cast<SidreLength>(sizeof(int)*15));

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

//------------------------------------------------------------------------------

TEST(sidre_datastore,destroy_buffer)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataBuffer * dbuff1 = ds->createBuffer()->allocate(INT_ID, BLEN);
  DataView * view1a = root->createView("view1a", INT_ID, BLEN)
                      ->attachBuffer(dbuff1);
  DataView * view1b = root->createView("view1b", INT_ID, BLEN)
                      ->attachBuffer(dbuff1);

  EXPECT_TRUE(checkViewValues(view1a, BUFFER, true, true, true, BLEN));
  EXPECT_TRUE(checkViewValues(view1b, BUFFER, true, true, true, BLEN));

  // destroyBuffer will detach from views.
  IndexType bindex = dbuff1->getIndex();
  ds->destroyBuffer(bindex);
  EXPECT_TRUE(ds->getBuffer(bindex) == NULL);

  // views no longer have buffers (but retain descriptions)
  EXPECT_TRUE(checkViewValues(view1a, EMPTY, true, false, false, BLEN));
  EXPECT_TRUE(checkViewValues(view1b, EMPTY, true, false, false, BLEN));

  EXPECT_FALSE(view1a->hasBuffer());
  EXPECT_FALSE(view1b->hasBuffer());

  delete ds;
}


//------------------------------------------------------------------------------
TEST(sidre_view,value_from_uninited_view)
{
    DataStore ds;
    DataView * view = ds.getRoot()->createView("empty");

    // check getScalar
    int val = view->getScalar();
    EXPECT_EQ(val,0);

    // check getArray
    int *aval_ptr =    view->getArray();
    EXPECT_TRUE( aval_ptr == NULL );

    int aval = view->getArray();
    EXPECT_EQ(aval,0);
    
    // check getData
    int *dval_ptr =    view->getData();
    EXPECT_TRUE( dval_ptr == NULL );

    int dval = view->getData();
    EXPECT_EQ(dval,0);
    
}
