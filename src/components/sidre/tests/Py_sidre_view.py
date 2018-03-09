###############################################################################
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
###############################################################################

import unittest
import sidre
#from sidre import InvalidIndex, nameIsValid

class SidreView(unittest.TestCase):

#------------------------------------------------------------------------------

    def test_create_views(self):
        ds = sidre.DataStore()
        root = ds.getRoot()
        self.assertEqual(type(root), sidre.Group)
        
        dv_0 = root.createViewAndAllocate("field0", sidre.INT_ID, 1)
        dv_1 = root.createViewAndAllocate("field1", sidre.INT_ID, 1)

        db_0 = dv_0.getBuffer()
        db_1 = dv_1.getBuffer()

        self.assertEqual(db_0.getIndex(), 0)
        self.assertEqual(db_1.getIndex(), 1)
        del ds

#------------------------------------------------------------------------------

    def test_int_buffer_from_view(self):
        ds = sidre.DataStore()
        root = ds.getRoot()

        dv = root.createViewAndAllocate("u0", sidre.INT_ID, 10)

        self.assertEqual(dv.getTypeID(), sidre.INT_ID)
#        int * data_ptr = dv.getData()

#        foreach i range(10):
#            data_ptr[i] = i*i

#        dv.print_()

#        self.assertEqual(dv.getTotalBytes(), sizeof(int) * 10)
        del ds

#------------------------------------------------------------------------------

    def XXtest_int_buffer_from_view_conduit_value(self):
        ds = sidre.DataStore()
        root = ds.getRoot()

        dv = root.createViewAndAllocate("u0", sidre.INT_ID, 10)
#        int * data_ptr = dv.getData()

#        for(int i=0  i<10  i++):
#            data_ptr[i] = i*i

        dv.print_()

        self.assertEqual(dv.getTotalBytes(), sizeof(int) * 10)
        del ds

#------------------------------------------------------------------------------

    def XXtest_int_array_strided_views(self):
        ds = sidre.DataStore()
        root = ds.getRoot()
        dbuff = ds.createBuffer()

        dbuff.declare(sidre.INT_ID, 10)
        dbuff.allocate()
#        int * data_ptr = static_cast<int *>(dbuff.getData())

#        for(int i=0  i<10  i++):
#            data_ptr[i] = i

        dbuff.print_()

        self.assertEqual(dbuff.getTotalBytes(), sizeof(int) * 10)

        dv_e = root.createView("even",dbuff)
        dv_o = root.createView("odd",dbuff)
        self.assertEqual(dbuff.getNumViews(), 2)

        # c_int(num_elems, offset [in bytes], stride [in bytes])
        dv_e.apply(sidre.INT_ID, 5, 0, 8)

        # c_int(num_elems, offset [in bytes], stride [in bytes])
        dv_o.apply(sidre.INT_ID, 5, 4, 8)

        dv_e.print_()
        dv_o.print_()

        # Check base pointer case:
#        int* v_e_ptr = dv_e.getData()
#        int* v_o_ptr = dv_o.getData()
#        for(int i=0  i<10  i += 2):
#            std::cout << "idx:" <<  i
#            << " e:" << v_e_ptr[i]
#            << " o:" << v_o_ptr[i]
#            << " em:" << v_e_ptr[i]  % 2
#            << " om:" << v_o_ptr[i]  % 2
#            << std::endl
#
#            self.assertEqual(v_e_ptr[i] % 2, 0)
#            self.assertEqual(v_o_ptr[i] % 2, 1)

        # Check Conduit mem-map struct case:
#        int_array dv_e_ptr = dv_e.getData()
#        int_array dv_o_ptr = dv_o.getData()
#        for(int i=0  i<5  ++i):
#            std::cout << "idx:" <<  i
#            << " e:" << dv_e_ptr[i]
#            << " o:" << dv_o_ptr[i]
#            << " em:" << dv_e_ptr[i]  % 2
#            << " om:" << dv_o_ptr[i]  % 2
#            << std::endl
#
#            self.assertEqual(dv_e_ptr[i] % 2, 0)
#            self.assertEqual(dv_o_ptr[i] % 2, 1)

        # Run similar test to above with different view apply method
        dv_e1 = root.createView("even1",dbuff)
        dv_o1 = root.createView("odd1",dbuff)
        self.assertEqual(dbuff.getNumViews(), 4)

        # (num_elems, offset [in # elems], stride [in # elems])
        dv_e1.apply(sidre.INT_ID, 5,0,2)

        # (num_elems, offset [in # elems], stride [in # elems])
        dv_o1.apply(sidre.INT_ID, 5,1,2)

        dv_e1.print_()
        dv_o1.print_()

        # Check base pointer case:
#        int* v_e1_ptr = dv_e1.getData()
#        int* v_o1_ptr = dv_o1.getData()
#        for(int i=0  i<10  i += 2):
#            std::cout << "idx:" <<  i
#            << " e1:" << v_e1_ptr[i]
#            << " oj:" << v_o1_ptr[i]
#            << " em1:" << v_e1_ptr[i]  % 2
#            << " om1:" << v_o1_ptr[i]  % 2
#            << std::endl
#
#            self.assertEqual(v_e1_ptr[i], v_e_ptr[i])
#            self.assertEqual(v_o1_ptr[i], v_o_ptr[i])

        # Check Conduit mem-map struct case:
#        int_array dv_e1_ptr = dv_e1.getData()
#        int_array dv_o1_ptr = dv_o1.getData()
#        for(int i=0  i<5  i++):
#            std::cout << "idx:" <<  i
#            << " e1:" << dv_e1_ptr[i]
#            << " o1:" << dv_o1_ptr[i]
#            << " em1:" << dv_e1_ptr[i]  % 2
#            << " om1:" << dv_o1_ptr[i]  % 2
#            << std::endl
#
#            self.assertEqual(dv_e1_ptr[i], dv_e_ptr[i])
#            self.assertEqual(dv_o1_ptr[i], dv_o_ptr[i])

        ds.print_()
        del ds

#------------------------------------------------------------------------------

    def XXtest_int_array_depth_view(self):
        ds = sidre.DataStore()
        root = ds.getRoot()
        dbuff = ds.createBuffer()

        depth_nelems = 10 

        # Allocate buffer to hold data for 4 "depth" views
        dbuff.declare(sidre.INT_ID, 4 * depth_nelems )
        dbuff.allocate()
#        int * data_ptr = static_cast<int *>(dbuff.getData())

#        for(size_t i = 0  i < 4 * depth_nelems  ++i):
#            data_ptr[i] = i / depth_nelems

        dbuff.print_()

        self.assertEqual(dbuff.getNumElements(), 4 * depth_nelems)

        # create 4 "depth" views and apply offsets into buffer
        views      = [ None, None, None, None ]
        view_names = [ "depth_0", "depth_1", "depth_2", "depth_3" ]

#        for (int id = 0 id < 2 ++id):
#            views[id] = root.createView(view_names[id], dbuff).apply(depth_nelems, id*depth_nelems)

        # call path including type
#        for (int id = 2 id < 4 ++id):
#            views[id] = root.createView(view_names[id], dbuff).apply(sidre.INT_ID, 
#                                                                     depth_nelems, 
#                                                                     id*depth_nelems)
        self.assertEqual(dbuff.getNumViews(), 4)

        # print depth views...
        for view in views:
            view.print_()

        # check values in depth views...
#        for view in views:
#            int* dv_ptr = view.getData()
#            for (size_t i = 0 i < depth_nelems ++i):
#                self.assertEqual(dv_ptr[i], id)

        ds.print_()
        del ds

#------------------------------------------------------------------------------

    def XXtest_int_array_view_attach_buffer(self):
        ds = sidre.DataStore()
        root = ds.getRoot()

        field_nelems = 10

        # create 2 "field" views with type and # elems
        elem_count = 0 
        field0 = root.createView("field0", sidre.INT_ID, field_nelems)
        elem_count += field0.getNumElements()
        field1 = root.createView("field1", sidre.INT_ID, field_nelems)
        elem_count += field1.getNumElements()
        self.assertEqual(elem_count, 2 * field_nelems)

        # create buffer to hold data for fields and allocate
        dbuff = ds.createBuffer().allocate(sidre.INT_ID, elem_count)
        self.assertEqual(dbuff.getNumElements(), elem_count)

        # Initilize buffer data for testing below.
#        int* b_ptr = dbuff.getData()
#        for(size_t i = 0  i < elem_count  ++i):
#            b_ptr[i] = i / field_nelems

        dbuff.print_()

        # attach field views to buffer and apply offsets into buffer
        field0.attachBuffer(dbuff).apply(field_nelems, 0 * field_nelems)
        field1.attachBuffer(dbuff).apply(field_nelems, 1 * field_nelems)
        self.assertEqual(dbuff.getNumViews(), 2)

        # print field views...
        field0.print_()
        field1.print_()

        # check values in field views...
#        int* f0_ptr = field0.getData()
#        for (size_t i = 0 i < field_nelems ++i):
#            self.assertEqual(f0_ptr[i], 0)
#        int* f1_ptr = field1.getData()
#        for (size_t i = 0 i < field_nelems ++i):
#            self.assertEqual(f1_ptr[i], 1)

        ds.print_()
        del ds

#------------------------------------------------------------------------------

    def XXtest_int_array_multi_view_resize(self):
        #
        # This example creates a 4 * 10 buffer of ints,
        # and 4 views that point the 4 sections of 10 ints
        #
        # We then create a new buffer to support 4*12 ints
        # and 4 views that point into them
        #
        # after this we use the old buffers to copy the values
        # into the new views
        #
        ds = sidre.DataStore()
        root = ds.getRoot()


        # create a group to hold the "old" or data we want to copy
        r_old = root.createGroup("r_old")
        # create a view to hold the base buffer and allocate
        base_old = r_old.createViewAndAllocate("base_data",  sidre.INT_ID, 40)

        # we will create 4 sub views of this array
#        int * data_ptr = base_old.getData()


        # init the buff with values that align with the
        # 4 subsections.
#        for(int i=0  i<10  i++):
#            data_ptr[i] = 1
#        for(int i=10  i<20  i++):
#            data_ptr[i] = 2
#        for(int i=20  i<30  i++):
#            data_ptr[i] = 3
#        for(int i=30  i<40  i++):
#            data_ptr[i] = 4

        # setup our 4 views
        buff_old = base_old.getBuffer()
        buff_old.print_()
        r0_old = r_old.createView("r0",buff_old)
        r1_old = r_old.createView("r1",buff_old)
        r2_old = r_old.createView("r2",buff_old)
        r3_old = r_old.createView("r3",buff_old)

        # each view is offset by 10 * the # of bytes in a int
        # c_int(num_elems, offset)
        offset =0
        r0_old.apply(sidre.INT_ID, 10, offset)

        offset += sizeof(int) * 10
        r1_old.apply(sidre.INT_ID, 10, offset)

        offset += sizeof(int) * 10
        r2_old.apply(sidre.INT_ID, 10, offset)

        offset += sizeof(int) * 10
        r3_old.apply(sidre.INT_ID, 10,offset)

        # check that our views actually point to the expected data
        #
#        int * r0_ptr = r0_old.getData()
#        for(int i=0  i<10  i++):
#            self.assertEqual(r0_ptr[i], 1)
#            # check pointer relation
#            self.assertEqual(&r0_ptr[i], &data_ptr[i])

#        int * r3_ptr = r3_old.getData()
#        for(int i=0  i<10  i++):
#            self.assertEqual(r3_ptr[i], 4)
#            # check pointer relation
#            self.assertEqual(&r3_ptr[i], &data_ptr[i+30])

        # create a group to hold the "old" or data we want to copy into
        r_new = root.createGroup("r_new")
        # create a view to hold the base buffer and allocate
        base_new = r_new.createViewAndAllocate("base_data", sidre.INT_ID, 4 * 12)

#        int * base_new_data = base_new.getData()
#        for (int i = 0  i < 4 * 12  ++i):
#            base_new_data[i] = 0

        buff_new = base_new.getBuffer()
        buff_new.print_()

        # create the 4 sub views of this array
        r0_new = r_new.createView("r0",buff_new)
        r1_new = r_new.createView("r1",buff_new)
        r2_new = r_new.createView("r2",buff_new)
        r3_new = r_new.createView("r3",buff_new)

        # apply views to r0,r1,r2,r3
        # each view is offset by 12 * the # of bytes in a int

        # c_int(num_elems, offset)
        offset =0
        r0_new.apply(sidre.INT_ID, 12, offset)

        offset += sizeof(int) * 12
        r1_new.apply(sidre.INT_ID, 12, offset)

        offset += sizeof(int) * 12
        r2_new.apply(sidre.INT_ID, 12, offset)

        offset += sizeof(int) * 12
        r3_new.apply(sidre.INT_ID, 12, offset)

        # update r2 as an example first
        buff_new.print_()
        r2_new.print_()

        # copy the subset of value
        r2_new.getNode().update(r2_old.getNode())
        r2_new.getNode().print_()
        buff_new.print_()


        # check pointer values
#        int * r2_new_ptr = r2_new.getData()

#        for(int i=0  i<10  i++):
#            self.assertEqual(r2_new_ptr[i], 3)

#        for(int i=10  i<12  i++):
#            self.assertEqual(r2_new_ptr[i], 0)     # assumes zero-ed alloc

        # update the other views
        r0_new.getNode().update(r0_old.getNode())
        r1_new.getNode().update(r1_old.getNode())
        r3_new.getNode().update(r3_old.getNode())

        buff_new.print_()


        ds.print_()
        del ds

#------------------------------------------------------------------------------

    def XXtest_int_array_realloc(self):
        ds = sidre.DataStore()
        root = ds.getRoot()
    
        # create a view to hold the base buffer
        a1 = root.createViewAndAllocate("a1", sidre.FLOAT_ID, 5)
        a2 = root.createViewAndAllocate("a2", sidre.INT_ID, 5)

#        float * a1_ptr = a1.getData()
#        int * a2_ptr = a2.getData()

#        for(int i=0  i<5  i++):
#            a1_ptr[i] =  5.0
#            a2_ptr[i] = -5

        self.assertEqual(a1.getTotalBytes(), sizeof(float)*5)
        self.assertEqual(a2.getTotalBytes(), sizeof(int)*5)


        a1.reallocate(sidre.FLOAT_ID, 10)
        a2.reallocate(sidre.INT_ID, 15)

        a1_ptr = a1.getData()
        a2_ptr = a2.getData()

#        for(int i=0  i<5  i++):
#            self.assertEqual(a1_ptr[i],5.0)
#            self.assertEqual(a2_ptr[i],-5)

#        for(int i=5  i<10  i++):
#            a1_ptr[i] = 10.0
#            a2_ptr[i] = -10

#        for(int i=10  i<15  i++):
#            a2_ptr[i] = -15

        self.assertEqual(a1.getTotalBytes(), sizeof(float)*10)
        self.assertEqual(a2.getTotalBytes(), sizeof(int)*15)

        # Try some errors
        # XXX  a1.reallocate(sidre.INT_ID(20))
        # XXX reallocate with a Schema
        
        ds.print_()
        del ds

#------------------------------------------------------------------------------

    def XXtest_simple_opaque(self):
        ds = sidre.DataStore() # create our main data store
        root = ds.getRoot()    # get access to our root data Group

#        int * src_data = new int[1]

        src_data[0] = 42

#        void * src_ptr = (void *)src_data

        opq_view = root.createView("my_opaque", src_ptr)

        # we have a buffer because an "external" view currently uses one
        self.assertEqual(ds.getNumBuffers(), 1)

        self.assertTrue(opq_view.isExternal())
        self.assertTrue(not opq_view.isApplied())
        self.assertTrue(opq_view.isOpaque())

#        void * opq_ptr = opq_view.getVoidPtr()

#        int * out_data = (int *)opq_ptr
        self.assertEqual(opq_ptr,src_ptr)
        self.assertEqual(out_data[0],42)

        ds.print_()
        del ds
#        delete [] src_data


if __name__ == '__main__':
    unittest.main()
