#!/bin/env python
#
#  Examples of PyBindgen
#

from pybindgen import *
import sys

mod = Module('sidre')
#mod.add_include('"my-module.h"')

# Declare all to allow for forward references
datastore  = mod.add_class('DataStore')
datagroup  = mod.add_class("DataGroup")
dataview   = mod.add_class("DataView")
databuffer = mod.add_class("DataBuffer")


datastore.add_constructor([])
datastore.add_method('SetInt', None, [param('int', 'value')])
datastore.add_method('getRoot', None, [])
datastore.add_method('getRoot2', retval('DataGroup *', caller_owns_return=True), [])
datastore.add_method('getRoot3', retval('DataGroup *', caller_owns_return=False), [])

#datastore.add_method('GetInt', retval('int'), [], is_const=True)

datagroup.add_method('createView',
                     retval('DataView *', caller_owns_return=False),
                     [param('char *', 'name'),
                      param('DataBuffer *', 'buf', transfer_ownership=False)])



mod.generate(sys.stdout)
