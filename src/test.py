#
# Sample use of datastore API without concern for the C++ implementation.
#

import datastore as ds
import numpy as np


# initializing the library creates a root directory
root = ds.init_library()
print "Global datastore", ds.the_datastore

print root
print root.name

data1 = root.new_data("name1").set_attr("dump", True)

print data1.name
print data1.contents_type
# print data1['foo']    data1 is not indexable by name

root.new_data("name2")
data2 = root["name2"]
print data2.name

# add existing item
data3 = ds.Data()
root.insert_item(data3, "data3")

dir1 = root.new_directory("dir1")
dir1.set_attr("dump", True)

print dir1.name
print dir1.contents_type
print "has dump1 =",  dir1.has_attr("dump")
print "dir1.dump =", dir1.get_attr("dump")
print "has nosuch =",  dir1.has_attr("nosuch")

try:
    print "dir1.nosuch =", dir1.get_attr("nosuch")
except:
    print "Exception"


dir1 = root.new_directory("dir1")
var = dir1.new_data("dname1").set_attr("dump", True)

# implicitly creates dir2
root.new_data("dir1/dir2/dname12")

root['dir1/dir2'].set_attr("dump", True)
# del_attr('dump')


# Add another Data to root after Directory node
root.new_data("lastname").set_attr("dump", True)


print "TREE"
root.print_tree()

print "Parents"
print "data1", data1.get_parent().name
print "dir1", root['dir1'].get_parent().name
print "dir1/dir2", root['dir1/dir2'].get_parent().name
print "dir1/dir2/dname12", root['dir1/dir2/dname12'].get_parent().name
print
print "iterate over root"
for item in root:
    print item.name, item.__class__
    

print
print "Dquery-1 - Access object by name"
print "Dquery-2 - Access object by index"
for i in [ 'name2', 'dir1/dir2/dname12', 1]:
    print 'root[', i, '] = ', root[i].name

print "Dquery-3 - Access set of objects by attribute"

print "items with dump"
newset = root.collect_attr('dump')
for item in newset:
    print item.name, item.__class__

# Dquery-4 - Access set of objects by type


# Dquery-5 - Return the number of objects that have a particular attribute.

# Dquery-6 - Answer whether an object with a given name has a particular attribute.
print
print "data1.has_attr('dump')", data1.has_attr('dump')
print "data1.has_attr('nosuch')", data1.has_attr('nosuch')

# Dquery-7 - Support boolean expression queries that combine multiple attributes,
#            attribute and type information, etc.




####### actual data
# using numpy's types for now  np.int32

print "before data1.descr = ", data1.get_shape(), data1.get_type(), data1.get_array()
data1.set_descr(shape=20, dtype=np.int32)
print "after data1.descr = ", data1.get_shape(), data1.get_type(), data1.get_array()
data1.allocate()
print "alloc data1.descr = ", data1.get_shape(), data1.get_type(), data1.get_array()

# change allocated array
#data1.set_descr(shape=10, dtype=np.int32)

b = np.arange(10)
data2.set_array(b)
print "after data2.descr = ", data2.get_shape(), data2.get_type(), data2.get_array()

# Index a Data object directly as an array
print "data2[3] = ", data2[3]
# data3[3]    exception since data3 has no array
