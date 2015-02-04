#
# test1.py - test datastore.
#

import numpy as np
import datastore as ds

myDS1 = ds.create_data_store("myDS1")
print("myDS1 =", myDS1)

#a = np.array([1,2,3,4])
a = 10

obj1 = ds.add_object(myDS1, "var1", a)
print("obj1", obj1)

obj1a = myDS1.var1
print("myDS1.var1", obj1a)

if obj1 is obj1a:
    print("Same object")
else:
    print("Different objects", id(obj1), id(obj1a))
