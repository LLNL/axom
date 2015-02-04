#
# testfunc.py - test functions
#

import datastore as ds

def worker1():
    print("In worker1")

worker1()

myDS1 = ds.create_data_store("myDS1")
print("myDS1 =", myDS1)

# Add to datastore
obj1 = ds.add_function(myDS1, "worker1", worker1)

# call object directly
print("Call object directly:")
obj1()

print("Call object via C:")
dstest.call_python_function(obj1)

