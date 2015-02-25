from collections import OrderedDict
import os.path

import numpy as np

class Node(object):
    """Contains attributes and connectivity information.
    """
    def __init__(self, name, parent, contents):
        self.name = name
        self.parent = parent
        self.attrs = {}
        self.contents = contents   # Data, Directory, ...

    def __getattr__(self, name):
        try:
            return getattr(self.contents, name)
        except AttributeError:
            raise AttributeError(name)

    def __getitem__(self, key):
        """Ask contents to index. """
        return self.contents[key]

    def set_attr(self, key, value):
        """ set attribute."""
        self.attrs[key] = value
        return self

    def get_attr(self, key):
        return self.attrs[key]

    def has_attr(self, key):
        return key in self.attrs

    def del_attr(self, key):
        del self.attrs[key]


    def get_parent(self):
        return self.parent

    def print_tree(self, depth=0):
        print "  "*depth, self.name

    def collect_attr_work(self, aname, rv):
        """collect all nodes which have attribute aname.
        """
        if aname in self.attrs:
            rv.add_item(self)


class Data(object):
    def __init__(self, shape=None, dtype=None, array=None):
        self.shape = shape
        self.dtype = dtype
        self.array = array
        self.contents_type = self.__class__

    def __getitem__(self, key):
        if self.array is not None:
            return self.array[key]
        raise IndexError

    def set_descr(self, shape, dtype):
        if self.array is not None:
            raise RuntimeError("Unable to change allocated array")
        self.shape = shape
        self.dtype = dtype

    def set_array(self, array):
        """Give an existing array to the object.
        """
        try:
            self.shape = array.shape
            self.dtype = array.dtype.type
            self.array = array
        except:
            raise

    def get_shape(self):
        return self.shape

    def get_type(self):
        return self.dtype

    def get_array(self):
        return self.array

    def allocate(self):
        if not self.shape or not self.dtype:
            raise RuntimeError("Must call set_descr first")
        self.array = np.zeros(self.shape, self.dtype)
        return self.array


class Directory(object):
    def __init__(self):
        self.items = OrderedDict()
        self.contents_type = self.__class__

    def __iter__(self):
        return iter(self.items.values())

    def add_item(self, item):
        """
        item must be a Node
        """
        self.items[item.name] = item

    def insert_item(self, item, fullname):
        """ Create a node wrapper around item and insert into the directory.
        """
        path, basename = os.path.split(fullname)
        node = Node(basename, self.node, item)
        item.node = node
        current = self
        if path:
            for part in fullname.split('/'):
                next = current.items.get(part)
                if next is None:
                    next = current.new_directory(part)
                current = next
        current.add_item(node)
        return node

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.items[  self.items.keys()[key] ]
        else:
            path = key.split('/')
            name = path.pop()
            current = self
            for node in path:
                current = current.items.get(node)
            return current.items[name]

    # helper functions
    def new_data(self, name):
        """ Create a Data node and insert into tree.
        """
        item = Data()
        return self.insert_item(item, name)

    def new_directory(self, name):
        item = Directory()
        return self.insert_item(item, name)

    def print_tree(self, depth=0):
        print "  "*depth, self.name
        depth += 1
        for item in self.items.values():
            item.print_tree(depth)

    def collect_attr_work(self, aname, rv):
        """collect all nodes which have attribute aname.
        May be called on any directory.
        """
        if hasattr(self, aname):
            rv.add_item(self)
        for item in self.items.values():
            item.collect_attr_work(aname, rv)

    def collect_attr(self, aname):
        """collect all nodes which have attribute aname.
        May be called on any directory.
        """
        rv = NodeItems(aname)
        self.collect_attr_work(aname, rv)
        return rv


class NodeItems(object):
    """  Not accessiable as a dictionary, only as a list
    since duplicate names may be in the list - unless indexed by full path name
    """
    def __init__(self, name):
        self.name = name
        self.items = []

    def __iter__(self):
        return iter(self.items)

    def add_item(self, item):
        self.items.append(item)



the_datastore = None


def init_library():

    # create the root object
    item = Directory()
    node = Node("/", None, item)
    item.node = node
    global the_datastore
    the_datastore = node
    return node

def fin_library():
    pass
