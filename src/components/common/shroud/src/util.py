#!/bin/env python
from __future__ import print_function


import collections
import string
import json

import parse_decl

fmt = string.Formatter()

default_template = dict(
    C_name_template='{C_prefix}{lower_class}_{underscore_name}{method_suffix}',
)


def wformat(template, d):
    # shorthand, wrap fmt.vformat
    return fmt.vformat(template, None, d)

def eval_templates(templates, node, fmt_dict):
    """ If variable is not already set in node, then
    set based on template.
    Update value in fmt_dict.
    """
    for tname in ( templates ) :
        if tname not in node:
            node[tname] = fmt.vformat(
                getattr(node['options'], tname + '_template'),
                None, fmt_dict)
        fmt_dict[tname] = node[tname]

def eval_template2(options, tname, fmt, vname, default):
    """ If a tname exists in options, use it; else use default.
    fmt[vname] = option[tname]
    """
    setattr(fmt, vname, wformat(options.get(tname, default), fmt))

def eval_template3(options, fmt, name, default):
    """ If a tname exists in options, use it; else use default.
    fmt[vname] = option[tname]
    """
    setattr(fmt, name, wformat(options.get(name + '_template', default), fmt))

def eval_template4(options, fmt, name):
    """ If a tname exists in options, use it; else use default.
    fmt[vname] = option[tname]
    """
    if hasattr(options, name):
        setattr(fmt, name, getattr(options, name))
    else:
        tname = name + '_template'
        setattr(fmt, name, wformat(options.get(tname,
                                               default_template[tname]), fmt))


# http://stackoverflow.com/questions/1175208/elegant-python-function-to-convert-camelcase-to-camel-case
def un_camel(text):
    """ Converts a CamelCase name into an under_score name. 

        >>> un_camel('CamelCase')
        'camel_case'
        >>> un_camel('getHTTPResponseCode')
        'get_http_response_code'
    """
    result = []
    pos = 0
    while pos < len(text):
        if text[pos].isupper():
            if pos-1 > 0 and text[pos-1].islower() or pos-1 > 0 and \
            pos+1 < len(text) and text[pos+1].islower():
                result.append("_%s" % text[pos].lower())
            else:
                result.append(text[pos].lower())
        else:
            result.append(text[pos])
        pos += 1
    return "".join(result)

# http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


class Typedef(object):
    """ Collect fields for an argument.
    This used to be a dict but a class has better access semantics: i.attr vs d['attr']
    It also initializes default values to avoid  d.get('attr', default)
    """
    # valid fields
    defaults = dict(
        base='unknown',       # base type: 'string'
        forward=None,         # forward declaration

        cpp_type=None,        # name of type in C++
        cpp_to_c='{var}',     # expression to convert from C++ to C

        c_type=None,          # name of type in C
        c_header=None,        # Name of C header file required for type
        c_to_cpp='{var}',     # expression to convert from C to C++
        c_fortran=None,       # expression to convert from C to Fortran

        f_type=None,         # name of type in Fortran
        fortran_derived=None,    # Fortran derived type name
        fortran_to_c='{var}', # expression to convert Fortran to C
        f_module=None,        # Fortran modules needed for type  (dictionary)
        f_return_code='{F_result} = {F_C_name}({arg_c_call})',
        )

    def __init__(self, name, **kw):
        self.name = name
#        for key, defvalue in self.defaults.items():
#            setattr(self, key, defvalue)
        self.__dict__.update(self.defaults) # set all default values
        self.update(kw)

    def update(self, d):
        """Add options from dictionary to self.
        """
        for key in d:
            if key in self.defaults:
                setattr(self, key, d[key])
            else:
                raise RuntimeError("Unknown key for Argument %s", key)

    def _to_dict(self):
        """Convert instance to a dictionary for json.
        """
        # only export non-default values
        a = {}
        for key, defvalue in self.defaults.items():
            value = getattr(self, key)
            if value is not defvalue:
                a[key] = value
        return a

    def __repr__(self):
        # only print non-default values
        args = []
        for key, defvalue in self.defaults.items():
            value = getattr(self, key)
            if value is not defvalue:
                if isinstance(value, str):
                    args.append("{0}='{1}'".format(key, value))
                else:
                    args.append("{0}={1}".format(key, value))
        return "Typedef('%s', " % self.name + ','.join(args) + ')'


class Options(object):
    """
    If attribute is not found, look in parent's.
    A replacement for a dictionary to allow obj.name syntax.
    It will automatically look in __parent to attribute if not found to allow
    A nesting of options.
    Use __attr to avoid xporting to json
    """
    def __init__(self, parent, **kw):
        self.__parent = parent
        self.__hidden = 43
        self.update(kw)

    def __getattr__(self, name):
        # we get here if the attribute does not exist in current instance
        if self.__parent:
            return getattr(self.__parent, name)
        else:
            raise AttributeError("%r object has no attribute %r" %
                                 (self.__class__.__name__, name))

    def __getitem__(self, key):
        """ Treat as dictionary for format command.
        """
        return getattr(self, key)

    def get(self, key, value=None):
        """ D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
        """
        try:
            return getattr(self, key)
        except AttributeError:
            return value

    def setdefault(self, key, value=None):
        """ D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D
        """
        if key not in self.__dict__:
            self.__dict__[key] = value
        return self.__dict__.get(key, value)

    def update(self, d, replace=True):
        """Add options from dictionary to self.
        """
        for key, value in d.items():
            if replace:
                setattr(self, key, value)
            elif not hasattr(self, key):
                setattr(self, key, value)

    def _to_dict(self):
        d = {}
        skip = '_' + self.__class__.__name__ + '__'   # __name is skipped
        for key, value in self.__dict__.items():
            if not key.startswith(skip):
                d[key] = value
        return d


class XXXClassNode(object):
    """Represent a class.  Usually a C++ class.
    It'd be nice if a group of related function which are used in an o-o manner
    would also be treated as a class.
    """
    def __init__(self, name):
        self.name = name
        self.options = None
        self.methods = []


class FunctionNode(object):
    def __init__(self):
        self.decl = None
        self.result = {}
        self.args = []
        self.method_suffix = ''
        self.arg_map = {}

    def set_decl(self, decl):
        """decl will compute result and args
        """
        self.decl = decl
        values = parse_decl.check_decl(decl)
        self._update_result_args(values)

    def update(self, d):
        """Update from a dictionary.
        """
        if 'decl' in d:
            self.set_decl(d['decl'])
        self._update_result_args(d)

    def _update_result_args(self, d):
        # Just assign if no value yet.
        if 'result' in d:
            if self.result:
                update(self.result, d['result'])
            else:
                self.result = d['result']
        if 'args' in d:
            if self.args:
                for arg in d['args']:
                    name = arg['name']
                    if name in self.arg_map:
                        # update existing arg
                        update(self.arg_map[name], arg)
                    else:
                        # append to current args
                        self.args.append(arg)
                        self.arg_map[name] = arg
            else:
                self.args = d['args']
                for arg in self.args:
                    name = arg['name']
                    self.arg_map[name] = arg

    def dump(self):
        print('FunctionNode:', self.decl)
        print(self.result)
        print(self.args)
        


class ExpandedEncoder(json.JSONEncoder):
    """Jason handler to convert objects into a dictionary when they have
    a _to_dict method.
    """
    def default(self, obj):
        if hasattr(obj, '_to_dict'):
            return obj._to_dict()
         # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


if __name__ == '__main__':
    print(un_camel('incrementCount'))
    print(un_camel('local_function1'))
    print(un_camel('getHTTPResponseCode'))


    a = dict(name='aa', type='int', attrs=dict(ptr=True))
    b = dict(attrs=dict(len='somefunction'))

    print("AAAA", a)
    print("BBBB", b)
#    b.update(a)
    update(b, a)
    print("CCCC", b)


    # Argument
    print("Test Typedef")
    a = Typedef('top', base='new_base') #, bird='abcd')
    print(json.dumps(a, cls=ExpandedEncoder, sort_keys=True))

    # Option
    print("Test Option")
    lev0 = Options(None, a=1, b=2, c=3)
    lev1 = Options(lev0, x=100, y=01, z=102)
    lev0.c2 = 32
    lev1.z2 = 103

    print(lev0.a)
    print(lev0.c2)
    #print(lev0.z)
    print(lev1.a)
    print(lev1.z)
    print(lev1.c2)
    print(lev1.z2)
    try:
        print(lev1.nosuch)
    except AttributeError:
        print("Passed nosuch attribute")

    lev1.setdefault('yyy', 'yyyvalue')

    print("GET",  lev1.get('a', 'notfound'))
    print("GET",  lev1.get('nosuch', 'notfound'))


    print("FORMAT:", wformat("{a} {z} {c2} {z2}", lev1))

    print(json.dumps(lev0, cls=ExpandedEncoder, sort_keys=True))
    print(json.dumps(lev1, cls=ExpandedEncoder, sort_keys=True))
