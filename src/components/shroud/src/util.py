#!/bin/env python
from __future__ import print_function


import collections
import copy
import string
import json
import os

import parse_decl

fmt = string.Formatter()

default_template = dict(
#    C_name='{C_prefix}{lower_class}_{underscore_name}{function_suffix}',

#    C_header_filename = 'wrap{cpp_class}.h',
#    C_impl_filename = 'wrap{cpp_class}.cpp',

#    F_name_impl = '{lower_class}_{underscore_name}{function_suffix}',
#    F_name_method = '{underscore_name}{function_suffix}',
#    F_name_generic = '{underscore_name}',
)


def wformat(template, dct):
    # shorthand, wrap fmt.vformat
    return fmt.vformat(template, None, dct)

def append_format(lst, template, dct):
    # shorthand, wrap fmt.vformat
    lst.append(fmt.vformat(template, None, dct))

def eval_template(node, name, tname='', fmt=None):
    """fmt[name] = node[name] or option[name + tname + '_template']
    """
    if fmt is None:
        fmt = node['fmt']
    if name in node:
        setattr(fmt, name, node[name])
    else:
        tname = name + tname + '_template'
        setattr(fmt, name, wformat(node['options'][tname], fmt))


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


class WrapperMixin(object):
    """Methods common to all wrapping classes.
    """

#####

    def _init_splicer(self, splicers):
        self.splicers = splicers
        self.splicer_stack = [ splicers ]
        self.splicer_names = [ ]
        self.splicer_path = ''

    def _push_splicer(self, name):
        level = self.splicer_stack[-1].setdefault(name, {})
        self.splicer_stack.append(level)
        self.splicer_names.append(name)
        self.splicer_path = '.'.join(self.splicer_names) + '.'

    def _pop_splicer(self, name):
        # XXX maybe use name for error checking, must pop in reverse order
        self.splicer_stack.pop()
        self.splicer_names.pop()
        if self.splicer_names:
            self.splicer_path = '.'.join(self.splicer_names) + '.'
        else:
            self.splicer_path = ''

    def _create_splicer(self, name, out, override=None, default=[]):
        # The prefix is needed when two different sets of output are being create
        # and they are not in sync.
        # Creating methods and derived types together.
        out.append('%s splicer begin %s%s' % (self.comment, self.splicer_path, name))
        if override:
            out.extend(override)
        else:
            out.extend(self.splicer_stack[-1].get(name, default))
        out.append('%s splicer end %s%s' % (self.comment, self.splicer_path, name))

#####

    def std_c_type(self, lang, arg, const=None, ptr=False):
        """
        Return the C type.
        pass-by-value default

        lang = c_type or cpp_type
        if const is None, use const from arg.

        attributes:
        ptr - True = pass-by-reference
        reference - True = pass-by-reference

        """
#        if lang not in [ 'c_type', 'cpp_type' ]:
#            raise RuntimeError
        t = []
        typedef = self.typedef.get(arg['type'], None)
        if typedef is None:
            raise RuntimeError("No such type %s" % arg['type'])

        if const is None:
            const = arg['attrs'].get('const', False)
        if const:
            t.append('const')

        t.append(getattr(typedef, lang))
        if arg['attrs'].get('ptr', ptr):
            t.append('*')
        elif arg['attrs'].get('reference', False):
            if lang == 'cpp_type':
                t.append('&')
            else:
                t.append('*')
        return ' '.join(t)

    def std_c_decl(self, lang, arg,
                   name=None, const=None, ptr=False):
        """
        Return the C declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ = self.std_c_type(lang, arg, const, ptr)
        return typ + ' ' + ( name or arg['name'] )

#####

    def namespace(self, node, position, output):
        options = node['options']
        namespace = options.namespace
        if position == 'begin':
            for name in namespace.split():
                output.append('namespace %s {' % name)
        else:
            for name in namespace.split():
                output.append('}  // namespace %s' % name)
#####

    def write_output_file(self, fname, directory, output):
        fp = open(os.path.join(directory, fname), 'w')
        fp.write('%s %s\n' % (self.comment, fname))
        fp.write(self.comment + ' This is generated code, do not edit\n')
        self.write_copyright(fp)
        self.indent = 0
        self.write_lines(fp, output)
        fp.close()
        self.log.write("Close %s\n" % fname)
        print("Wrote", fname)

    def write_copyright(self, fp):
        for line in self.tree.get('copyright', []):
            if line:
                fp.write(self.comment + ' ' + line + '\n')
            else:
                fp.write(self.comment + '\n')

    def write_lines(self, fp, lines):
        """ Write lines with indention and newlines.
        """
        for line in lines:
            if isinstance(line, int):
                self.indent += int(line)
            else:
                for subline in line.split("\n"):
                    fp.write('    ' * self.indent)
                    fp.write(subline)
                    fp.write('\n')

    def write_doxygen_file(self, output, fname, node, cls):
        """ Write a doxygen comment block for a file.
        """
        output.append(self.doxygen_begin)
        output.append(self.doxygen_cont + ' \\file %s' % fname)
        if cls:
            output.append(self.doxygen_cont +
                          ' \\brief Shroud generated wrapper for %s class'
                          % node['name'])
        else:
            output.append(self.doxygen_cont +
                          ' \\brief Shroud generated wrapper for %s library'
                          % node['options'].library)
        output.append(self.doxygen_end)

    def write_doxygen(self, output, docs):
        """Write a doxygen comment block for a function.
        Uses brief, description, and return from docs.
        """
        output.append(self.doxygen_begin)
        if 'brief' in docs:
            output.append(self.doxygen_cont + ' \\brief %s' % docs['brief'])
            output.append(self.doxygen_cont)
        if 'description' in docs:
            desc = docs['description']
            if desc.endswith('\n'):
                lines = docs['description'].split('\n')
                lines.pop()  # remove trailing newline
            else:
                lines = [desc]
            for line in lines:
                output.append(self.doxygen_cont + ' ' + line)
        if 'return' in docs:
            output.append(self.doxygen_cont)
            output.append(self.doxygen_cont + ' \\return %s' % docs['return'])
        output.append(self.doxygen_end)


class Typedef(object):
    """ Collect fields for an argument.
    This used to be a dict but a class has better access semantics: i.attr vs d['attr']
    It also initializes default values to avoid  d.get('attr', default)
    """
    # valid fields
    defaults = dict(
        base='unknown',       # Base type: 'string'
        forward=None,         # Forward declaration
        typedef=None,         # Initialize from existing type

        cpp_type=None,        # Name of type in C++
        cpp_to_c='{var}',     # Expression to convert from C++ to C
        cpp_header=None,      # Name of C++ header file required for implementation
                              # For example, if cpp_to_c was a function

        c_type=None,          # Name of type in C
        c_header=None,        # Name of C header file required for type
        c_to_cpp='{var}',     # Expression to convert from C to C++
        c_fortran=None,       # Expression to convert from C to Fortran
        c_argdecl=None,       # List of argument declarations for C wrapper, None=match declaration
                              # used with string_from_buffer 
        c_pre_call = None,    # Statement to execute before call
        c_post_call = None,   # Statement to execute after call
        c_return_code=None,

        f_c_args=None,        # List of argument names to F_C routine
        f_c_argdecl=None,     # List of declarations to F_C routine

        f_type=None,          # Name of type in Fortran
        f_derived_type=None,  # Fortran derived type name
        f_args=None,          # Argument in Fortran wrapper to call C.
        f_module=None,        # Fortran modules needed for type  (dictionary)
        f_return_code=None,
        f_kind = None,        # Fortran kind of type
        f_cast = '{var}',     # Expression to convert to type
                              # e.g. intrinsics such as int and real
        f_use_tmp = False,    # Pass {tmp_var} to C routine instead of {var}
        f_argsdecl = None,    # List of declarations need by argument.
        f_pre_call = None,    # Statement to execute before call, often to coerce types
        f_post_call = None,   # Statement to execute after call - cleanup, coerce result

# XXX - maybe later.  For not in wrapping routines
#        f_attr_len_trim = None,
#        f_attr_len = None,
#        f_attr_size = None,

        result_as_arg = None, # override fields when result should be treated as an argument

        PY_format='O',        # 'format unit' for PyArg_Parse
        PY_PyTypeObject=None, # variable name of PyTypeObject instance
        PY_PyObject=None,     # typedef name of PyObject instance
        PY_ctor=None,         # expression to create object.
                              # ex. PyBool_FromLong({rv})
        PY_to_object=None,    # PyBuild - object = converter(address)
        PY_from_object=None,  # PyArg_Parse - status = converter(object, address);
        PY_post_parse='KKK',  # Used if PY_PyTypeObject is set
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

    def XXXcopy(self):
        n = Typedef(self.name)
        n.update(self._to_dict())
        return n

    def clone_as(self, name):
        n = Typedef(name)
        n.update(self._to_dict())
        return n

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

    def __contains__(self, item):
        return hasattr(self, item)

    def __repr__(self):
        return str(self._to_dict())

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

    def inlocal(self, key):
        """ Return true if key is defined locally
        i.e. does not check parent.
        """
        return key in self.__dict__

    def _to_dict(self):
        d = {}
        skip = '_' + self.__class__.__name__ + '__'   # __name is skipped
        for key, value in self.__dict__.items():
            if not key.startswith(skip):
                d[key] = value
        return d

    def _to_full_dict(self, d=None):
        if d is None:
            d = self._to_dict()
        else:
            d.update( self._to_dict())
        if self.__parent:
            self.__parent._to_full_dict(d)
        return d


def copy_function_node(node):
    """Create a copy of a function node to use with C++ template.
    """
    # Shallow copy everything
    new = node.copy()

    # Deep copy dictionaries
    for field in [ 'args', 'attrs', 'result' ]:
        new[field] = copy.deepcopy(node[field])

    # Add new Options in chain
    for field in [ 'fmt', 'options' ]:
        new[field] = Options(node[field])

    return new

class XXXClassNode(object):
    """Represent a class.  Usually a C++ class.
    It'd be nice if a group of related function which are used in an o-o manner
    would also be treated as a class.
    """
    def __init__(self, name):
        self.name = name
        self.options = None
        self.methods = []


class XXXFunctionNode(object):
    def __init__(self):
        self.decl = None
        self.result = {}
        self.args = []
        self.function_suffix = ''
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

    print("IN  a", 'a' in lev1)
    print("IN  z", 'z' in lev1)
    print("IN  c2", 'c2' in lev1)
    print("IN  z2", 'z2' in lev1)
    print("IN nosuch", 'nosuch' in lev1)


    print("FORMAT:", wformat("{a} {z} {c2} {z2}", lev1))

    print(json.dumps(lev0, cls=ExpandedEncoder, sort_keys=True))
    print(json.dumps(lev1, cls=ExpandedEncoder, sort_keys=True))
