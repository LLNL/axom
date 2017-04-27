#!/bin/env python
from __future__ import print_function
from __future__ import absolute_import


import collections
import copy
import string
import json
import os

from . import parse_decl

fmt = string.Formatter()

default_template = dict(
    # C_name='{C_prefix}{class_lower}_{underscore_name}{function_suffix}',

    # C_header_filename = 'wrap{cpp_class}.h',
    # C_impl_filename = 'wrap{cpp_class}.cpp',

    # F_name_impl = '{class_lower}_{underscore_name}{function_suffix}',
    # F_name_function = '{underscore_name}{function_suffix}',
    # F_name_generic = '{underscore_name}',
)


def wformat(template, dct):
    # shorthand, wrap fmt.vformat
    try:
        return fmt.vformat(template, None, dct)
    except AttributeError as e:
        print(e)
        raise SystemExit('Error with template: ' + template)


def append_format(lst, template, dct):
    # shorthand, wrap fmt.vformat
    lst.append(wformat(template, dct))


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


def as_yaml(obj, order, indent, output):
    """Write out obj in YAML syntax
    obj    - a dictionary or an instance with attributes to dump.
    order  - order of keys to dump
    indent - indention level.
    output - list of output lines.

    This is not really intendent to be a general routine.
    It has some knowledge of what it expects in order to create
    a YAML file similar to what a user may write.
    """

    prefix = "  " * indent
    for key in order:
        if isinstance(obj, collections.Mapping):
            value = obj[key]
        else:
            value = getattr(obj, key)

        if not value:
            # skip empty values such as None or {}
            pass
        elif isinstance(value, basestring):
            # avoid treating strings as a sequence
            # quote strings which start with { to avoid treating them
            # as a dictionary.
            if value.startswith('{'):
                output.append('{}{}: "{}"'.format(prefix, key, value))
            else:
                output.append('{}{}: {}'.format(prefix, key, value))
        elif isinstance(value, collections.Sequence):
            # Keys which are are an array of string (code templates)
            if key in ('declare', 'pre_call', 'pre_call_trim', 'post_call',
                       'post_parse', 'ctor',
                   ):
                output.append('{}{}: |'.format(prefix, key))
                for i in value:
                    output.append('{}  {}'.format(prefix, i))
            else:
                output.append('{}{}:'.format(prefix, key))
                for i in value:
                    output.append('{}- {}'.format(prefix, i))
        elif isinstance(value, collections.Mapping):
            output.append('{}{}:'.format(prefix, key))
            order0 = value.keys()
            order0.sort()
            as_yaml(value, order0, indent + 1, output)
        else:
            # numbers or booleans
            output.append('{}{}: {}'.format(prefix, key, value))


def typedef_wrapped_defaults(typedef):
    """Add some defaults to typedef.
    When dumping typedefs to a file, only a subset is written
    since the rest are boilerplate.  This function restores
    the boilerplate.
    """
    if typedef.base != 'wrapped':
        return

    typedef.cpp_to_c=('static_cast<{c_const}%s *>('
                      'static_cast<{c_const}void *>({cpp_var}))' %
                      typedef.c_type)

    # opaque pointer -> void pointer -> class instance pointer
    typedef.c_to_cpp=('static_cast<{c_const}%s{c_ptr}>('
                      'static_cast<{c_const}void *>({c_var}))' %
                      typedef.cpp_type)

    typedef.f_type='type(%s)' % typedef.f_derived_type
    typedef.f_c_type='type(C_PTR)'

    # XXX module name may not conflict with type name
#    typedef.f_module={fmt_class.F_module_name:[unname]}

    # return from C function
    # f_c_return_decl='type(CPTR)' % unname,
    typedef.f_return_code=('{F_result}%{F_derived_member} = '
                           '{F_C_call}({F_arg_c_call_tab})')
    typedef.f_c_module={ 'iso_c_binding': ['C_PTR']}

    typedef.py_statements=dict(
        intent_in=dict(
            post_parse=[
                '{cpp_var} = {py_var} ? {py_var}->{BBB} : NULL;',
            ],
        ),
        intent_out=dict(
            ctor=[
                ('{PyObject} * {py_var} = '
                 'PyObject_New({PyObject}, &{PyTypeObject});'),
                '{py_var}->{BBB} = {cpp_var};',
            ]
        ),
    )
    # typedef.PY_ctor='PyObject_New({PyObject}, &{PyTypeObject})'

    typedef.LUA_type='LUA_TUSERDATA'
    typedef.LUA_pop=('({LUA_userdata_type} *)luaL_checkudata'
                     '({LUA_state_var}, 1, "{LUA_metadata}")')
    # typedef.LUA_push=None  # XXX create a userdata object with metatable
    # typedef.LUA_statements={}

    # allow forward declarations to avoid recursive headers
    typedef.forward=typedef.cpp_type


def extern_C(output, position):
    """Create extern "C" guards for C++
    """
    if position == 'begin':
        output.extend([
                '#ifdef __cplusplus',
                'extern "C" {',
                '#endif'
                ])
    else:
        output.extend([
                '#ifdef __cplusplus',
                '}',
                '#endif'
                ])


class WrapperMixin(object):
    """Methods common to all wrapping classes.
    """

#####

    def _init_splicer(self, splicers):
        self.splicers = splicers
        self.splicer_stack = [splicers]
        self.splicer_names = []
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

    def _create_splicer(self, name, out, default=[]):
        """Insert a splicer with *name* into list *out*.
        Use the splicer from the splicer_stack if it exists.
        This allows the user to replace the default text.
        TODO:
          Option to ignore splicer stack to generate original code
        """
        # The prefix is needed when two different sets of output
        # are being create and they are not in sync.
        # Creating methods and derived types together.
        show_splicer_comments = self.tree['options'].show_splicer_comments
        if show_splicer_comments:
            out.append('%s splicer begin %s%s' % (
                self.comment, self.splicer_path, name))
        out.extend(self.splicer_stack[-1].get(name, default))
        if show_splicer_comments:
            out.append('%s splicer end %s%s' % (
                self.comment, self.splicer_path, name))

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
        return typ + ' ' + (name or arg['name'])

#####

    def namespace(self, library, cls, position, output):
        if cls and 'namespace' in cls:
            namespace = cls['namespace']
            if namespace.startswith('-'):
                return
        else:
            namespace = library['namespace']
        if not namespace:
            return
        output.append('')
        if position == 'begin':
            for name in namespace.split():
                output.append('namespace %s {' % name)
        else:
            lst = namespace.split()
            lst.reverse()
            for name in lst:
                output.append('}  // namespace %s' % name)
#####

    def write_output_file(self, fname, directory, output):
        """
        fname  - file name
        directory - output directory
        output - list of lines to write
        """
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
        """
        Write the copyright from the input YAML file.
        """
        for line in self.tree.get('copyright', []):
            if line:
                fp.write(self.comment + ' ' + line + '\n')
            else:
                # convert None to blank line
                fp.write(self.comment + '\n')

    def write_lines(self, fp, lines):
        """ Write lines with indention and newlines.
        """
        for line in lines:
            if isinstance(line, int):
                self.indent += int(line)
            else:
                for subline in line.split("\n"):
                    if len(subline) == 0:
                        fp.write('\n')
                    elif subline[0] == '#':
                        # preprocessing directives work better in column 1
                        fp.write(subline)
                        fp.write('\n')
                    else:
                        fp.write('    ' * self.indent)
                        fp.write(subline)
                        fp.write('\n')

    def write_doxygen_file(self, output, fname, library, cls):
        """ Write a doxygen comment block for a file.
        """
        node = cls or library
        output.append(self.doxygen_begin)
        output.append(self.doxygen_cont + ' \\file %s' % fname)
        if cls:
            output.append(self.doxygen_cont +
                          ' \\brief Shroud generated wrapper for {} class'
                          .format(node['name']))
        else:
            output.append(self.doxygen_cont +
                          ' \\brief Shroud generated wrapper for {} library'
                          .format(node['library']))
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
    This used to be a dict but a class has better access semantics:
       i.attr vs d['attr']
    It also initializes default values to avoid  d.get('attr', default)
    """

    # Array of known keys with default values
    _order = (
        ('base', 'unknown'),      # Base type: 'string'
        ('forward', None),        # Forward declaration
        ('typedef', None),        # Initialize from existing type

        ('cpp_type', None),       # Name of type in C++
        ('cpp_to_c', '{cpp_var}'), # Expression to convert from C++ to C
        ('cpp_header', None),     # Name of C++ header file required for implementation
                                  # For example, if cpp_to_c was a function
        ('cpp_local_var', False), # True if c_to_cpp requires a local C variable

        ('c_type', None),         # Name of type in C
        ('c_header', None),       # Name of C header file required for type
        ('c_to_cpp', '{c_var}'),  # Expression to convert from C to C++
        ('c_statements', {}),
        ('c_return_code', None),

        ('f_c_args', None),       # List of argument names to F_C routine
        ('f_c_argdecl', None),    # List of declarations to F_C routine
        ('f_c_module', None),     # Fortran modules needed for interface  (dictionary)

        ('f_type', None),         # Name of type in Fortran
        ('f_c_type', None),       # Type for C interface
        ('f_to_c', None),         # Expression to convert from Fortran to C
        ('f_derived_type', None), # Fortran derived type name
        ('f_args', None),         # Argument in Fortran wrapper to call C.
        ('f_module', None),       # Fortran modules needed for type  (dictionary)
        ('f_return_code', None),
        ('f_cast', '{f_var}'),    # Expression to convert to type
                                  # e.g. intrinsics such as int and real
        ('f_statements', {}),
        ('f_helper', {}),         # helper functions to insert into module as PRIVATE

        ('result_as_arg', None),  # override fields when result should be treated as an argument

        # Python
        ('PY_format', 'O'),       # 'format unit' for PyArg_Parse
        ('PY_PyTypeObject', None), # variable name of PyTypeObject instance
        ('PY_PyObject', None),    # typedef name of PyObject instance
        ('PY_ctor', None),        # expression to create object.
                                  # ex. PyBool_FromLong({rv})
        ('PY_to_object', None),   # PyBuild - object'=converter(address)
        ('PY_from_object', None), # PyArg_Parse - status=converter(object, address);
        ('py_statements', {}),

        # Lua
        ('LUA_type', 'LUA_TNONE'),
        ('LUA_pop', 'POP'),
        ('LUA_push', 'PUSH'),
        ('LUA_statements', {}),
    )


    _keyorder, _valueorder = zip(*_order)

    # valid fields
    defaults = dict(_order)

    def __init__(self, name, **kw):
        self.name = name
#        for key, defvalue in self.defaults.items():
#            setattr(self, key, defvalue)
        self.__dict__.update(self.defaults)  # set all default values
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

    def __as_yaml__(self, indent, output):
        """Write out entire typedef as YAML.
        """
        as_yaml(self, self._keyorder, indent, output)

    def __export_yaml__(self, indent, output):
        """Write out a subset of a wrapped type.
        Other fields are set with typedef_wrapped_defaults.
        """
        as_yaml(self, [
            'base',
            'cpp_header',
            'cpp_type',
            'c_type',
            'c_header',
            'f_derived_type',
            'f_to_c',
            'f_module',
        ], indent, output)


class Options(object):
    """
    If attribute is not found, look in parent's.
    A replacement for a dictionary to allow obj.name syntax.
    It will automatically look in __parent for attribute if not found to allow
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
            d.update(self._to_dict())
        if self.__parent:
            self.__parent._to_full_dict(d)
        return d


def copy_function_node(node):
    """Create a copy of a function node to use with C++ template.
    """
    # Shallow copy everything
    new = node.copy()

    # Deep copy dictionaries
    for field in ['args', 'attrs', 'result']:
        new[field] = copy.deepcopy(node[field])

    # Add new Options in chain
    for field in ['fmt', 'options']:
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
    # Argument
    print("Test Typedef")
    a = Typedef('top', base='new_base')  # , bird='abcd')
    print(json.dumps(a, cls=ExpandedEncoder, sort_keys=True))

    print("FORMAT:", wformat("{a} {z} {c2} {z2}", lev1))

    print(json.dumps(lev0, cls=ExpandedEncoder, sort_keys=True))
    print(json.dumps(lev1, cls=ExpandedEncoder, sort_keys=True))
