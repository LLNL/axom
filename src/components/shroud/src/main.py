#!/bin/env python3

"""
generate language bindings
"""

#
# Annotate the YAML tree with additional internal fields
#  _decl            - generated declaration.  Includes computed attributes
#  _function_index  - sequence number function, used in lieu of a pointer
#  _generated       - who generated this function
#  _PTR_F_C_index   - Used by fortran wrapper to find index of C function to call
#  _PTR_C_CPP_index - Used by C wrapper to find index of C++ function to call
#  _subprogram      - subroutine or function
#
#
from __future__ import print_function

import argparse
import copy
import json
import os
import sys
import yaml

import util
import parse_decl
import splicer
import wrapc
import wrapf
import wrapp
import wrapl

# char functions cannot be wrapped directly in intel 15.  Instead the result is passed down
# as an argument from the Fortran wrapper to the C wrapper.
# Similar to how char * funtions are handled.
intel_15_fix = True

wformat = util.wformat

class Config(object):
    """A class to stash configuration values.
    """
    pass


class Schema(object):
    """
    Verify that the input dictionary has the correct fields.
    Create defaults for missing fields.


    check_schema
      check_classes
        check_class
          check_function
        check_class_depedencies
          check_function_dependencies
      check_functions
        check_function
    """
    def __init__(self, config):
        self.config = config
        self.fmt_stack = []

    def push_options(self, node):
        """ Push a new set of options.
        Copy current options, then update with new options.
        """
        new = util.Options(parent=self.options_stack[-1])
        if 'options' in node and \
                node['options'] is not None:
            if not isinstance(node['options'], dict):
                raise TypeError("options must be a dictionary")
            new.update(node['options'])
        self.options_stack.append(new)
        node['options'] = new
        return new

    def pop_options(self):
        self.options_stack.pop()

    def check_options_only(self, node):
        """Process an options only entry in a list.

        Return True if node only has options.
        node is assumed to be a dictionary.
        Update current set of options from node['options'].
        """
        if len(node) != 1:
            return False
        options = node.get('options', None)
        if not options:
            return False
        if not isinstance(options, dict):
            raise TypeError("options must be a dictionary")

        # replace current options
        new = util.Options(parent=self.options_stack[-1])
        new.update(node['options'])
        self.options_stack[-1] = new
        return True

    def push_fmt(self, node):
        fmt = util.Options(self.fmt_stack[-1])
        self.fmt_stack.append(fmt)
        node['fmt'] = fmt
        return fmt

    def pop_fmt(self):
        self.fmt_stack.pop()

    def check_schema(self):
        """ Routine to check entire schema of input tree"
        """
        node = self.config

        # default options
        def_options = util.Options(
            parent=None,
            debug=False,   # print additional debug info

            F_module_per_class=True,
            F_string_len_trim=True,
            F_force_wrapper=False,

            wrap_c       = True,
            wrap_fortran = True,
            wrap_python  = False,
            wrap_lua     = False,

            doxygen = True,       # create doxygen comments

            # blank for functions, set in classes.
            class_name_template = '{class_lower}_',

            C_header_filename_library_template = 'wrap{library}.h',
            C_impl_filename_library_template = 'wrap{library}.cpp',

            C_header_filename_class_template = 'wrap{cpp_class}.h',
            C_impl_filename_class_template = 'wrap{cpp_class}.cpp',

            C_name_template = '{C_prefix}{class_name}{underscore_name}{function_suffix}',

            # Fortran's names for C functions
            F_C_prefix = 'c_',
            F_C_name_template = '{F_C_prefix}{class_name}{underscore_name}{function_suffix}',

            F_name_impl_template ='{class_name}{underscore_name}{function_suffix}',

            F_name_method_template = '{underscore_name}{function_suffix}',
            F_name_generic_template = '{underscore_name}',

            F_module_name_library_template = '{library_lower}_mod',
            F_impl_filename_library_template = 'wrapf{library_lower}.f',

            F_module_name_class_template = '{class_lower}_mod',
            F_impl_filename_class_template = 'wrapf{cpp_class}.f',

            F_name_instance_get = 'get_instance',
            F_name_instance_set = 'set_instance',
            F_name_associated = 'associated',

            )
        wrapp.add_templates(def_options)
        wrapl.add_templates(def_options)

        if 'options' in node:
            def_options.update(node['options'])
        self.options_stack = [ def_options ]
        node['options'] = def_options

        fmt_library = node['fmt'] = util.Options(None)
        fmt_library.library       = node['library']
        fmt_library.library_lower = fmt_library.library.lower()
        fmt_library.library_upper = fmt_library.library.upper()
        fmt_library.function_suffix = ''   # assume no suffix
        fmt_library.overloaded    = False
        fmt_library.class_name = ''
        fmt_library.C_prefix      = def_options.get('C_prefix', fmt_library.library_upper[:3] + '_')
        fmt_library.F_C_prefix    = def_options['F_C_prefix']
        fmt_library.rv            = 'rv'  # return value
        if node['namespace']:
            fmt_library.namespace_scope = '::'.join(node['namespace'].split()) + '::'
        else:
            fmt_library.namespace_scope = ''
        util.eval_template(node, 'C_header_filename', '_library')
        util.eval_template(node, 'C_impl_filename', '_library')
        self.fmt_stack.append(fmt_library)

        # default some options based on other options
        # All class/methods and functions may go into this file or just functions.
        util.eval_template(node, 'F_module_name', '_library')
        util.eval_template(node, 'F_impl_filename', '_library')

        def_types = dict(
            void    = util.Typedef('void',
                c_type   = 'void',
                cpp_type = 'void',
#                fortran = 'subroutine',
                c_fortran = 'type(C_PTR)',
                f_type = 'type(C_PTR)',
                PY_ctor = 'PyCapsule_New({cpp_var}, NULL, NULL)',
                ),
            int    = util.Typedef('int',
                c_type    = 'int',
                cpp_type  = 'int',
                f_kind    = 'C_INT',
                f_cast    = 'int({f_var}, C_INT)',
                c_fortran = 'integer(C_INT)',
                f_type    = 'integer(C_INT)',
                f_module  = dict(iso_c_binding=['C_INT']),
                PY_format = 'i',
                LUA_type  = 'LUA_TNUMBER',
                LUA_pop   = 'lua_tointeger({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushinteger({LUA_state_var}, {c_var})',
                ),
            long   = util.Typedef('long',
                c_type    = 'long',
                cpp_type  = 'long',
                f_kind    = 'C_LONG',
                f_cast    = 'int({f_var}, C_LONG)',
                c_fortran = 'integer(C_LONG)',
                f_type    = 'integer(C_LONG)',
                f_module  = dict(iso_c_binding=['C_LONG']),
                PY_format = 'l',
                LUA_type  = 'LUA_TNUMBER',
                LUA_pop   = 'lua_tointeger({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushinteger({LUA_state_var}, {c_var})',
                ),
            size_t   = util.Typedef('size_t',
                c_type    = 'size_t',
                cpp_type  = 'size_t',
                c_header  = 'stdlib.h',
                f_kind    = 'C_SIZE_T',
                f_cast    = 'int({f_var}, C_SIZE_T)',
                c_fortran = 'integer(C_SIZE_T)',
                f_type    = 'integer(C_SIZE_T)',
                f_module  = dict(iso_c_binding=['C_SIZE_T']),
                PY_ctor   = 'PyInt_FromLong({rv})',
                LUA_type  = 'LUA_TNUMBER',
                LUA_pop   = 'lua_tointeger({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushinteger({LUA_state_var}, {c_var})',
                ),

            float   = util.Typedef('float',
                c_type    = 'float',
                cpp_type  = 'float',
                f_kind    = 'C_FLOAT',
                f_cast    = 'real({f_var}, C_FLOAT)',
                c_fortran = 'real(C_FLOAT)',
                f_type    = 'real(C_FLOAT)',
                f_module  = dict(iso_c_binding=['C_FLOAT']),
                PY_format = 'f',
                LUA_type  = 'LUA_TNUMBER',
                LUA_pop   = 'lua_tonumber({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushnumber({LUA_state_var}, {c_var})',
                ),
            double   = util.Typedef('double',
                c_type    = 'double',
                cpp_type  = 'double',
                f_kind    = 'C_DOUBLE',
                f_cast    = 'real({f_var}, C_DOUBLE)',
                c_fortran = 'real(C_DOUBLE)',
                f_type    = 'real(C_DOUBLE)',
                f_module  = dict(iso_c_binding=['C_DOUBLE']),
                PY_format = 'd',
                LUA_type  = 'LUA_TNUMBER',
                LUA_pop   = 'lua_tonumber({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushnumber({LUA_state_var}, {c_var})',
                ),

            bool   = util.Typedef('bool',
                c_type    = 'bool',
                cpp_type  = 'bool',
                f_kind    = 'C_BOOL',
                c_fortran = 'logical(C_BOOL)',

                f_type    = 'logical',
                f_module  = dict(iso_c_binding=['C_BOOL']),
                f_statements = dict(
                    intent_in = dict(
                        declare = [
                            'logical(C_BOOL) {c_var}',
                            ],
                        pre_call = [
                            '{c_var} = {f_var}  ! coerce to C_BOOL',
                            ],
                        ),
                    result = dict(
                        # The wrapper is needed to convert bool to logical
                        need_wrapper = True,
                        ),
                    ),

                py_statements = dict(
                    intent_in = dict(
                        post_parse = [
                            '{cpp_var} = PyObject_IsTrue({py_var});',
                            ],
                        ),
                    ),

#XXX            PY_format = 'p',  # Python 3.3 or greater
                PY_ctor   = 'PyBool_FromLong({rv})',
                PY_PyTypeObject = 'PyBool_Type',
                LUA_type  = 'LUA_TBOOLEAN',
                LUA_pop   = 'lua_toboolean({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushboolean({LUA_state_var}, {c_var})',
                ),

            # implies null terminated string
            char = util.Typedef('char',
                cpp_type = 'char',
#                cpp_header = '<string>',
#                cpp_to_c = '{cpp_var}.c_str()',  # . or ->

                c_type   = 'char',    # XXX - char *

                c_statements = dict(
                    intent_in = dict(
                        cpp_local_var = True,
                        pre_call = [
                            'int {c_var_len} = strlen({c_var});',
                            'char * {cpp_var} = new char [{c_var_len} + 1];',
                            'std::strncpy({cpp_var}, {c_var}, {c_var_len});',
                            '{cpp_var}[{c_var_len}] = \'\\0\';'
                            ],
                        pre_call_trim = [
                            'char * {cpp_var} = new char [{c_var_len} + 1];',
                            'std::strncpy({cpp_var}, {c_var}, {c_var_len});',
                            '{cpp_var}[{c_var_len}] = \'\\0\';'
                            ],
                        post_call = [
                            'delete [] {cpp_var};'
                            ],
                        ),
                    intent_out = dict(
                        cpp_local_var = True,
                        pre_call = [
                            'char * {cpp_var} = new char [{c_var_num} + 1];',
                            ],
                        post_call = [
                            'asctoolkit::shroud::FccCopy({c_var}, {c_var_len}, {cpp_val});',
                            'delete [] {cpp_var};',
                            ],
                        cpp_header = 'shroudrt.hpp',
                        ),
                    result = dict(
                        post_call = [
                            'asctoolkit::shroud::FccCopy({c_var}, {c_var_len}, {cpp_val});',
                            ],
                        cpp_header = 'shroudrt.hpp',
                        ),
                    ),

                c_fortran  = 'character(kind=C_CHAR)',
                f_type     = 'character(*)',
##                f_args = 'trim({var}) // C_NULL_CHAR',
#                f_module = dict(iso_c_binding = [ 'C_NULL_CHAR' ]),
                f_module = dict(iso_c_binding=None),
#                f_return_code = '{F_result} = fstr({F_C_name}({F_arg_c_call_tab}))',
                PY_format = 's',
                PY_ctor = 'PyString_FromString({c_var})',
                LUA_type  = 'LUA_TSTRING',
                LUA_pop   = 'lua_tostring({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushstring({LUA_state_var}, {c_var})',
                base = 'string',
                ),

            # char scalar
            char_scalar = util.Typedef('char_scalar',
                cpp_type = 'char',
#                cpp_header = '<string>',
#                cpp_to_c = '{cpp_var}.c_str()',  # . or ->

                c_type   = 'char',    # XXX - char *

                c_fortran  = 'character(kind=C_CHAR)',
                f_type     = 'character',
##                f_args = 'trim({var}) // C_NULL_CHAR',
#                f_module = dict(iso_c_binding = [ 'C_NULL_CHAR' ]),
                f_module = dict(iso_c_binding=None),
#                f_return_code = '{F_result} = fstr({F_C_name}({F_arg_c_call_tab}))',
                PY_format = 's',
                PY_ctor = 'PyString_FromString({c_var})',
                LUA_type  = 'LUA_TSTRING',
                LUA_pop   = 'lua_tostring({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushstring({LUA_state_var}, {c_var})',
##                base = 'string',
                ),

            # C++ std::string
            string = util.Typedef('string',
                cpp_type = 'std::string',
                cpp_header = '<string>',
                cpp_to_c = '{cpp_var}.c_str()',  # . or ->

                c_type   = 'char',    # XXX - char *

                c_statements = dict(
                    intent_in = dict(
                        cpp_local_var = True,
                        pre_call = [
                            '{C_const}std::string {cpp_var}({c_var});'
                            ],
                        pre_call_trim = [
                            '{C_const}std::string {cpp_var}({c_var}, {c_var_trim});'
                            ],
                    ),
                    intent_out = dict(
                        post_call = [
                            'asctoolkit::shroud::FccCopy({c_var}, {c_var_len}, {cpp_val});',
                            ],
                        cpp_header = 'shroudrt.hpp'
                        ),
                    result = dict(
                        post_call = [
                            'asctoolkit::shroud::FccCopy({c_var}, {c_var_len}, {cpp_val});',
                            ],
                        cpp_header = 'shroudrt.hpp'
                        ),
                    ),

                c_fortran  = 'character(kind=C_CHAR)',
                f_type     = 'character(*)',
##                f_args = 'trim({var}) // C_NULL_CHAR',
#                f_module = dict(iso_c_binding = [ 'C_NULL_CHAR' ]),
                f_module = dict(iso_c_binding=None),
#                f_return_code = '{F_result} = fstr({F_C_name}({F_arg_c_call_tab}))',

                py_statements = dict(
                    intent_in = dict(
                        cpp_local_var = True,
                        post_parse = [
                            '{C_const}std::string {cpp_var}({c_var});'
                            ],
                        ),
                    ),
                PY_format = 's',
                PY_ctor = 'PyString_FromString({c_var})',
                LUA_type  = 'LUA_TSTRING',
                LUA_pop   = 'lua_tostring({LUA_state_var}, {LUA_index})',
                LUA_push  = 'lua_pushstring({LUA_state_var}, {c_var})',
                base = 'string',
                ),

            MPI_Comm = util.Typedef(
                'MPI_Comm',
                cpp_type = 'MPI_Comm',
                c_header = 'mpi.h',
                c_type = 'MPI_Fint',
                c_fortran = 'integer(C_INT)',  # usually, MPI_Fint will be equivalent to int
                f_type = 'integer',
                cpp_to_c = 'MPI_Comm_c2f({cpp_var})',
                c_to_cpp = 'MPI_Comm_f2c({c_var})',
                ),
            )

        if intel_15_fix:
            # Copy C++ function result into C result argument.
            def_types['char_scalar'].c_statements = dict(
                result = dict(
                    post_call = [
                        '// {c_var_len} is always 1, test to silence warning about unused variable',
                        'if ({c_var_len} == 1) *{c_var} = {cpp_val};',
                        ],
                    ),
                )

        # aliases
        def_types['std::string']     = def_types['string']
        def_types['integer(C_INT)']  = def_types['int']
        def_types['integer(C_LONG)'] = def_types['long']
        def_types['real(C_FLOAT)']   = def_types['float']
        def_types['real(C_DOUBLE)']  = def_types['double']

        # pure fortran string
        tmp = def_types['string'].clone_as('string_result_fstr')
        tmp.update(dict(
                f_return_code = '{F_result} = fstr({F_C_name}({F_arg_c_call_tab}))',
                f_helper = dict(f_return_code=dict(fstr=True)),
                ))
        def_types[tmp.name] = tmp

        types_dict = node.get('types', None)
        if types_dict is not None:
            if not isinstance(types_dict, dict):
                raise TypeError("types must be a dictionary")
            for key, value in types_dict.items():
                if not isinstance(value, dict):
                    raise TypeError("types '%s' must be a dictionary" % key)

                if 'typedef' in value:
                    copy_type = value['typedef']
                    orig = def_types.get(copy_type, None)
                    if not orig:
                        raise RuntimeError("No type for typedef %s" % copy_type)
                    def_types[key] = util.Typedef(key)
                    def_types[key].update(def_types[copy_type]._to_dict())


                if key in def_types:
                    def_types[key].update(value)
                else:
                    def_types[key] = util.Typedef(key, **value)

        patterns = node.setdefault('patterns', [])

        node['types'] = def_types
        self.typedef = node['types']

        classes = node.setdefault('classes', [])
        self.check_classes(classes)
        self.check_functions(node, 'functions')

    def check_classes(self, node):
        if not isinstance(node, list):
            raise TypeError("classes must be a list")
        for cls in node:
            if not isinstance(cls, dict):
                raise TypeError("classes[n] must be a dictionary")
            self.check_class(cls)

    def check_class(self, node):
        if 'name' not in node:
            raise RuntimeError('Expected name for class')
        name = node['name']

        # default cpp_header to blank
        if 'cpp_header' not in node:
            node['cpp_header'] = ''
        if node['cpp_header'] is None:
            # YAML turns blank strings into None
            node['cpp_header'] = ''

        options = self.push_options(node)
        fmt_class = self.push_fmt(node)
        fmt_class.cpp_class = name
        fmt_class.class_lower = name.lower()
        fmt_class.class_upper = name.upper()
        if 'C_prefix' in options:
            fmt_class.C_prefix = options.C_prefix
        if 'F_C_prefix' in options:
            fmt_class.F_C_prefix = options.F_C_prefix
        util.eval_template(node, 'class_name')
        
        if options.F_module_per_class:
            util.eval_template(node, 'F_module_name', '_class')
            util.eval_template(node, 'F_impl_filename', '_class')

        # Only one file per class for C.
        util.eval_template(node, 'C_header_filename', '_class')
        util.eval_template(node, 'C_impl_filename', '_class')

        self.check_functions(node, 'methods')
        self.pop_fmt()
        self.pop_options()

    def check_function(self, node):
        """ Make sure necessary fields are present for a function.
        """
        self.push_options(node)
        fmt_func = self.push_fmt(node)

#        func = util.FunctionNode()
#        func.update(node)
#        func.dump()

        node['attrs'] = {}
        node['args'] = []

        if 'decl' in node:
            # parse decl and add to dictionary
            values = parse_decl.check_decl(node['decl'])
            util.update(node, values)  # recursive update
        if 'function_suffix' in node and node['function_suffix'] is None:
            # YAML turns blanks strings into None
            node['function_suffix'] = ''
        if 'default_arg_suffix' in node:
            default_arg_suffix = node['default_arg_suffix']
            if not isinstance(default_arg_suffix, list):
                raise RuntimeError('default_arg_suffix must be a list')
            for i, value in enumerate(node['default_arg_suffix']):
                if value is None:
                    node['default_arg_suffix'][i] = '' # YAML turns blanks strings to None
        if 'result' not in node:
            raise RuntimeError("Missing result")
        result = node['result']
        if 'name' not in result:
            raise RuntimeError("Missing result.name")
        if 'type' not in result:
            raise RuntimeError("Missing result.type")
        if 'attrs' not in result:
            result['attrs'] = {}

        result = node['result']

        fmt_func.function_name   = result['name']
        fmt_func.underscore_name = util.un_camel(result['name'])

        # docs
        self.pop_fmt()
        self.pop_options()

    def check_functions(self, node, member):
        """ check functions.

        Create a new list without the options only entries.
        """
        functions = node.get(member, [])

        if not isinstance(functions, list):
            raise TypeError("functions must be a list")
        only_functions = []
        for func in functions:
            if self.check_options_only(func):
                continue
            self.check_function(func)
            only_functions.append(func)
        node[member] = only_functions


class GenFunctions(object):
    """
    Generate types from class.
    Generate functions based on overload/template/generic/attributes
    Computes fmt.function_suffix.
    """

    def __init__(self, tree, config):
        self.tree = tree    # json tree
        self.config = config
        self.typedef = tree['types']

    def gen_library(self):
        tree = self.tree

        # Order of creating.
        # Each is given a _function_index when created.
        tree['function_index'] = self.function_index = []

        for cls in tree['classes']:
            cls['methods'] = self.define_function_suffix(cls['methods'])
        tree['functions'] = self.define_function_suffix(tree['functions'])

        for cls in tree['classes']:
            self.check_class_dependencies(cls)

    def append_function_index(self, node):
        """append to function_index, set index into node.
        """
        ilist = self.function_index
        node['_function_index'] = len(ilist)
#        node['fmt'].function_index = str(len(ilist)) # debugging
        ilist.append(node)

    def define_function_suffix(self, functions):
        """ Look for functions with the same name.
        Return a new list with overloaded function inserted.
        """
        cpp_overload = {}
        for function in functions:
            if 'function_suffix' in function:
                function['fmt'].function_suffix = function['function_suffix']
            self.append_function_index(function)
            cpp_overload. \
                setdefault(function['result']['name'], []). \
                append(function['_function_index'])

        # keep track of which function are overloaded in C++.
        for key, value in cpp_overload.items():
            if len(value) > 1:
                for index in value:
                    self.function_index[index]['_cpp_overload'] = value

        # Create additional functions needed for wrapping
        ordered_functions = []
        for method in functions:
            if '_has_default_arg' in method:
                self.has_default_args(method, ordered_functions)
            ordered_functions.append(method)
            if 'cpp_template' in method:
                method['_overloaded'] = True
                self.template_function(method, ordered_functions)

        # Look for overloaded functions
        overloaded_functions = {}
        for function in ordered_functions:
#            if not function['options'].wrap_c:
#                continue
            if 'cpp_template' in function:
                continue
            overloaded_functions.setdefault(
                function['result']['name'], []).append(function)

        # look for function overload and compute function_suffix
        for mname, overloads in overloaded_functions.items():
            if len(overloads) > 1:
                for i, function in enumerate(overloads):
                    function['_overloaded'] = True
                    if not function['fmt'].inlocal('function_suffix'):
                        function['fmt'].function_suffix =  '_%d' % i

        # Create additional C bufferify functions.
        ordered3 = []
        for method in ordered_functions:
            ordered3.append(method)
            self.string_to_buffer_and_len(method, ordered3)

        # Create multiple generic Fortran wrappers to call a
        # single C functions
        ordered4 = []
        for method in ordered3:
            ordered4.append(method)
            if not method['options'].wrap_fortran:
                continue
            if 'fortran_generic' in method:
                method['_overloaded'] = True
                self.generic_function(method, ordered4)

        self.gen_functions_decl(ordered4)

        return ordered4

    def template_function(self, node, ordered_functions):
        """ Create overloaded functions for each templated argument.
        """
        if len(node['cpp_template']) != 1:
            # In the future it may be useful to have multiple templates
            # That the would start creating more permutations
            raise NotImplemented("Only one cpp_templated type for now")
        for typename, types in node['cpp_template'].items():
            for type in types:
                new = util.copy_function_node(node)
                ordered_functions.append(new)
                self.append_function_index(new)

                new['_generated'] = 'cpp_template'
                fmt = new['fmt']
                fmt.function_suffix = fmt.function_suffix + '_' + type
                del new['cpp_template']
                options = new['options']
                options.wrap_c = True
                options.wrap_fortran = True
                options.wrap_python = False
                options.wrap_lua = False
                # Convert typename to type
                fmt.CPP_template = '<%s>' %  type
                if new['result']['type'] == typename:
                    new['result']['type'] = type
                    fmt.CPP_return_templated = True
                for arg in new['args']:
                    if arg['type'] == typename:
                        arg['type'] = type

        # Do not process templated node, instead process
        # generated functions above.
        options = node['options']
        options.wrap_c = False
        options.wrap_fortran = False
        options.wrap_python = False
        options.wrap_lua = False

    def generic_function(self, node, ordered_functions):
        """ Create overloaded functions for each generic method.
        """
        if len(node['fortran_generic']) != 1:
            # In the future it may be useful to have multiple generic arguments
            # That the would start creating more permutations
            raise NotImplemented("Only one generic arg for now")
        for argname, types in node['fortran_generic'].items():
            for type in types:
                new = util.copy_function_node(node)
                ordered_functions.append(new)
                self.append_function_index(new)

                new['_generated'] = 'fortran_generic'
                new['_PTR_F_C_index'] = node['_function_index']
                fmt = new['fmt']
                # XXX append to existing suffix
                fmt.function_suffix = fmt.function_suffix + '_' + type
                del new['fortran_generic']
                options = new['options']
                options.wrap_c = False
                options.wrap_fortran = True
                options.wrap_python = False
                options.wrap_lua = False
                # Convert typename to type
                for arg in new['args']:
                    if arg['name'] == argname:
                        # Convert any typedef to native type with f_type
                        argtype = arg['type']
                        typedef = self.typedef[argtype]
                        typedef = self.typedef[typedef.f_type]
                        if not typedef.f_cast:
                            raise RuntimeError("unable to cast type %s in fortran_generic" % arg['type'])
                        arg['type'] = type

        # Do not process templated node, instead process
        # generated functions above.
        options = node['options']
#        options.wrap_c = False
        options.wrap_fortran = False
#        options.wrap_python = False

    def has_default_args(self, node, ordered_functions):
        """
        For each function which has a default argument, generate
        a version for each possible call.
          void func(int i = 0, int j = 0)
        generates
          void func()
          void func(int i)
          void func(int i, int j)
        """
        default_funcs = []

        default_arg_suffix = node.get('default_arg_suffix', [])
        ndefault = 0

        min_args = 0
        for i, arg in enumerate(node['args']):
            if 'default' not in arg['attrs']:
                min_args += 1
                continue
            new = util.copy_function_node(node)
            self.append_function_index(new)
            new['_generated'] = 'has_default_arg'
            del new['args'][i:]  # remove trailing arguments
#            try:
            del new['_has_default_arg']
#            except:
#                pass
            options = new['options']
            options.wrap_c = True
            options.wrap_fortran = True
            options.wrap_python = False
            options.wrap_lua = False
            fmt = new['fmt']
            try:
                fmt.function_suffix = default_arg_suffix[ndefault]
            except IndexError:
#               fmt.function_suffix = fmt.function_suffix + '_nargs%d' % (i + 1)
                pass
            default_funcs.append(new['_function_index'])
            ordered_functions.append(new)
            ndefault += 1

        # keep track of generated default value functions
        node['_default_funcs'] = default_funcs
        node['_nargs'] = ( min_args, len(node['args']) )
        # The last name calls with all arguments (the original decl)
        try:
            node['fmt'].function_suffix = default_arg_suffix[ndefault]
        except IndexError:
            pass

    def string_to_buffer_and_len(self, node, ordered_functions):
        """ Check if function has any string arguments and will be wrapped by Fortran.
        If so then create a new C function that will convert string arguments into 
        a buffer and length.
        """
        options = node['options']

        # If a C++ function returns a std::string instance, the default wrapper
        # will not compile since the wrapper will be declared as char.
        # It will also want to return the c_str of a stack variable.
        # Warn and turn off the wrapper.
        result = node['result']
        result_type = result['type']
        try:
            result_typedef = self.typedef[result_type]
        except KeyError:
            # wrapped classes have not been added yet.  Only care about string here.
            result_typedef = None
        attrs = result['attrs']
        result_is_ptr = (attrs.get('ptr', False) or
                         attrs.get('reference', False))
        if result_typedef and result_typedef.base == 'string' and \
                result_type != 'char' and \
                not result_is_ptr:
            options.wrap_c = False
#            options.wrap_fortran = False
            self.config.log.write("Skipping %s, unable to create C wrapper for function returning std::string instance (must return a pointer or reference). \n" % ( result['name']) )

        if options.wrap_fortran is False:
            return
        if options.F_string_len_trim is False:
            return

        # Is result or any argument a string?
        has_string_arg = False
        for arg in node['args']:
            argtype = arg['type']
            if self.typedef[argtype].base == 'string':
                attrs = arg['attrs']
                is_ptr = (attrs.get('ptr', False) or
                          attrs.get('reference', False))
                if is_ptr:
                    has_string_arg = True
                    # Force len attribute when intent is OUT
                    # so the wrapper will know how much space can be written to.
                    intent = attrs['intent']
                    if intent in ['out', 'inout']:
                        attrs['len'] = 'N' + arg['name']
                else:
                    arg['type'] = 'char_scalar'

        has_string_result = False
        result_as_arg = ''  # only applies to string functions
        is_pure = node['attrs'].get('pure', False)
        if result_typedef.base == 'string':
            if result_type == 'char' and not result_is_ptr:
                if intel_15_fix:
                    result['attrs']['len'] = 1
                    has_string_result = True
                    result_as_arg = options.get('F_string_result_as_arg', '')
                    result_name = result_as_arg or 'SH_F_rv'
                result['type'] = 'char_scalar'
            else:
                has_string_result = True
                result_as_arg = options.get('F_string_result_as_arg', '')
                result_name = result_as_arg or 'SH_F_rv'

        if not (has_string_result or has_string_arg):
            return

#        options = node['options']
#        options.wrap_fortran = False
        # Preserve wrap_c.  This keep a version which accepts char * arguments.

        # Create a new C function and change arguments to add len_trim attribute
        C_new = util.copy_function_node(node)
        ordered_functions.append(C_new)
        self.append_function_index(C_new)

        C_new['_generated'] = 'string_to_buffer_and_len'
        C_new['_error_pattern_suffix'] = '_as_buffer'
        fmt = C_new['fmt']
        fmt.function_suffix = fmt.function_suffix + '_bufferify'

        options = C_new['options']
        options.wrap_c = True
        options.wrap_fortran = False
        options.wrap_python = False
        options.wrap_lua = False
        C_new['_PTR_C_CPP_index'] = node['_function_index']

        newargs = []
        for arg in C_new['args']:
            argtype = arg['type']
            if self.typedef[argtype].base == 'string':
                # strings passed in need len_trim
                # strings returned need len
                intent = arg['attrs']['intent']
                if intent in ['in', 'inout']:
                    arg['attrs']['len_trim'] = 'L' + arg['name']
                if intent in ['out', 'inout']:
                    arg['attrs']['len'] = 'N' + arg['name']

        if has_string_result:
            # Add additional argument to hold result
            result_as_string = copy.deepcopy(result)
            result_as_string['name'] = result_name
            attrs = result_as_string['attrs']
            attrs['const'] = False
            attrs['len'] = 'L' + result_name
            attrs['intent'] = 'out'
            attrs['_is_result'] = True
            if not result_is_ptr:
                attrs['ptr'] = True
                attrs['reference'] = False
            C_new['args'].append(result_as_string)

            # convert to subroutine
            C_new['_subprogram'] = 'subroutine'
            result = C_new['result']
            result['type'] = 'void'
            attrs = result['attrs']
            attrs['const'] = False
            attrs['ptr'] = False
            attrs['reference'] = False

        if is_pure:
            # Return a character(*) function
            # Create a fortran function with a different result type
            # so fstr will be called on result.
            F_new = util.copy_function_node(node)
            ordered_functions.append(F_new)
            self.append_function_index(F_new)

            F_new['result']['type'] = 'string_result_fstr'
            F_new['_PTR_F_C_index'] = node['_function_index']
            options = F_new['options']
            options.wrap_c       = False
            options.wrap_fortran = True
            options.wrap_python = False
            options.wrap_lua = False

            # Do not wrap original function (Has a different result type)
            node['options'].wrap_fortran = False

        elif result_as_arg:
            # Create Fortran function without bufferify function_suffix but
            # with len attributes on string arguments.
            F_new = util.copy_function_node(C_new)
            ordered_functions.append(F_new)
            self.append_function_index(F_new)

            # Fortran function should wrap the new C function
            F_new['_PTR_F_C_index'] = C_new['_function_index']
            options = F_new['options']
            options.wrap_c       = False
            options.wrap_fortran = True
            options.wrap_python = False
            options.wrap_lua = False
            # Do not add '_bufferify'
            F_new['fmt'].function_suffix = node['fmt'].function_suffix

            # Do not wrap original function (does not have result argumument)
            node['options'].wrap_fortran = False
        else:
            # Fortran function may call C subroutine if string result
            node['_PTR_F_C_index'] = C_new['_function_index']

    def check_class_dependencies(self, node):
        """
        Check used_types and find which header and module files
        to use for this class
        """
        # keep track of types which are used by methods arguments
        used_types = {}
        for method in node['methods']:
            self.check_function_dependencies(method, used_types)

        modules = {}
        for typ in used_types.values():
            if typ.f_module:
                for mname, only in typ.f_module.items():
                    module = modules.setdefault(mname, {})
                    if only:  # Empty list means no ONLY clause
                        for oname in only:
                            module[oname] = True

        # Always add C_PTR, needed for class F_derived_member
        modules.setdefault('iso_c_binding', {})['C_PTR'] = True

        F_modules = []  # array of tuples ( name, (only1, only2) )
        for mname in sorted(modules):
            F_modules.append((mname, sorted(modules[mname])))
        node['F_module_dependencies'] = F_modules

    def check_function_dependencies(self, node, used_types):
        """Record which types are used by a function.
        """
        if 'cpp_template' in node:
            # The templated type will raise an error.
            # XXX - Maybe dummy it out
            # XXX - process templated types
            return
        result = node['result']
        rv_type = result['type']
        if rv_type not in self.typedef:
            raise RuntimeError("Unknown type {} for function decl: {}".format(rv_type, node['decl']))
        result_typedef = self.typedef[rv_type]
        # XXX - make sure it exists
        used_types[result['type']] = result_typedef
        for arg in node['args']:
            argtype = arg['type']
            if argtype in self.typedef:
                used_types[arg['type']] = self.typedef[argtype]
            else:
                raise RuntimeError("%s not defined" % argtype)

    _skip_annotations = [ 'const', 'ptr', 'reference' ]
    def gen_annotations_decl(self, attrs, decl):
        """Append annotations from attrs onto decl in sorted order.
        Skip some that are already handled.
        """
        keys = attrs.keys()
        keys.sort()
        for key in keys:
            if key[0] == '_':  # internal attribute
                continue
            if key in self._skip_annotations:
                continue
            value = attrs[key]
            if value is True:
                decl.append('+' + key)
            elif value is False:
                pass
#                decl.append('-' + key)
            else:
                decl.append('+%s(%s)' % (key, value) )

    def gen_arg_decl(self, arg, decl):
        """ Generate declaration for a single arg (or result)
        """
        attrs = arg['attrs']
        if attrs.get('const', False):
            decl.append('const ')
        decl.append(arg['type'] + ' ')
        if attrs.get('ptr', False):
            decl.append('* ')
        if attrs.get('reference', False):
            decl.append('& ')
        decl.append(arg['name'])

    def gen_functions_decl(self, functions):
        """ Generate _decl for generated all functions.
        """
        for node in functions:
            decl = [ ]
            self.gen_arg_decl(node['result'], decl)

            if node['args']:
                decl.append('(')
                for arg in node['args']:
                    self.gen_arg_decl(arg, decl)
                    self.gen_annotations_decl(arg['attrs'], decl)
                    decl.append(', ')
                decl[-1] = ')'
            else:
                decl.append('()')

            attrs = node['attrs']
            if attrs.get('const', False):
                decl.append(' const')
            self.gen_annotations_decl(attrs, decl)

            node['_decl'] = ''.join(decl)


class VerifyAttrs(object):
    """
    Check attributes and set some defaults.
    Generate types for classes.
    """
    def __init__(self, tree, config):
        self.tree = tree    # json tree
        self.config = config
        self.typedef = tree['types']

    def verify_attrs(self):
        tree = self.tree

        for cls in tree['classes']:
            self.create_class_typedef(cls)

        for cls in tree['classes']:
            for func in cls['methods']:
                self.check_arg_attrs(func)

        for func in tree['functions']:
            self.check_arg_attrs(func)

    def create_class_typedef(self, cls):
        # create typedef for each class before generating code
        # this allows classes to reference each other
        name = cls['name']
        fmt_class = cls['fmt']

        if name not in self.typedef:
#            unname = util.un_camel(name)
            unname = name.lower()
            cname = fmt_class.C_prefix + unname
            self.typedef[name] = util.Typedef(
                name,
                cpp_type = name,
                cpp_to_c = 'static_cast<{C_const}%s *>(static_cast<{C_const}void *>({cpp_var}))' % cname,
                c_type = cname,
                # opaque pointer -> void pointer -> class instance pointer
                c_to_cpp = 'static_cast<{C_const}%s{ptr}>(static_cast<{C_const}void *>({c_var}))' % name,
                c_fortran = 'type(C_PTR)',
                f_type = 'type(%s)' % unname,
                f_derived_type = unname,
                f_args = '{c_var}%{F_derived_member}',
                # XXX module name may not conflict with type name
                f_module = {fmt_class.F_module_name:[unname]},

                # return from C function
#                f_c_return_decl = 'type(CPTR)' % unname,
                f_return_code = '{F_result}%{F_derived_member} = {F_C_name}({F_arg_c_call_tab})',

                py_statements = dict(
                    intent_in = dict(
                        post_parse = [
                            '{cpp_var} = {py_var} ? {py_var}->{BBB} : NULL;',
                            ],
                        ),
                    intent_out = dict(
                        ctor = [
                            '{PyObject} * {py_var} = PyObject_New({PyObject}, &{PyTypeObject});',
                            '{py_var}->{BBB} = {cpp_var};',
                            ]
                        ),
                    ),
#                PY_ctor = 'PyObject_New({PyObject}, &{PyTypeObject})',

                LUA_type  = 'LUA_TUSERDATA',
                LUA_pop = '({LUA_userdata_type} *)luaL_checkudata({LUA_state_var}, 1, "{LUA_metadata}")',
#                LUA_push = None,  # XXX create a userdata object with metatable
#                LUA_statements={},

                # allow forward declarations to avoid recursive headers
                forward = name,
                base = 'wrapped',
                )

        typedef = self.typedef[name]
        fmt_class.C_type_name = typedef.c_type

    def check_arg_attrs(self, node):
        """Regularize attributes
        intent: lower case, no parens, must be in, out, or inout
        value: if pointer, default to False (pass-by-reference;
               else True (pass-by-value).
        """
        options = node['options']
        if not options.wrap_fortran and not options.wrap_c:
            return

        # cache subprogram type
        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)
        #  'void'=subroutine   'void *'=function
        if result_type == 'void' and not result_is_ptr:
            node['_subprogram'] = 'subroutine'
        else:
            node['_subprogram'] = 'function'

        found_default = False
        for arg in node['args']:
            argname = arg['name']
            typename = arg['type']
            typedef = self.typedef.get(typename, None)
            if typedef is None:
                cpp_template = node.get('cpp_template', {})
                if typename not in cpp_template:
                    raise RuntimeError("No such type %s: %s %s" % (
                            typename, typename, argname))
                    # XXX print decl

            attrs = arg['attrs']
            is_ptr = (attrs.get('ptr', False) or
                      attrs.get('reference', False))

            # intent
            intent = attrs.get('intent', None)
            if intent is None:
                if not is_ptr:
                    attrs['intent'] = 'in'
                elif attrs.get('const', False):
                    attrs['intent'] = 'in'
                elif typedef.base == 'string':
                    attrs['intent'] = 'inout'
                else:
                    # void *
                    attrs['intent'] = 'in' # XXX must coordinate with VALUE
            else:
                intent = intent.lower()
                if intent[0] == '(' and intent[-1] == ')':
                    intent = intent[1:-1]
                if intent in ['in', 'out', 'inout']:
                    attrs['intent'] = intent
                else:
                    raise RuntimeError(
                        "Bad value for intent: " + attrs['intent'])

            # value
            value = attrs.get('value', None)
            if value is None:
                if is_ptr:
                    if typedef.c_fortran == 'type(C_PTR)':
                        # This causes Fortran to dereference the C_PTR
                        # Otherwise a void * argument becomes void **
                        attrs['value'] = True
                    else:
                        attrs['value'] = False
                else:
                    attrs['value'] = True

            # dimension
            dimension = attrs.get('dimension', None)
            if dimension:
                if attrs.get('value', False):
                    raise RuntimeError("argument must not have value=True")
                if not is_ptr:
                    raise RuntimeError("dimension attribute can only be used on pointer and references")
                if dimension is True:
                    # No value was provided, provide default
                    attrs['dimension'] = '(*)'

            if 'default' in attrs:
                found_default = True
                node['_has_default_arg'] = True
            elif found_default is True:
                raise RuntimeError("Expected default value for %s" % argname)

            # compute argument names for some attributes
            # XXX make sure they don't conflict with other names
            len_name = attrs.get('len', False)
            if len_name is True:
                attrs['len'] = 'L' + argname
            len_name = attrs.get('len_trim', False)
            if len_name is True:
                attrs['len'] = 'L' + argname
#        if typedef.base == 'string':


class Namify(object):
    """Compute names of functions in library.
    Need to compute F_name and C_F_name since they interact.
    Compute all C names first, then Fortran.
    A Fortran function may call a generated C function via
    _PTR_F_C_index
    Also compute number which may be controlled by options.    

    C_name - Name of C function
    F_C_name - Fortran function for C interface
    F_name_impl - Name of Fortran function implementation
    """
    def __init__(self, tree, config):
        self.tree = tree    # json tree
        self.config = config

    def name_library(self):
        options = self.tree['options']
        fmt_library = self.tree['fmt']
        fmt_library.C_this = options.get('C_this', 'self')

        fmt_library.F_this = options.get('F_this', 'obj')
        fmt_library.F_result = options.get('F_result', 'rv')
        fmt_library.F_derived_member = options.get('F_derived_member', 'voidptr')

        self.name_language(self.name_function_c)
        self.name_language(self.name_function_fortran)

    def name_language(self, handler):
        tree = self.tree
        for cls in tree['classes']:
            for func in cls['methods']:
                handler(cls, func)

            options = cls['options']
            fmt_class = cls['fmt']
            if 'F_this' in options:
                fmt_class.F_this = options.F_this

        for func in tree['functions']:
            handler(None, func)

    def name_function_c(self, cls, node):
        options = node['options']
        if not options.wrap_c:
            return
        fmt_func = node['fmt']
        
        util.eval_template(node, 'C_name')
        util.eval_template(node, 'F_C_name')
        fmt_func.F_C_name = fmt_func.F_C_name.lower()

        if 'C_this' in options:
            fmt_func.C_this = options.C_this

    def name_function_fortran(self, cls, node):
        """ Must process C functions for to generate their names.
        """
        options = node['options']
        if not options.wrap_fortran:
            return
        fmt_func = node['fmt']

        util.eval_template(node, 'F_name_impl')
        util.eval_template(node, 'F_name_method')
        util.eval_template(node, 'F_name_generic')

        if 'F_this' in options:
            fmt_func.F_this = options.F_this
        if 'F_result' in options:
            fmt_func.F_result = options.F_result


def main():

    appname = 'shroud'
    appver = '0.1'

    parser = argparse.ArgumentParser(
        prog=appname,
        description="""Create Fortran or Python wrapper for a C++ library.
""")
    parser.add_argument('--version', action='version', version=appver)
    parser.add_argument('--outdir', default='',
                        help='Directory for output files')
    parser.add_argument('--outdir-c-fortran', default='',
                        dest='outdir_c_fortran',
                        help='Directory for C/Fortran wrapper output files, ' +
                        'overrides --outdir')
    parser.add_argument('--outdir-python', default='',
                        dest='outdir_python',
                        help='Directory for Python wrapper output files, ' +
                        'overrides --outdir')
    parser.add_argument('--outdir-lua', default='',
                        dest='outdir_lua',
                        help='Directory for Lua wrapper output files, ' +
                        'overrides --outdir')
    parser.add_argument('--logdir', default='',
                        help='Directory for log files')
    parser.add_argument('--cfiles', default='',
                        help='Output file with list of C and C++ files created')
    parser.add_argument('--ffiles', default='',
                        help='Output file with list of Fortran created')
    parser.add_argument('--path', default=[], action='append',
                        help='Colon delimited paths to search for splicer files, may be supplied multiple times to create path')
    parser.add_argument('filename', nargs='*',
                        help='Input file to process')

    args = parser.parse_args()

    # check command line options
    if len(args.filename) == 0:
        raise SystemExit("Must give at least one input file")
    if args.outdir and not os.path.isdir(args.outdir):
        raise SystemExit("outdir %s does not exist" % args.outdir)
    if args.outdir_c_fortran and not os.path.isdir(args.outdir_c_fortran):
        raise SystemExit("outdir-fortran %s does not exist" % args.outdir_c_fortran)
    if args.outdir_python and not os.path.isdir(args.outdir_python):
        raise SystemExit("outdir-python %s does not exist" % args.outdir_python)
    if args.outdir_lua and not os.path.isdir(args.outdir_lua):
        raise SystemExit("outdir-luan %s does not exist" % args.outdir_lua)
    if args.logdir and not os.path.isdir(args.logdir):
        raise SystemExit("logdir %s does not exist" % args.logdir)

    # append all paths together
    if args.path:
        search_path = []
        for pth in args.path:
            search_path.extend(pth.split(':'))
    else:
        search_path = ['.']

    basename = os.path.splitext(os.path.basename(args.filename[0]))[0]
    logpath = os.path.join(args.logdir, basename + '.log')
    log = open(logpath, 'w')

    # pass around a configuration object
    config = Config()
    config.c_fortran_dir = args.outdir_c_fortran or args.outdir
    config.python_dir = args.outdir_python or args.outdir
    config.lua_dir = args.outdir_lua or args.outdir
    config.log = log
    config.cfiles = []  # list of C/C++ files created
    config.ffiles = []  # list of Fortran files created

    # accumulated input
    all = dict(
        library = 'default_library',
        cpp_header = '',
        namespace = '',
        )
    splicers = dict(c={}, f={}, py={}, lua={})

    for filename in args.filename:
        root, ext = os.path.splitext(filename)
        if ext == '.yaml':
            log.write("Read yaml %s\n" % os.path.basename(filename))
            fp = open(filename, 'r')
            d = yaml.load(fp.read())
            fp.close()
            all.update(d)
#            util.update(all, d)  # recursive update
        elif ext == '.json':
            raise NotImplemented("Can not deal with json input for now")
        else:
            # process splicer file on command line, search path is not used
            splicer.get_splicer_based_on_suffix(filename, splicers)

#    print(all)


    Schema(all).check_schema()
    VerifyAttrs(all, config).verify_attrs()
    GenFunctions(all, config).gen_library()
    Namify(all, config).name_library()

    if 'splicer' in all:
        # read splicer files defined in input yaml file
        for suffix, names in all['splicer'].items():
            # suffix = 'c', 'f', 'py', 'lua'
            subsplicer = splicers.setdefault(suffix, {})
            for name in names:
                for pth in search_path:
                    fullname = os.path.join(pth, name)
#                    log.write("Try splicer %s\n" % fullname)
                    if os.path.isfile(fullname):
                        break
                else:
                    fullname = None
                if not fullname:
                    raise RuntimeError("File not found: %s" % name)
                log.write("Read splicer %s\n" % name)
                splicer.get_splicers(fullname, subsplicer)

                
    if all['options'].wrap_c:
        wrapc.Wrapc(all, config, splicers['c']).wrap_library()

    if all['options'].wrap_fortran:
        wrapf.Wrapf(all, config, splicers['f']).wrap_library()

    if all['options'].wrap_python:
        wrapp.Wrapp(all, config, splicers['py']).wrap_library()

    if all['options'].wrap_lua:
        wrapl.Wrapl(all, config, splicers['lua']).wrap_library()

    # when dumping json, remove function_index to avoid duplication
    del all['function_index']

    jsonpath = os.path.join(args.logdir, basename + '.json')
    fp = open(jsonpath, 'w')
    json.dump(all, fp, cls=util.ExpandedEncoder, sort_keys=True, indent=4)
    fp.close()

    # Write list of output files.  May be useful for build systems
    if args.cfiles:
        with open(args.cfiles, 'w') as fp:
            if config.cfiles:
                fp.write(' '.join(config.cfiles))
            fp.write('\n')
    if args.ffiles:
        with open(args.ffiles, 'w') as fp:
            if config.ffiles:
                fp.write(' '.join(config.ffiles))
            fp.write('\n')


    log.close()

# This helps when running with a pipe, like CMake's execute_process
# It doesn't fix the error, but it does report a better error message
# http://www.thecodingforums.com/threads/help-with-a-piping-error.749747/
    sys.stdout.flush()

#    sys.stderr.write("Some useful message")  # example error message
    sys.exit(0)  # set status for errors


if __name__ == '__main__':
    main()
