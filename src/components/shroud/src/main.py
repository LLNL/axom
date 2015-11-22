#!/bin/env python3

"""
generate language bindings
"""
from __future__ import print_function

import os
import json
import argparse
import yaml
import sys

import util
import parse_decl
import splicer
import wrapc
import wrapf
import wrapp

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

            library='default_library',
            namespace='',
            cpp_header='',

            F_module_per_class=True,
            F_string_len_trim=True,

            wrap_c       = True,
            wrap_fortran = True,
            wrap_python  = True,

            C_header_filename_library_template = 'wrap{library}.h',
            C_impl_filename_library_template = 'wrap{library}.cpp',

            C_header_filename_class_template = 'wrap{cpp_class}.h',
            C_impl_filename_class_template = 'wrap{cpp_class}.cpp',

            C_name_method_template = '{C_prefix}{lower_class}_{underscore_name}{function_suffix}',
            C_name_function_template = '{C_prefix}{underscore_name}{function_suffix}',

            F_name_impl_method_template = '{lower_class}_{underscore_name}{function_suffix}',
            F_name_impl_function_template ='{underscore_name}{function_suffix}',

            F_name_method_template = '{underscore_name}{function_suffix}',
            F_name_generic_template = '{underscore_name}',

            F_module_name_library_template = '{lower_library}_mod',
            F_impl_filename_library_template = 'wrapf{lower_library}.f',

            F_module_name_class_template = '{lower_class}_mod',
            F_impl_filename_class_template = 'wrapf{cpp_class}.f',

            )
        wrapp.add_templates(def_options)

        if 'options' in node:
            def_options.update(node['options'])
        self.options_stack = [ def_options ]
        node['options'] = def_options

        fmt_library = node['fmt'] = util.Options(None)
        fmt_library.library       = def_options['library']
        fmt_library.lower_library = fmt_library.library.lower()
        fmt_library.upper_library = fmt_library.library.upper()
        fmt_library.function_suffix = ''   # assume no suffix
        fmt_library.overloaded    = False
        fmt_library.C_prefix      = def_options.get('C_prefix', fmt_library.upper_library[:3] + '_')
        fmt_library.rv            = 'rv'  # return value
        util.eval_template2(node, 'C_header_filename', '_library')
        util.eval_template2(node, 'C_impl_filename', '_library')
        self.fmt_stack.append(fmt_library)

        # default some options based on other options
        # All class/methods and functions may go into this file or just functions.
        util.eval_template2(node, 'F_module_name', '_library')
        util.eval_template2(node, 'F_impl_filename', '_library')

        def_types = dict(
            void    = util.Typedef('void',
                c_type   = 'void',
                cpp_type = 'void',
#                fortran = 'subroutine',
                c_fortran = 'type(C_PTR)',
                f_type = 'type(C_PTR)',
                ),
            int    = util.Typedef('int',
                c_type    = 'int',
                cpp_type  = 'int',
                f_kind    = 'C_INT',
                f_cast    = 'int({var}, C_INT)',
                c_fortran = 'integer(C_INT)',
                f_type    = 'integer(C_INT)',
                f_module  = dict(iso_c_binding=['C_INT']),
                PY_format = 'i',
                ),
            long   = util.Typedef('long',
                c_type    = 'long',
                cpp_type  = 'long',
                f_kind    = 'C_LONG',
                f_cast    = 'int({var}, C_LONG)',
                c_fortran = 'integer(C_LONG)',
                f_type    = 'integer(C_LONG)',
                f_module  = dict(iso_c_binding=['C_LONG']),
                PY_format = 'l',
                ),
            size_t   = util.Typedef('size_t',
                c_type    = 'size_t',
                cpp_type  = 'size_t',
                c_header  = 'stdlib.h',
                f_kind    = 'C_SIZE_T',
                f_cast    = 'int({var}, C_SIZE_T)',
                c_fortran = 'integer(C_SIZE_T)',
                f_type    = 'integer(C_SIZE_T)',
                f_module  = dict(iso_c_binding=['C_SIZE_T']),
                PY_ctor   = 'PyInt_FromLong({rv})',
                ),

            float   = util.Typedef('float',
                c_type    = 'float',
                cpp_type  = 'float',
                f_kind    = 'C_FLOAT',
                f_cast    = 'real({var}, C_FLOAT)',
                c_fortran = 'real(C_FLOAT)',
                f_type    = 'real(C_FLOAT)',
                f_module  = dict(iso_c_binding=['C_FLOAT']),
                PY_format = 'f',
                ),
            double   = util.Typedef('double',
                c_type    = 'double',
                cpp_type  = 'double',
                f_kind    = 'C_DOUBLE',
                f_cast    = 'real({var}, C_DOUBLE)',
                c_fortran = 'real(C_DOUBLE)',
                f_type    = 'real(C_DOUBLE)',
                f_module  = dict(iso_c_binding=['C_DOUBLE']),
                PY_format = 'd',
                ),

            bool   = util.Typedef('bool',
                c_type    = 'bool',
                cpp_type  = 'bool',
                f_kind    = 'C_BOOL',
                c_fortran = 'logical(C_BOOL)',

                f_type    = 'logical',
                f_use_tmp  = True,
                f_pre_decl = 'logical(C_BOOL) {tmp_var}',
                f_pre_call = '{tmp_var} = {var}  ! coerce to C_BOOL',

                PY_ctor   = 'PyBool_FromLong({rv})',
#                PY_PyTypeObject = 'PyBool_Type',
#  after parsearg, expectArgs = PyObject_IsTrue(py_expectArgs);
                ),
            # implies null terminated string
            string = util.Typedef('string',
                c_type   = 'char',    # XXX - char *
                cpp_type = 'std::string',
                cpp_to_c = '{var}.c_str()',  # . or ->
                c_fortran  = 'character(kind=C_CHAR)',
                f_type     = 'character(*)',
                fortran_to_c = 'trim({var}) // C_NULL_CHAR',
#                f_module = dict(iso_c_binding = [ 'C_NULL_CHAR' ]),
                f_module = dict(iso_c_binding=None),
                f_return_code = '{F_result} = fstr({F_C_name}({F_arg_c_call_tab}))',
                PY_format = 's',
                PY_ctor = 'PyString_FromString({var})',
                base = 'string',
                ),
            # create std::string from buffer and length
            string_from_buffer = util.Typedef('string_from_buffer',
                c_type   = 'char',    # XXX - char *
                c_argdecl = ['const char *{var}', 'int len_{var}'],
                c_to_cpp = 'std::string({var}, len_{var})',
                cpp_type = 'std::string',
                cpp_to_c = '{var}.c_str()',  # . or ->
                c_fortran  = 'character(kind=C_CHAR)',
                f_c_args   = [ '{var}', 'len_{var}'],
                f_c_argdecl = [ 'type(C_PTR), intent(IN), value :: {var}',
                                'integer(C_INT), intent(IN), value :: len_{var}' ],
                f_type     = 'character(*)',
                fortran_to_c = '{var}, len_trim({var})',
#                f_module = dict(iso_c_binding = [ 'C_NULL_CHAR' ]),
                f_module = dict(iso_c_binding=None),
                f_return_code = '{F_result} = fstr({F_C_name}({F_arg_c_call_tab}))',
                base = 'string',
                ),
            )

        # aliases
        def_types['std::string']     = def_types['string']
        def_types['integer(C_INT)']  = def_types['int']
        def_types['integer(C_LONG)'] = def_types['long']
        def_types['real(C_FLOAT)']   = def_types['float']
        def_types['real(C_DOUBLE)']  = def_types['double']

        # result_as_arg
        tmp = def_types['string'].clone_as('string_result_as_arg')
        tmp.update(dict(
                f_rv_decl     = 'character(*), intent(OUT) :: {result_arg}',
                f_pre_decl    = 'type(C_PTR) :: {F_result}',
                f_return_code = '{F_result} = {F_C_name}({F_arg_c_call_tab})',
                f_post_call   = 'call FccCopyPtr({result_arg}, len({result_arg}), {F_result})',
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

        functions = node.setdefault('functions', [])
        self.check_functions(functions)

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

        options = self.push_options(node)
        fmt_class = self.push_fmt(node)
        fmt_class.cpp_class = name
        fmt_class.lower_class = name.lower()
        fmt_class.upper_class = name.upper()
        if 'C_prefix' in options:
            fmt_class.C_prefix = options.C_prefix

        if options.F_module_per_class:
            util.eval_template2(node, 'F_module_name', '_class')
            util.eval_template2(node, 'F_impl_filename', '_class')

        # Only one file per class for C.
        util.eval_template2(node, 'C_header_filename', '_class')
        util.eval_template2(node, 'C_impl_filename', '_class')

        methods = node.setdefault('methods', [])
        for method in methods:
            if not isinstance(method, dict):
                raise TypeError("classes[n]['methods'] must be a dictionary")
            self.check_function(method)

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
            # YAML turns blanks strings to None
            node['function_suffix'] = ''
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

        fmt_func.method_name =     result['name']
        fmt_func.underscore_name = util.un_camel(result['name'])
        if 'function_suffix' in node:
            fmt_func.function_suffix = node['function_suffix']

        # docs
        self.pop_fmt()
        self.pop_options()

    def check_functions(self, func_list):
        """ check functions which are not in a class.
        """
        if not isinstance(func_list, list):
            raise TypeError("functions must be a list")
        for func in func_list:
            self.check_function(func)


class GenFunctions(object):
    """
    Generate types from class.
    Generate functions based on overload/template/generic/attributes
    """

    def __init__(self, tree, config):
        self.tree = tree    # json tree
        self.config = config
        self.typedef = tree['types']

    def gen_library(self):
        tree = self.tree
        tree['function_index'] = self.function_index = []

        for cls in tree['classes']:
            self.create_class_typedef(cls)

        for cls in tree['classes']:
            cls['methods'] = self.define_function_suffix(cls['methods'])
        tree['functions'] = self.define_function_suffix(tree['functions'])

        for cls in tree['classes']:
            self.check_class_dependencies(cls)

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
                cpp_to_c = 'static_cast<{C_const}%s *>(static_cast<{C_const}void *>({var}))' % cname,
                c_type = cname,
                # opaque pointer -> void pointer -> class instance pointer
                c_to_cpp = 'static_cast<{C_const}%s{ptr}>(static_cast<{C_const}void *>({var}))' % name,
                c_fortran = 'type(C_PTR)',
                f_type = 'type(%s)' % unname,
                fortran_derived = unname,
                fortran_to_c = '{var}%{F_derived_member}',
                # XXX module name may not conflict with type name
                f_module = {fmt_class.F_module_name:[unname]},

                # return from C function
#                f_c_return_decl = 'type(CPTR)' % unname,
                f_return_code = '{F_result}%{F_derived_member} = {F_C_name}({F_arg_c_call_tab})',

                # allow forward declarations to avoid recursive headers
                forward = name,
                base = 'wrapped',
                )

        typedef = self.typedef[name]
        fmt_class.C_type_name = typedef.c_type

    def append_function_index(self, node):
        """append to function_index, set index into node.
        """
        ilist = self.function_index
        node['function_index'] = len(ilist)
#        node['fmt'].function_index = str(len(ilist)) # debugging
        ilist.append(node)

    def define_function_suffix(self, functions):
        """ look for functions with the same name
        Return a new list if overloaded function inserted.
        """
        overloaded_functions = {}
        for function in functions:
            self.append_function_index(function)
            overloaded_functions.setdefault(function['result']['name'], []).append(function)

        # look for function overload and compute function_suffix
        for mname, overloads in overloaded_functions.items():
            if len(overloads) > 1:
                for i, function in enumerate(overloads):
#                    method['fmt'].overloaded = True
                    if 'function_suffix' not in function:
                        function['fmt'].function_suffix =  '_%d' % i

        # Create additional functions needed for wrapping
        ordered_functions = []
        for method in functions:
            ordered_functions.append(method)
            if 'cpp_template' in method:
                self.template_function(method, ordered_functions)
            if 'fortran_generic' in method:
                self.generic_function(method, ordered_functions)
            self.string_to_buffer_and_len(method, ordered_functions)
        return ordered_functions

    def template_function(self, node, ordered_functions):
        """ Create overloaded functions for each templated method.
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

                new['generated'] = 'cpp_template'
                fmt = new['fmt']
                fmt.function_suffix = '_' + type
                del new['cpp_template']
                options = new['options']
                options.wrap_c = True
                options.wrap_fortran = True
                options.wrap_python = False
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

                new['generated'] = 'fortran_generic'
                fmt = new['fmt']
                fmt.function_suffix = '_' + type
                fmt.PTR_F_C_index = node['function_index']
                del new['fortran_generic']
                options = new['options']
                options.wrap_c = False
                options.wrap_fortran = True
                options.wrap_python = False
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

    def string_to_buffer_and_len(self, node, ordered_functions):
        """ Check if function has any string arguments and will be wrapped by Fortran.
        If so then create a new C function that will convert string arguments into 
        a buffer and length.
        """
        options = node['options']
        if options.wrap_fortran is False:
            return
        if options.F_string_len_trim is False:
            return

        has_strings = False
        for arg in node['args']:
            argtype = arg['type']
            if self.typedef[argtype].base == 'string':
                has_strings = True
                break
        if has_strings is False:
            return

        new = util.copy_function_node(node)
        ordered_functions.append(new)
        self.append_function_index(new)

        new['generated'] = 'string_to_buffer_and_len'
        fmt = new['fmt']
        fmt.function_suffix = '_bufferify'

        options = new['options']
        options.wrap_c = True
        options.wrap_fortran = False
        options.wrap_python = False

        options = node['options']
        #        options.wrap_fortran = False
#        # Current Fortran function should use this new C function
        node['fmt'].PTR_F_C_index = new['function_index']

        newargs = []
        for arg in new['args']:
            argtype = arg['type']
            if self.typedef[argtype].base == 'string':
                # Add len_trim attribute
                arg['attrs']['len_trim'] = True

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

class VerifyAttrs(object):
    """
    This must be called after GenFunctions has generated typedefs
    for classes.
    """
    def __init__(self, tree, config):
        self.tree = tree    # json tree
        self.config = config

    def verify_attrs(self):
        tree = self.tree
        self.typedef = tree['types']

        for cls in tree['classes']:
            for func in cls['methods']:
                self.check_arg_attrs(func)

        for func in tree['functions']:
            self.check_arg_attrs(func)

    def check_arg_attrs(self, node):
        """Regularize attributes
        intent: lower case, no parens, must be in, out, or inout
        value: if pointer, default to False (pass-by-reference;
               else True (pass-by-value).
        """
        options = node['options']
        if not options.wrap_fortran and not options.wrap_c:
            return

        for arg in node['args']:
            typedef = self.typedef.get(arg['type'], None)
            if typedef is None:
                raise RuntimeError("No such type %s" % arg['type'])

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
                else:
                    attrs['intent'] = 'inout'  # Fortran default
                    attrs['intent'] = 'in' # must coordinate with VALUE
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

#        if typedef.base == 'string':


class Namify(object):
    """Compute names of functions in library.
    Need to compute F_name and C_F_name since they interact.
    Compute all C names first, then Fortran.
    A Fortran function may call a generated C function via PTR_F_C_index
    Also compute number which may be controlled by options.    

    C_name - Name of C function
    F_C_name - Fortran function for C interface
    F_name - Name of Fortran function
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
        
        if cls:
            util.eval_template2(node, 'C_name', '_method')
        else:
            util.eval_template2(node, 'C_name', '_function')

        if 'F_C_name' in node:
            fmt_func.F_C_name = node['F_C_name']
        else:
            fmt_func.F_C_name = fmt_func.C_name.lower()

        if 'C_this' in options:
            fmt_func.C_this = options.C_this

    def name_function_fortran(self, cls, node):
        """ Must process C functions for to generate their names.
        """
        options = node['options']
        if not options.wrap_fortran:
            return
        fmt_func = node['fmt']

        if cls:
            util.eval_template2(node, 'F_name_impl', '_method')
        else:
            util.eval_template2(node, 'F_name_impl', '_function')

        util.eval_template(options, fmt_func, 'F_name_method')
        util.eval_template2(node, 'F_name_generic')

        if 'F_this' in options:
            fmt_func.F_this = options.F_this
        if 'F_result' in options:
            fmt_func.F_result = options.F_result


if __name__ == '__main__':

    appname = 'modulator3'
    appver = '0.1'

    parser = argparse.ArgumentParser(prog=appname)
    parser.add_argument('--version', action='version', version=appver)
    parser.add_argument('--outdir', default='',
                        help='directory for output files')
    parser.add_argument('--logdir', default='',
                        help='directory for log files')
    parser.add_argument('--path', default=[], action='append',
                        help='colon delimited paths to search for splicer files, may be supplied multiple times to create path')
    parser.add_argument('filename', nargs='*',
                        help='Input file to process')

    args = parser.parse_args()

    # check command line options
    if len(args.filename) == 0:
        raise SystemExit("Must give at least one input file")
    if args.outdir and not os.path.isdir(args.outdir):
        raise SystemExit("outdir %s does not exist" % args.outdir)
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
    config.binary_dir = args.outdir
    config.log = log

    # accumulated input
    all = {}
    splicers = dict(c={}, f={}, py={})

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
    GenFunctions(all, config).gen_library()
    VerifyAttrs(all, config).verify_attrs()
    Namify(all, config).name_library()

    if 'splicer' in all:
        # read splicer files defined in input yaml file
        for suffix, names in all['splicer'].items():
            # suffix = 'c', 'f', 'py'
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


    wrapc.Wrapc(all, config, splicers['c']).wrap_library()

    wrapf.Wrapf(all, config, splicers['f']).wrap_library()

    wrapp.Wrapp(all, config, splicers['py']).wrap_library()

    # when dumping json, remove function_index to avoid duplication
    del all['function_index']

    jsonpath = os.path.join(args.logdir, basename + '.json')
    fp = open(jsonpath, 'w')
    json.dump(all, fp, cls=util.ExpandedEncoder, sort_keys=True, indent=4)
    fp.close()

    log.close()

# This helps when running with a pipe, like CMake's execute_process
# It doesn't fix the error, but it does report a better error message
# http://www.thecodingforums.com/threads/help-with-a-piping-error.749747/
    sys.stdout.flush()

#    sys.stderr.write("Some useful message")  # example error message
    sys.exit(0)  # set status for errors

