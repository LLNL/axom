#!/bin/env python3
"""
Generate Fortran bindings for C++ code.

module {F_module_name}

  type {F_derived_name}
    type(C_PTR) {F_derived_member}
  contains
    procedure :: {F_name_method} => {F_name_impl}
    generic :: {F_name_generic} => {F_name_method}, ...
  end type {F_derived_name}

  ! interface for C functions
  interface
    {F_C_pure_clause}{F_C_subprogram} {F_C_name}({F_C_arguments}) &
        {F_C_result_clause} &
        bind(C, name="{C_name}")
      {arg_c_decl}
    end {F_C_subprogram} {F_C_name}
  end interface

  interface {F_name_generic}
    module procedure {F_name_impl}
  end interface {F_name_generic}

contains

  {F_pure_clause} {F_subprogram} {F_name_impl}({F_arguments}){F_result_clause}
      {F_C_name}({F_arg_c_call_tab})
     {F_code}
  end {F_subprogram} {F_name_impl}

end module {F_module_name}
----------
"""
from __future__ import print_function
from __future__ import absolute_import

import os

from . import whelpers
from . import util
from .util import wformat, append_format


class Wrapf(util.WrapperMixin):
    """Generate Fortran bindings.
    """

    def __init__(self, tree, config, splicers):
        self.tree = tree    # json tree
        self.patterns = tree['patterns']
        self.config = config
        self.log = config.log
        self.typedef = tree['types']
        self._init_splicer(splicers)
        self.comment = '!'
        self.doxygen_begin = '!>'
        self.doxygen_cont = '!!'
        self.doxygen_end = '!<'

    def _begin_output_file(self):
        """Start a new class for output"""
        self.use_stmts = []
        self.f_type_decl = []
        self.c_interface = []
        self.generic_interface = []
        self.impl = []          # implementation, after contains
        self.operator_impl = []
        self.operator_map = {}  # list of function names by operator
        # {'.eq.': [ 'abc', 'def'] }
        self.c_interface.append('interface')
        self.c_interface.append(1)
        self.f_helper = {}

    def _end_output_file(self):
        self.c_interface.append(-1)
        self.c_interface.append('end interface')

    def _begin_class(self):
        self.f_type_generic = {}  # look for generic methods
        self.type_bound_part = []

    def _c_type(self, arg):
        """
        Return the Fortran type, and dimension
        pass-by-value default

        attributes:
        ptr - True = pass-by-reference
        reference - True = pass-by-reference
        dimension - True = '(*)'

        """
        t = []
        typedef = self.typedef[arg['type']]
        attrs = arg['attrs']
        intent = attrs.get('intent', None)

        typ = typedef.c_fortran
        if typ is None:
            raise RuntimeError("Type {} has no value for c_fortran".format(arg['type']))
        t.append(typ)
        if attrs.get('value', False):
            t.append('value')
        if intent:
            t.append('intent(%s)' % intent.upper())
        if typedef.base == 'string':
            dimension = '(*)'  # is array
        else:
            # XXX should C always have dimensions of '(*)'?
            dimension = attrs.get('dimension', '')
        return (', '.join(t), dimension)

    def _c_decl(self, arg, name=None):
        """
        Return the Fortran declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ, dimension = self._c_type(arg)
        rv = typ + ' :: ' + (name or arg['name']) + dimension
        return rv

    def _f_type(self, arg, default=None, local=False):
        """
        Return the Fortran type, and array attribute
        pass-by-value default

        attributes:
        ptr - True = pass-by-reference
        reference - True = pass-by-reference

        If local==True, this is a local variable, skip attributes
          OPTIONAL, VALUE, and INTENT
        """
        t = []
        typedef = self.typedef[arg['type']]
        attrs = arg['attrs']
        intent = attrs.get('intent', None)

        typ = typedef.f_type
        t.append(typ)
        if not local:  # must be dummy argument
            if attrs.get('value', False):
                t.append('value')
            if intent:
                t.append('intent(%s)' % intent.upper())
            # --            if default is None:
            # --                default = attrs.get('default', '')
            # --            if default != '':
            # --                t.append('optional')
        dimension = attrs.get('dimension', '')
        return (', '.join(t), dimension)

    def _f_decl(self, arg, name=None, default=None, local=False):
        """
        Return the Fortran declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ, dimension = self._f_type(arg, default=default, local=local)
        rv = typ + ' :: ' + (name or arg['name']) + dimension
        return rv

    def wrap_library(self):
        options = self.tree['options']
        fmt_library = self.tree['fmt']
        fmt_library.F_result_clause = ''
        fmt_library.F_pure_clause = ''
        fmt_library.F_C_result_clause = ''
        fmt_library.F_C_pure_clause = ''

        self._begin_output_file()
        self._push_splicer('class')
        for node in self.tree['classes']:
            self._begin_class()

            name = node['name']
            # how to decide module name, module per class
#            module_name = node['options'].setdefault('module_name', name.lower())
            self.wrap_class(node)
            if options.F_module_per_class:
                self._pop_splicer('class')
                self._end_output_file()
                self.write_module(self.tree, node)
                self._begin_output_file()
                self._push_splicer('class')
        self._pop_splicer('class')

        if self.tree['functions']:
            self._begin_class()  # clear out old class info
            if options.F_module_per_class:
                self._begin_output_file()
            self.tree['F_module_dependencies'] = []
            for node in self.tree['functions']:
                self.wrap_function(None, node)
            self.c_interface.append('')
            self._create_splicer('additional_interfaces', self.c_interface)
            self.impl.append('')
            self._create_splicer('additional_functions', self.impl)

            # Look for generics
            # splicer to extend generic
            self._push_splicer('generic')
            iface = self.generic_interface
            for key in sorted(self.f_type_generic.keys()):
                generics = self.f_type_generic[key]
                if len(generics) > 1:
                    self._push_splicer(key)
                    iface.append('')
                    iface.append('interface ' + key)
                    iface.append(1)
                    for genname in generics:
                        iface.append('module procedure ' + genname)
                    iface.append(-1)
                    iface.append('end interface ' + key)
                    self._pop_splicer(key)
            self._pop_splicer('generic')

            if options.F_module_per_class:
                # library module
                self._end_output_file()
                self._create_splicer('module_use', self.use_stmts)
                self.write_module(self.tree, None)

        if not options.F_module_per_class:
            # put all functions and classes into one module
            self.tree['F_module_dependencies'] = []
            self._end_output_file()
            self.write_module(self.tree, None)

        self.write_c_helper()

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        unname = util.un_camel(name)
        typedef = self.typedef[name]

        fmt_class = node['fmt']

        fmt_class.F_derived_name = typedef.f_derived_type

        # wrap methods
        self._push_splicer(fmt_class.cpp_class)
        self._create_splicer('module_use', self.use_stmts)
        self._push_splicer('method')
        for method in node['methods']:
            self.wrap_function(node, method)
        self._pop_splicer('method')
        self.write_object_get_set(node, fmt_class)
        self.impl.append('')
        self._create_splicer('additional_functions', self.impl)
        self._pop_splicer(fmt_class.cpp_class)

        # type declaration
        self.f_type_decl.append('')
        self._push_splicer(fmt_class.cpp_class)
        self._create_splicer('module_top', self.f_type_decl)
        self.f_type_decl.extend([
                '',
                wformat('type {F_derived_name}', fmt_class),
                1,
                wformat('type(C_PTR), private :: {F_derived_member}', fmt_class),
                ])
        self._create_splicer('component_part', self.f_type_decl)
        self.f_type_decl.extend([
                -1, 'contains', 1,
                ])
        self.f_type_decl.extend(self.type_bound_part)

        # Look for generics
        # splicer to extend generic
        self._push_splicer('generic')
        for key in sorted(self.f_type_generic.keys()):
            methods = self.f_type_generic[key]
            if len(methods) > 1:
                self.f_type_decl.append('generic :: %s => &' % key)
                self.f_type_decl.append(1)
                self._create_splicer(key, self.f_type_decl)
                for genname in methods[:-1]:
                    self.f_type_decl.append(genname + ',  &')
                self.f_type_decl.append(methods[-1])
                self.f_type_decl.append(-1)
        self._pop_splicer('generic')

        self._create_splicer('type_bound_procedure_part', self.f_type_decl)
        self.f_type_decl.extend([
                 -1,
                 wformat('end type {F_derived_name}', fmt_class),
                 ])

        self.c_interface.append('')
        self._create_splicer('additional_interfaces', self.c_interface)

        self._pop_splicer(fmt_class.cpp_class)

        # overload operators
        self.overload_compare(
            fmt_class, '.eq.', fmt_class.class_lower + '_eq',
            wformat('c_associated(a%{F_derived_member}, b%{F_derived_member})',
                    fmt_class))
#        self.overload_compare(fmt_class, '==', fmt_class.class_lower + '_eq', None)
        self.overload_compare(
            fmt_class, '.ne.', fmt_class.class_lower + '_ne',
            wformat(
                '.not. c_associated'
                '(a%{F_derived_member}, b%{F_derived_member})',
                fmt_class))
#        self.overload_compare(fmt_class, '/=', fmt_class.class_lower + '_ne', None)

    def write_object_get_set(self, node, fmt_class):
        """Write get and set methods for instance pointer.

        node = class dictionary
        """
        options = node['options']
        impl = self.impl
        fmt = util.Options(fmt_class)

        fmt.F_instance_ptr = wformat('{F_this}%{F_derived_member}', fmt)

        # get
        fmt.underscore_name = options['F_name_instance_get']
        if fmt.underscore_name:
            fmt.underscore_name = options['F_name_instance_get']
            fmt.F_name_method = wformat(options['F_name_method_template'], fmt)
            fmt.F_name_impl = wformat(options['F_name_impl_template'], fmt)

            self.type_bound_part.append('procedure :: %s => %s' % (
                    fmt.F_name_method, fmt.F_name_impl))

            impl.append('')
            append_format(
                impl, 'function {F_name_impl}({F_this}) '
                'result ({F_derived_member})', fmt)
            impl.append(1)
            impl.append('use iso_c_binding, only: C_PTR')
            impl.append('implicit none')
            append_format(
                impl, 'class({F_derived_name}), intent(IN) :: {F_this}', fmt)
            append_format(impl, 'type(C_PTR) :: {F_derived_member}', fmt)
            append_format(impl, '{F_derived_member} = {F_instance_ptr}', fmt)
            impl.append(-1)
            append_format(impl, 'end function {F_name_impl}', fmt)

        # set
        fmt.underscore_name = options['F_name_instance_set']
        if fmt.underscore_name:
            fmt.F_name_method = wformat(options['F_name_method_template'], fmt)
            fmt.F_name_impl = wformat(options['F_name_impl_template'], fmt)

            self.type_bound_part.append('procedure :: %s => %s' % (
                    fmt.F_name_method, fmt.F_name_impl))

            impl.append('')
            append_format(
                impl, 'subroutine {F_name_impl}'
                '({F_this}, {F_derived_member})', fmt)
            impl.append(1)
            impl.append('use iso_c_binding, only: C_PTR')
            impl.append('implicit none')
            append_format(
                impl, 'class({F_derived_name}), intent(INOUT) :: {F_this}',
                fmt)
            append_format(
                impl, 'type(C_PTR), intent(IN) :: {F_derived_member}', fmt)
            append_format(impl, '{F_instance_ptr} = {F_derived_member}', fmt)
            impl.append(-1)
            append_format(impl, 'end subroutine {F_name_impl}', fmt)

        # associated
        fmt.underscore_name = options['F_name_associated']
        if fmt.underscore_name:
            fmt.F_name_method = wformat(options['F_name_method_template'], fmt)
            fmt.F_name_impl = wformat(options['F_name_impl_template'], fmt)

            self.type_bound_part.append('procedure :: %s => %s' % (
                    fmt.F_name_method, fmt.F_name_impl))

            impl.append('')
            append_format(
                impl, 'function {F_name_impl}({F_this}) result (rv)', fmt)
            impl.append(1)
            impl.append('use iso_c_binding, only: c_associated')
            impl.append('implicit none')
            append_format(
                impl, 'class({F_derived_name}), intent(IN) :: {F_this}', fmt)
            impl.append('logical rv')
            append_format(impl, 'rv = c_associated({F_instance_ptr})', fmt)
            impl.append(-1)
            append_format(impl, 'end function {F_name_impl}', fmt)

    def overload_compare(self, fmt_class, operator, procedure, predicate):
        """ Overload .eq. and .eq.
        """
        fmt = util.Options(fmt_class)
        fmt.procedure = procedure
        fmt.predicate = predicate

        ops = self.operator_map.setdefault(operator, [])
        ops.append(procedure)

        if predicate is None:
            # .eq. and == use same function
            return

        operator = self.operator_impl
        operator.append('')
        append_format(operator, 'function {procedure}(a,b) result (rv)', fmt)
        operator.append(1)
        operator.append('use iso_c_binding, only: c_associated')
        operator.append('implicit none')
        append_format(operator,
                      'type({F_derived_name}), intent(IN) ::a,b', fmt)
        operator.append('logical :: rv')
        append_format(operator, 'if ({predicate}) then', fmt)
        operator.append(1)
        operator.append('rv = .true.')
        operator.append(-1)
        operator.append('else')
        operator.append(1)
        operator.append('rv = .false.')
        operator.append(-1)
        operator.append('endif')
        operator.append(-1)
        append_format(operator, 'end function {procedure}', fmt)

    def wrap_function(self, cls, node):
        """
        cls  - class node or None for functions
        node - function/method node

        Wrapping involves both a C interface and a Fortran wrapper.
        For some generic functions there may be single C method with
        multiple Fortran wrappers.

        """
        if cls:
            cls_function = 'method'
        else:
            cls_function = 'function'

        options = node['options']
        wrap = []
        if options.wrap_c:
            wrap.append('C-interface')
        if options.wrap_fortran:
            wrap.append('Fortran')
        if not wrap:
            return

        self.log.write(', '.join(wrap))
        self.log.write(" {0} {1[_decl]}\n".format(cls_function, node))

        # Create fortran wrappers first.
        # If no real work to do, call the C function directly.
        if options.wrap_fortran:
            self.wrap_function_impl(cls, node)
        if options.wrap_c:
            self.wrap_function_interface(cls, node)

    def update_f_module(self, modules, f_module):
        """aggragate the information from f_module into modules.
        modules is a dictionary of dictionaries:
            modules['iso_c_bindings']['C_INT'] = True
        f_module is a dictionary of lists:
            dict(iso_c_binding=['C_INT'])
        """
        if f_module is not None:
            for mname, only in f_module.items():
                module = modules.setdefault(mname, {})
                if only:  # Empty list means no ONLY clause
                    for oname in only:
                        module[oname] = True

    def set_f_module(self, modules, mname, *only):
        """Add a module to modules.
        """
        module = modules.setdefault(mname, {})
        if only:  # Empty list means no ONLY clause
            for oname in only:
                module[oname] = True

    def sort_module_info(self, modules, module_name):
        """Return USE statements based on modules.
        Skip module_name.  Usually the current module.
        """
        arg_f_use = []
        for mname in sorted(modules):
            if mname == module_name:
                continue
            only = modules[mname]
            if only:
                arg_f_use.append('use %s, only : %s' % (
                        mname, ', '.join(only)))
            else:
                arg_f_use.append('use %s' % mname)
        return arg_f_use

    def wrap_function_interface(self, cls, node):
        """
        Write Fortran interface for C function
        cls  - class node or None for functions
        node - function/method node

        Wrapping involves both a C interface and a Fortran wrapper.
        For some generic functions there may be single C method with
        multiple Fortran wrappers.
        """
        options = node['options']
        fmt_func = node['fmt']
        fmt = util.Options(fmt_func)

        func_is_const = node['attrs'].get('const', False)

        result = node['result']
        result_type = result['type']
        subprogram = node['_subprogram']

        if node.get('return_this', False):
            result_type = 'void'
            subprogram = 'subroutine'

        result_typedef = self.typedef[result_type]
        is_ctor = node['attrs'].get('constructor', False)
        is_const = node['attrs'].get('const', False)
        is_pure = node['attrs'].get('pure', False)

        arg_c_names = []  # argument names for functions
        arg_c_decl = []   # declaraion of argument names

        # find subprogram type
        # compute first to get order of arguments correct.
        # Add
        if subprogram == 'subroutine':
            fmt.F_C_subprogram = 'subroutine'
        else:
            fmt.F_C_subprogram = 'function'
            fmt.F_C_result_clause = 'result(%s)' % fmt.F_result

        if cls:
            # Add 'this' argument
            if not is_ctor:
                arg_c_names.append(fmt.C_this)
                arg_c_decl.append(
                    'type(C_PTR), value, intent(IN) :: ' + fmt.C_this)

        args_all_in = True   # assume all arguments are intent(in)
        for arg in node['args']:
            # default argument's intent
            # XXX look at const, ptr
            arg_typedef = self.typedef[arg['type']]
            fmt.c_var = arg['name']
            attrs = arg['attrs']

            if attrs.get('intent', 'inout') != 'in':
                args_all_in = False

            # argument names
            if arg_typedef.f_c_args:
                for argname in arg_typedef.f_c_args:
                    arg_c_names.append(argname)
            else:
                arg_c_names.append(arg['name'])

            # argument declarations
            if arg_typedef.f_c_argdecl:
                for argdecl in arg_typedef.f_c_argdecl:
                    append_format(arg_c_decl, argdecl, fmt)
            else:
                arg_c_decl.append(self._c_decl(arg))

            len_trim = attrs.get('len_trim', None)
            if len_trim:
                arg_c_names.append(len_trim)
                arg_c_decl.append(
                    'integer(C_INT), value, intent(IN) :: %s' % len_trim)
            len_arg = attrs.get('len', None)
            if len_arg:
                arg_c_names.append(len_arg)
                arg_c_decl.append(
                    'integer(C_INT), value, intent(IN) :: %s' % len_arg)

        if (subprogram == 'function' and
                (is_pure or (func_is_const and args_all_in))):
            fmt.F_C_pure_clause = 'pure '

        fmt.F_C_arguments = options.get(
            'F_C_arguments', ', '.join(arg_c_names))

        if fmt.F_C_subprogram == 'function':
            if result_typedef.base == 'string':
                arg_c_decl.append('type(C_PTR) %s' % fmt.F_result)
            else:
                # XXX - make sure ptr is set to avoid VALUE
                arg_dict = dict(name=fmt.F_result,
                                type=result_type,
                                attrs=dict(ptr=True))
                arg_c_decl.append(self._c_decl(arg_dict))

        c_interface = self.c_interface
        c_interface.append('')
        c_interface.append(wformat(
            '{F_C_pure_clause}{F_C_subprogram} {F_C_name}'
            '({F_C_arguments}) &',
            fmt))
        c_interface.append(2)  # extra indent for continued line
        if fmt.F_C_result_clause:
            c_interface.append(wformat(
                    '{F_C_result_clause} &',
                    fmt))
        c_interface.append(wformat(
                'bind(C, name="{C_name}")',
                fmt))
        c_interface.append(-1)
        c_interface.append('use iso_c_binding')
        c_interface.append('implicit none')
        c_interface.extend(arg_c_decl)
        c_interface.append(-1)
        c_interface.append(wformat('end {F_C_subprogram} {F_C_name}', fmt))

    def wrap_function_impl(self, cls, node):
        """
        Wrap implementation of Fortran function
        """
        options = node['options']
        fmt_func = node['fmt']

        # Assume that the C function can be called directly.
        # If the wrapper does any work, then set need_wraper to True
        need_wrapper = options['F_force_wrapper']

        # Look for C routine to wrap
        # Usually the same node unless it is a generic function
        C_node = node
        generated = []
        if '_generated' in C_node:
            generated.append(C_node['_generated'])
        while '_PTR_F_C_index' in C_node:
            C_node = self.tree['function_index'][C_node['_PTR_F_C_index']]
            if '_generated' in C_node:
                generated.append(C_node['_generated'])
#  #This is no longer true with the result as an argument
#        if len(node['args']) != len(C_node['args']):
#            raise RuntimeError("Argument mismatch between Fortran and C functions")

        fmt_func.F_C_call = C_node['fmt'].F_C_name

        func_is_const = node['attrs'].get('const', False)

        # Fortran return type
        result = node['result']
        result_type = result['type']
        subprogram = node['_subprogram']
        c_subprogram = C_node['_subprogram']

        if node.get('return_this', False):
            result_type = 'void'
            subprogram = 'subroutine'
            c_subprogram = 'subroutine'

        result_typedef = self.typedef[result_type]
        is_ctor = node['attrs'].get('constructor', False)
        is_dtor = node['attrs'].get('destructor', False)
        is_const = result['attrs'].get('const', False)

        # this catches stuff like a bool to logical conversion which
        # requires the wrapper
        if result_typedef.f_statements.get('result', {}) \
                                      .get('need_wrapper', False):
            need_wrapper = True

        fmt_func.F_instance_ptr = wformat('{F_this}%{F_derived_member}', fmt_func)

        arg_c_call = []      # arguments to C function

        arg_f_names = []
        arg_f_decl = []
        modules = {}   # indexed as [module][variable]

        if subprogram == 'function':
            fmt_func.F_result_clause = ' result(%s)' % fmt_func.F_result
        fmt_func.F_subprogram = subprogram

        if cls:
            need_wrapper = True
            # Add 'this' argument
            if not is_ctor:
                arg_c_call.append(fmt_func.F_instance_ptr)
                arg_f_names.append(fmt_func.F_this)
                arg_f_decl.append(wformat(
                        'class({F_derived_name}) :: {F_this}',
                        fmt_func))

        #
        # Fortran and C arguments may have different types (fortran generic)
        #
        # f_var - argument to Fortran function (wrapper function)
        # c_var - argument to C function (wrapped function)
        #
        # May be one more argument to C function than Fortran function
        # (the result)
        #
        pre_call = []
        f_args = node['args']
        f_index = -1       # index into f_args
        for c_index, c_arg in enumerate(C_node['args']):
            fmt_arg = c_arg.setdefault('fmtf', util.Options(fmt_func))
            fmt_arg.f_var = c_arg['name']
            fmt_arg.c_var = fmt_arg.f_var

            f_arg = True   # assume C and Fortran arguments match
            c_attrs = c_arg['attrs']
            if c_attrs.get('_is_result', False):
                result_as_arg = options.get('F_string_result_as_arg', '')
                if not result_as_arg:
                    # passing Fortran function result variable down to C
                    f_arg = False
                    fmt_arg.c_var = fmt_func.F_result
                    fmt_arg.f_var = fmt_func.F_result

            if f_arg:
                f_index += 1
                f_arg = f_args[f_index]
                arg_f_names.append(fmt_arg.f_var)
                arg_f_decl.append(self._f_decl(f_arg))

                arg_type = f_arg['type']
                arg_typedef = self.typedef[arg_type]

                f_statements = arg_typedef.f_statements
                for intent in ['intent_in']:
                    cmd_list = f_statements.get(intent, {}).get('declare', [])
                    if cmd_list:
                        need_wrapper = True
                        fmt_arg.c_var = 'SH_' + fmt_arg.f_var
                        for cmd in cmd_list:
                            append_format(arg_f_decl, cmd, fmt_arg)

                    cmd_list = f_statements.get(intent, {}).get('pre_call', [])
                    if cmd_list:
                        need_wrapper = True
                        for cmd in cmd_list:
                            append_format(pre_call, cmd, fmt_arg)

                    self.update_f_module(modules, arg_typedef.f_module)

            # Now C function arguments
            # May have different types, like generic
            # or different attributes, like adding +len to string args
            arg_typedef = self.typedef[c_arg['type']]

            # Attributes   None=skip, True=use default, else use value
            len_trim = c_arg['attrs'].get('len_trim', None)
            len_arg = c_arg['attrs'].get('len', None)
            if arg_typedef.f_args:
                need_wrapper = True
                append_format(arg_c_call, arg_typedef.f_args, fmt_arg)
            elif arg_typedef.f_to_c:
                need_wrapper = True
                append_format(arg_c_call, arg_typedef.f_to_c, fmt_arg)
            elif f_arg and c_arg['type'] != f_arg['type']:
                need_wrapper = True
                append_format(arg_c_call, arg_typedef.f_cast, fmt_arg)
                self.update_f_module(modules, arg_typedef.f_module)
            else:
                append_format(arg_c_call, '{c_var}', fmt_arg)

            len_trim = c_arg['attrs'].get('len_trim', None)
            if len_trim:
                need_wrapper = True
                append_format(arg_c_call, 'len_trim({f_var}, kind=C_INT)', fmt_arg)
                self.set_f_module(modules, 'iso_c_binding', 'C_INT')
            len_arg = c_arg['attrs'].get('len', None)
            if len_arg:
                need_wrapper = True
                append_format(arg_c_call, 'len({f_var}, kind=C_INT)', fmt_arg)
                self.set_f_module(modules, 'iso_c_binding', 'C_INT')

        fmt_func.F_arg_c_call = ', '.join(arg_c_call)
        # use tabs to insert continuations
        fmt_func.F_arg_c_call_tab = '\t' + '\t'.join(arg_c_call)
        fmt_func.F_arguments = options.get('F_arguments', ', '.join(arg_f_names))

        # declare function return value after arguments
        # since arguments may be used to compute return value
        # (for example, string lengths)
        if subprogram == 'function':
            # if func_is_const:
            #     fmt_func.F_pure_clause = 'pure '
            if result_typedef.base == 'string':
                # special case returning a string
                rvlen = result['attrs'].get('len', None)
                if rvlen is None:
                    rvlen = wformat(
                        'strlen_ptr({F_C_call}({F_arg_c_call_tab}))',
                        fmt_func)
                else:
                    rvlen = str(rvlen)  # convert integers
                fmt_func.rvlen = wformat(rvlen, fmt_func)
                line1 = wformat(
                    'character(kind=C_CHAR, len={rvlen}) :: {F_result}',
                    fmt_func)
                self.append_method_arguments(arg_f_decl, line1)
                self.set_f_module(modules, 'iso_c_binding', 'C_CHAR')
            else:
                arg_f_decl.append(self._f_decl(result, name=fmt_func.F_result))
            self.update_f_module(modules, result_typedef.f_module)

        if not is_ctor:
            # Add method to derived type
            if '_overloaded' in node:
                need_wrapper = True
            if '_CPP_return_templated' not in node:
                # if return type is templated in C++,
                # then do not set up generic since only the
                # return type may be different (ex. getValue<T>())
                if cls:
                    gname = fmt_func.F_name_method
                else:
                    gname = fmt_func.F_name_impl
                self.f_type_generic.setdefault(
                    fmt_func.F_name_generic, []).append(gname)
            self.type_bound_part.append('procedure :: %s => %s' % (
                    fmt_func.F_name_method, fmt_func.F_name_impl))

        # body of function
        # XXX sname = fmt_func.F_name_impl
        sname = fmt_func.F_name_method
        splicer_code = self.splicer_stack[-1].get(sname, None)
        if 'F_code' in node:
            need_wrapper = True
            F_code = [wformat(node['F_code'], fmt_func)]
        elif splicer_code:
            need_wrapper = True
            F_code = splicer_code
        else:
            F_code = []
            if is_ctor:
                line1 = wformat(
                    '{F_result}%{F_derived_member} = '
                    '{F_C_call}({F_arg_c_call_tab})', fmt_func)
                self.append_method_arguments(F_code, line1)
            elif c_subprogram == 'function':
                f_return_code = result_typedef.f_return_code
                if f_return_code is None:
                    f_return_code = (
                        '{F_result} = {F_C_call}({F_arg_c_call_tab})')
                else:
                    self.f_helper.update(
                        result_typedef.f_helper.get('f_return_code', {}))
                    need_wrapper = True
                line1 = wformat(f_return_code, fmt_func)
                self.append_method_arguments(F_code, line1)
            else:
                line1 = wformat('call {F_C_call}({F_arg_c_call_tab})', fmt_func)
                self.append_method_arguments(F_code, line1)

#            if result_typedef.f_post_call:
#                need_wrapper = True
#                # adjust return value or cleanup
#                append_format(F_code, result_typedef.f_post_call, fmt_func)
            if is_dtor:
                F_code.append(wformat(
                    '{F_this}%{F_derived_member} = C_NULL_PTR', fmt_func))
                self.set_f_module(modules, 'iso_c_binding', 'C_NULL_PTR')

        arg_f_use = self.sort_module_info(modules, fmt_func.F_module_name)

        if need_wrapper:
            impl = self.impl
            impl.append('')
            if options.debug:
                impl.append('! %s' % node['_decl'])
                if generated:
                    impl.append('! %s' % ' - '.join(generated))
                impl.append('! function_index=%d' % node['_function_index'])
                if options.doxygen and 'doxygen' in node:
                    self.write_doxygen(impl, node['doxygen'])
            impl.append(wformat(
                '{F_subprogram} {F_name_impl}'
                '({F_arguments}){F_result_clause}', fmt_func))
            impl.append(1)
            impl.extend(arg_f_use)
            impl.append('implicit none')
            impl.extend(arg_f_decl)
            impl.extend(pre_call)
            self._create_splicer(sname, impl, F_code)
            impl.append(-1)
            impl.append(wformat('end {F_subprogram} {F_name_impl}', fmt_func))
        else:
            fmt_func.F_C_name = fmt_func.F_name_impl

    def append_method_arguments(self, F_code, line1):
        """Append each argment in arg_c_call as a line in the function.
        Must account for continuations
        Replace tabs in line1 with continuations.
        """
        # part[0] = beginning
        # part[1] = first argument
        parts = line1.split('\t')

        if len(parts) < 3:
            F_code.append(''.join(parts))
        else:
            F_code.append(parts[0] + '  &')
            F_code.append(1)
            for arg in parts[1:-1]:
                F_code.append(arg + ',  &')
            F_code.append(parts[-1])
            F_code.append(-1)

    def write_module(self, library, cls):
        """ Write Fortran wrapper module.
        """
        node = cls or library
        options = node['options']
        fmt_node = node['fmt']
        fname = fmt_node.F_impl_filename
        module_name = fmt_node.F_module_name

        output = []

        if options.doxygen:
            self.write_doxygen_file(output, fname, library, cls)
        self._create_splicer('file_top', output)

        output.append('module %s' % module_name)
        output.append(1)

        if options.F_module_per_class:
            # XXX this will have some problems because of forward declarations
            for mname, only in node['F_module_dependencies']:
                if mname == module_name:
                    continue
                if only:
                    output.append('use %s, only : %s' % (
                            mname, ', '.join(only)))
                else:
                    output.append('use %s' % mname)
            output.extend(self.use_stmts)
            output.append('implicit none')
            output.append('')
        else:
            output.append('use, intrinsic :: iso_c_binding, only : C_PTR')
            self._create_splicer('module_use', output)
            output.extend(self.use_stmts)
            output.append('implicit none')
            output.append('')
            self._create_splicer('module_top', output)

        # XXX output.append('! splicer push class')
        output.extend(self.f_type_decl)
        # XXX  output.append('! splicer pop class')
        output.append('')

        # Interfaces for operator overloads
        if self.operator_map:
            ops = sorted(self.operator_map)
            for op in ops:
                output.append('')
                output.append('interface operator (%s)' % op)
                output.append(1)
                for opfcn in self.operator_map[op]:
                    output.append('module procedure %s' % opfcn)
                output.append(-1)
                output.append('end interface')
            output.append('')

        output.extend(self.c_interface)
        output.extend(self.generic_interface)

        # Insert any helper functions needed
        # (They are duplicated in each module)
        helper_source = []
        if self.f_helper:
            find_all_helpers(self.f_helper)

            helpers = sorted(self.f_helper)
            private_names = []
            interface_lines = []
            for helper in helpers:
                helper_info = whelpers.FHelpers[helper]
                private_names.extend(helper_info.get('private', []))
                lines = helper_info.get('interface', None)
                if lines:
                    interface_lines.append(lines)
                lines = helper_info.get('source', None)
                if lines:
                    helper_source.append(lines)
            if private_names:
                output.append('')
                output.append('private ' + ', '.join(private_names))
            output.extend(interface_lines)

        output.append(-1)
        output.append('')
        output.append('contains')
        output.append(1)

        output.extend(self.impl)

        output.extend(self.operator_impl)

        output.extend(helper_source)

        output.append(-1)
        output.append('')
        output.append('end module %s' % module_name)

        self.config.ffiles.append(
            os.path.join(self.config.c_fortran_dir, fname))
        self.write_output_file(fname, self.config.c_fortran_dir, output)

    def write_c_helper(self):
        """ Write C helper functions that will be used by the wrappers.
        """
        pass


def find_all_helpers(helpers, check=None):
    """Find all helper functions recursively.
    A helper function is required by some argument/result conversions.
    """

    if check is None:
        # do all top level helpers
        # Copy initial keys since helpers may change
        keys = list(helpers.keys())
        for check in keys:
            for name in whelpers.FHelpers[check].get('f_helper', []):
                find_all_helpers(helpers, name)
    else:
        if check not in helpers:
            helpers[check] = True
            for name in whelpers.FHelpers[check].get('f_helper', []):
                find_all_helpers(helpers, name)
