#!/bin/env python3
"""
Generate Fortran bindings for C++ code.

module {F_module_name}

type {F_derived_name}
  type(C_PTR) {F_derived_member}

constains
  procedure :: {F_name_method} => {F_name_impl}
  generic :: {F_name_generic} => {F_name_method}, ...
end type

interface
  {F_C_pure_clause}{F_C_subprogram} {F_C_name}({F_C_arguments}){F_C_result_clause} &
      bind(C, name="{C_name}")
    {arg_c_decl}
  end {F_C_subprogram} {F_C_name}

end interface

contains

 {F_pure_clause} {F_subprogram} {F_name_impl}({F_arguments}){F_result_clause}
     {F_C_name}({F_arg_c_call_tab})
 end {F_subprogram} {F_name_impl}

----------
TODO:
  intent is kludged for now.  They're all intent(IN) because ifort
  requires them for pure functions

"""
from __future__ import print_function

import util
from util import wformat, append_format

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

    def _begin_output_file(self):
        """Start a new class for output"""
        self.f_type_decl = []
        self.c_interface = []
        self.generic_interface = []
        self.impl = []         # implementation, after contains
        self.operator_impl = []
        self.operator_map = {} # list of function names by operator
                               # {'.eq.': [ 'abc', 'def'] }
        self.c_interface.append('interface')
        self.c_interface.append(1)

    def _end_output_file(self):
        self.c_interface.append(-1)
        self.c_interface.append('end interface')

    def _begin_class(self):
        self.f_type_generic = {} # look for generic methods
        self.type_bound_part = []


    def _c_type(self, arg):
        """
        Return the Fortran type, and array attribute
        pass-by-value default

        attributes:
        ptr - True = pass-by-reference
        reference - True = pass-by-reference

        """
        t = []
        typedef = self.typedef.get(arg['type'], None)
        if typedef is None:
            raise RuntimeError("No such type %s" % arg['type'])
        is_ptr = (arg['attrs'].get('ptr', False) or
                  arg['attrs'].get('reference', False))
        intent = arg['attrs'].get('intent', None)
        if intent:
            intent_str = ', intent(%s)' % intent.upper()
        else:
            intent_str = ''
        is_value = arg['attrs'].get('value', False)

        typ = typedef.c_fortran
        if typedef.base == 'string':
            return (typ + intent_str, True)  # is array
        else:
            #        if arg['attrs'].get('const', False):
            #            t.append('const')
            t.append(typ)
            if is_value:
                t.append(', value')
            t.append(intent_str)
            return (''.join(t), arg['attrs'].get('array', False))

    def _c_decl(self, arg, name=None):
        """
        Return the Fortran declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ, arr = self._c_type(arg)
        rv = typ + ' :: ' + ( name or arg['name'] )
        if arr:
            rv += '(*)'
        return rv

    def _f_type(self, arg, default=None):
        """
        Return the Fortran type, and array attribute
        pass-by-value default

        attributes:
        ptr - True = pass-by-reference
        reference - True = pass-by-reference

        """
        t = []
        typedef = self.typedef.get(arg['type'], None)
        if typedef is None:
            raise RuntimeError("No such type %s" % arg['type'])

        typ = typedef.f_type
        if typedef.base == 'string':
            return (typ, False)  # not array
        else:
            #        if arg['attrs'].get('const', False):
            #            t.append('const')
            t.append(typ)
            if default is None:
                default = arg['attrs'].get('default', '')
            if default != '':
                t.append('optional')
#            if not (arg['attrs'].get('ptr', False) or
#                    arg['attrs'].get('reference', False)):
#                t.append(', value')
            return (', '.join(t), arg['attrs'].get('array', False))

    def _f_decl(self, arg, name=None, default=None):
        """
        Return the Fortran declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ, arr = self._f_type(arg, default=default)
        rv = typ + ' :: ' + ( name or arg['name'] )
        if arr:
            rv += '(*)'
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
                self._end_output_file()
                self.write_module(node)
                self._begin_output_file()
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
                self._end_output_file()
                self.write_module(self.tree)

        if not options.F_module_per_class:
            # put all functions and classes into one module
            self.tree['F_module_dependencies'] = []
            self._end_output_file()
            self.write_module(self.tree)

        self.write_c_helper()

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        unname = util.un_camel(name)
        typedef = self.typedef[name]

        fmt_class = node['fmt']

        fmt_class.F_derived_name = typedef.fortran_derived

        # wrap methods
        self._push_splicer(fmt_class.cpp_class)
        self._push_splicer('method')
        for method in node['methods']:
            self.wrap_function(node, method)
        self._pop_splicer('method')
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
                wformat('type(C_PTR) {F_derived_member}', fmt_class),
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
        self.overload_compare(fmt_class, '.eq.', fmt_class.lower_class + '_eq',
                              wformat('c_associated(a%{F_derived_member}, b%{F_derived_member})', fmt_class))
#        self.overload_compare(fmt_class, '==', fmt_class.lower_class + '_eq', None)
        self.overload_compare(fmt_class, '.ne.', fmt_class.lower_class + '_ne',
                              wformat('.not. c_associated(a%{F_derived_member}, b%{F_derived_member})', fmt_class))
#        self.overload_compare(fmt_class, '/=', fmt_class.lower_class + '_ne', None)
        
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
        append_format(operator, 'type({F_derived_name}), intent(IN) ::a,b', fmt)
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
        if 'decl' in node:
            self.log.write(" {0} {1[decl]}\n".format(cls_function, node))
        else:
            self.log.write(" {0} {1[result][name]}\n".format(cls_function, node))

        if options.wrap_c:
            self.wrap_function_interface(cls, node)
        if options.wrap_fortran:
            self.wrap_function_impl(cls, node)

    def wrap_function_interface(self, cls, node):
        """
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
        result_is_ptr = result['attrs'].get('ptr', False)

        if node.get('return_this', False):
            result_type = 'void'
            result_is_ptr = False

        result_typedef = self.typedef[result_type]
        is_ctor  = node['attrs'].get('constructor', False)
        is_const = node['attrs'].get('const', False)
        is_pure  = node['attrs'].get('pure', False)

        arg_c_names = [ ]  # argument names for functions
        arg_c_decl = [ ]   # declaraion of argument names

        # find subprogram type
        # compute first to get order of arguments correct.
        # Add 
        if result_type == 'void' and not result_is_ptr:
            #  void=subroutine   void *=function
            fmt.F_C_subprogram = 'subroutine'
        else:
            fmt.F_C_subprogram = 'function'
            fmt.F_C_result_clause = ' result(%s)' % fmt.F_result
            if is_pure or func_is_const:
                fmt.F_C_pure_clause = 'pure '

        if cls:
            # Add 'this' argument
            if not is_ctor:
                arg_c_names.append(fmt.C_this)
                arg_c_decl.append('type(C_PTR), value, intent(IN) :: ' + fmt.C_this)

        for arg in node['args']:
            # default argument's intent
            # XXX look at const, ptr
            arg_typedef = self.typedef[arg['type']]
            fmt.var = arg['name']
            attrs = arg['attrs']
            if 'intent' not in attrs:
                attrs['intent'] = 'in'
            if 'value' not in attrs:
                attrs['value'] = True

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

            len_trim = arg['attrs'].get('len_trim', None)
            if len_trim:
                if len_trim is True:
                    len_trim = 'L' + arg['name']
                arg_c_names.append(len_trim)
                arg_c_decl.append('integer(C_INT), value, intent(IN) :: %s' % len_trim)

        fmt.F_C_arguments = options.get('F_C_arguments', ', '.join(arg_c_names))

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
                '{F_C_pure_clause}{F_C_subprogram} {F_C_name}({F_C_arguments}){F_C_result_clause} &',
                fmt))
        c_interface.append(2)  # extra indent for continued line
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
        fmt = util.Options(fmt_func)

        # look for C routine to wrap
        # usually the same node unless it is a generic function
        if 'PTR_F_C_index' in fmt_func:
            C_node = self.tree['function_index'][fmt_func.PTR_F_C_index]
            if len(node['args']) != len(C_node['args']):
                raise RuntimeError("Argument mismatch between Fortran and C functions")
        else:
            C_node = node
        fmt.F_C_name = C_node['fmt'].F_C_name

        func_is_const = node['attrs'].get('const', False)

        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)

        if node.get('return_this', False):
            result_type = 'void'
            result_is_ptr = False

        result_typedef = self.typedef[result_type]
        is_ctor  = node['attrs'].get('constructor', False)
        is_dtor  = node['attrs'].get('destructor', False)
        is_const = result['attrs'].get('const', False)

        # Special case some string handling
        if result_typedef.base == 'string' and \
                options.get('F_string_result_as_arg', False):
            # convert function into subroutine with argument for result
            result_string = True

            # Use the result_as_arg typedef
            result_typedef = self.typedef[result_typedef.name + '_result_as_arg']
            fmt.result_arg = options.F_string_result_as_arg
        else:
            result_string = False

        fmt.F_instance_ptr = wformat('{F_this}%{F_derived_member}', fmt)

        arg_c_call = []      # arguments to C function

        arg_f_names = [ ]
        arg_f_decl = [ ]
        arg_f_use  = [ 'use iso_c_binding' ]  # XXX totally brain dead for now

        # find subprogram type
        # compute first to get order of arguments correct.
        # Add 
        if result_type == 'void' and not result_is_ptr:
            #  void=subroutine   void *=function
            subprogram = 'subroutine'
        elif result_string:
            subprogram = 'subroutine'
        else:
            subprogram = 'function'
            fmt.F_result_clause = ' result(%s)' % fmt.F_result
        fmt.F_subprogram    = subprogram

        if cls:
            # Add 'this' argument
            if not is_ctor:
                arg_c_call.append(fmt.F_instance_ptr)
                arg_f_names.append(fmt.F_this)
                arg_f_decl.append(wformat(
                        'class({F_derived_name}) :: {F_this}',
                        fmt))

        optional = []
        c_args = C_node['args']
        for i, arg in enumerate(node['args']):
            # process Fortran function arguments first
            # default argument's intent
            # XXX look at const, ptr
            attrs = arg['attrs']
            if 'intent' not in attrs:
                attrs['intent'] = 'in'
            if 'value' not in attrs:
                attrs['value'] = True

            fmt.var = arg['name']
            fmt.tmp_var = 'tmp_' + fmt.var
            arg_f_names.append(fmt.var)
            arg_f_decl.append(self._f_decl(arg))

            if 'default' in attrs:
                arg_f_decl.append(self._f_decl(arg, name=fmt.tmp_var, default=''))
                fmt.default_value = attrs['default']
                optional.extend([
                        wformat('if (present({var})) then', fmt),
                        1,
                        wformat('{tmp_var} = {var}', fmt),
                        -1,
                        'else',
                        1,
                        wformat('{tmp_var} = {default_value}', fmt),
                        -1,
                        'endif'])
                fmt.var = fmt.tmp_var  # pass tmp to C function

            arg_typedef = self.typedef[arg['type']]

            if arg_typedef.f_pre_decl:
                append_format(optional, arg_typedef.f_pre_decl, fmt)
            if arg_typedef.f_pre_call:
                append_format(optional, arg_typedef.f_pre_call, fmt)
            if arg_typedef.f_use_tmp:
                fmt.var = fmt.tmp_var

            # Then C function arguments
            # match to corresponding C argument -- must have same number of args.
            # may have different types, like generic
            # or different attributes, like adding +len to string args
            c_arg = c_args[i]
            arg_typedef = self.typedef[c_arg['type']]

            # Attributes   None=skip, True=use default, else use value
            len_trim = c_arg['attrs'].get('len_trim', None)
            if len_trim:
#                fmt.len_trim_var = 'L' + arg['name']
                if len_trim is True:
                    len_trim = 'len_trim({var})'
                append_format(arg_c_call, '{var}', fmt)
                append_format(arg_c_call, len_trim, fmt)
# XXX need both, cast then fortran_to_c
            elif c_arg['type'] != arg['type']:
                append_format(arg_c_call, arg_typedef.f_cast, fmt)
            else:
                append_format(arg_c_call, arg_typedef.fortran_to_c, fmt)

        if result_string:
            arg_f_names.append(fmt.result_arg)
            if result_typedef.f_rv_decl:
                append_format(arg_f_decl, result_typedef.f_rv_decl, fmt)
            if result_typedef.f_pre_decl:
                append_format(arg_f_decl, result_typedef.f_pre_decl, fmt)

        fmt.F_arg_c_call = ', '.join(arg_c_call)
        fmt.F_arg_c_call_tab = '\t' + '\t'.join(arg_c_call) # use tabs to insert continuations
        fmt.F_arguments = options.get('F_arguments', ', '.join(arg_f_names))

        # declare function return value after arguments
        # since arguments may be used to compute return value
        # (for example, string lengths)
        if subprogram == 'function':
#            if func_is_const:
#                fmt.F_pure_clause = 'pure '
            if result_typedef.base == 'string':
                # special case returning a string
                rvlen = result['attrs'].get('len', None)
                if rvlen is None:
                    rvlen = wformat('strlen_ptr({F_C_name}({F_arg_c_call}))', fmt)
                fmt.rvlen = wformat(rvlen, fmt)
                arg_f_decl.append(
                    wformat('character(kind=C_CHAR, len={rvlen}) :: {F_result}',
                            fmt))
            else:
                arg_f_decl.append(self._f_decl(result, name=fmt.F_result))

        if not is_ctor:
            # Add method to derived type
            F_name_method = fmt.F_name_method
#            if not fmt.get('CPP_template', None):
            if not fmt.get('CPP_return_templated', False):
                # if return type is templated in C++, then do not set up generic
                # since only the return type may be different (ex. getValue<T>())
                self.f_type_generic.setdefault(fmt.F_name_generic,[]).append(F_name_method)
            self.type_bound_part.append('procedure :: %s => %s' % (
                    F_name_method, fmt.F_name_impl))

        # body of function
        splicer_code = self.splicer_stack[-1].get(fmt_func.F_name_method, None)
        if 'F_code' in options:
            F_code = [   wformat(options.F_code, fmt) ]
        elif splicer_code:
            F_code = splicer_code
        else:
            F_code = []
            if is_ctor:
                line1 = wformat('{F_result}%{F_derived_member} = {F_C_name}({F_arg_c_call_tab})', fmt)
                self.append_method_arguments(F_code, line1)
            elif result_string:
                line1 = wformat(result_typedef.f_return_code, fmt)
                self.append_method_arguments(F_code, line1)
            elif subprogram == 'function':
                line1 = wformat(result_typedef.f_return_code, fmt)
                self.append_method_arguments(F_code, line1)
            else:
                line1 = wformat('call {F_C_name}({F_arg_c_call_tab})', fmt)
                self.append_method_arguments(F_code, line1)

            if result_typedef.f_post_call:
                # adjust return value or cleanup
                append_format(F_code, result_typedef.f_post_call, fmt)
            if is_dtor:
                F_code.append(wformat('{F_this}%{F_derived_member} = C_NULL_PTR', fmt))

        impl = self.impl
        impl.append('')
        impl.append(wformat('{F_subprogram} {F_name_impl}({F_arguments}){F_result_clause}', fmt))
        impl.append(1)
        impl.extend(arg_f_use)
        impl.append('implicit none')
        impl.extend(arg_f_decl)
        impl.extend(optional)
        self._create_splicer(fmt.F_name_method, impl, F_code)
        impl.append(-1)
        impl.append(wformat('end {F_subprogram} {F_name_impl}', fmt))

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

    def write_module(self, node):
        options = node['options']
        fmt_class = node['fmt']
        fname = fmt_class.F_impl_filename
        module_name = fmt_class.F_module_name

        output = []
        output.append('module %s' % module_name)
        output.append(1)

        output.append('use fstr_mod')
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
            output.append('implicit none')
            output.append('')
        else:
            output.append('use, intrinsic :: iso_c_binding, only : C_PTR')
            output.append('implicit none')
            output.append('')
            self._create_splicer('module_top', output)

#X        output.append('! splicer push class')
        output.extend(self.f_type_decl)
#X        output.append('! splicer pop class')
        output.append('')

        # Interfaces for operator overloads
        if self.operator_map:
            ops = self.operator_map.keys()
            ops.sort()
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

        output.append(-1)
        output.append('')
        output.append('contains')
        output.append(1)

        output.extend(self.impl)

        output.extend(self.operator_impl)

        output.append(-1)
        output.append('')
        output.append('end module %s' % module_name)

        self.write_output_file(fname, self.config.binary_dir, output)

    def write_c_helper(self):
        """ Write C helper functions that will be used by the wrappers.
        """
        pass
