#!/bin/env python3
"""
Generate Fortran bindings for C++ code.

module {F_module_name}

type {F_derived_name}
  type(C_PTR) {F_this}
constains
  procedure :: {F_name_method} => {F_name_impl}
  generic :: {F_name_generic} => {F_name_method}, ...
end type


TODO:
  intent is kludged for now.  They're all intent(IN) because ifort
  requires them for pure functions

"""
from __future__ import print_function

import util
import fwrap_util

import os

wformat = util.wformat

class Wrapf(object):
    """Generate Fortran bindings.
    """

    def __init__(self, tree, config, splicers):
        self.tree = tree    # json tree
        self.config = config
        self.splicers = splicers
        self.log = config.log
        self.typedef = tree['typedef']
        self.splicer_stack = [ splicers ]

    def _begin_output_file(self):
        """Start a new class for output"""
        self.f_type_decl = []
        self.c_interface = []
        self.impl = []         # implementation

    def _begin_class(self):
        self.f_type_generic = {} # look for generic methods

    def _push_splicer(self, name):
        level = self.splicer_stack[-1].setdefault(name, {})
        self.splicer_stack.append(level)

    def _pop_splicer(self, name):
        # XXX maybe use name for error checking
        self.splicer_stack.pop()

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

        typ = typedef.c_fortran
        if typedef.base == 'string':
            return (typ, True)   # is array
        else:
            #        if arg['attrs'].get('const', False):
            #            t.append('const')
            t.append(typ)
            if not is_ptr:
                t.append(', value')
            if intent:
                t.append(', intent(%s)' % intent.upper())
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

    def _f_type(self, arg):
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
#            if not (arg['attrs'].get('ptr', False) or
#                    arg['attrs'].get('reference', False)):
#                t.append(', value')
            return (''.join(t), arg['attrs'].get('array', False))

    def _f_decl(self, arg, name=None):
        """
        Return the Fortran declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ, arr = self._f_type(arg)
        rv = typ + ' :: ' + ( name or arg['name'] )
        if arr:
            rv += '(*)'
        return rv

    def wrap_library(self):
        options = self.tree['options']
        fmt_library = self.tree['fmt']
        fmt_library.F_this = options.get('F_this', 'obj')
        fmt_library.F_result = options.get('F_result', 'rv')
        fmt_library.F_result_clause = ''
        fmt_library.F_pure_clause = ''

        self._begin_output_file()
        self._push_splicer('class')
        for node in self.tree['classes']:
            self._begin_class()
            self.c_interface.append('interface')
            self.c_interface.append(1)

            name = node['name']
            # how to decide module name, module per class
#            module_name = node['options'].setdefault('module_name', name.lower())
            self.wrap_class(node)
            self.c_interface.append(-1)
            self.c_interface.append('end interface')
            if options.F_module_per_class:
                self.write_module(node)
                self._begin_output_file()
        self._pop_splicer('class')

        if not options.F_module_per_class:
            # put all classes into one module
            self.tree['F_module_dependencies'] = []
            self.write_module(self.tree)

        for node in self.tree['functions']:
            self.wrap_function(node)

        fwrap_util.write_runtime(self.tree, self.config)

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        typedef = self.typedef[name]

        options = node['options']
        fmt_class = node['fmt']
        self._push_splicer(fmt_class.lower_class)
        fmt_class.F_derived_name = typedef.fortran_derived
        if 'F_this' in options:
            fmt_class.F_this = options.F_this

        unname = util.un_camel(name)
        self.f_type_decl.extend([
                '',
                wformat('type {F_derived_name}', fmt_class),
                1,
                wformat('type(C_PTR) {F_this}', fmt_class),
                -1, 'contains', 1,
                ])

        self.impl.append(wformat('! splicer push class.{lower_class}.method', fmt_class))
        self._push_splicer('method')
        for method in node['methods']:
            self.wrap_method(node, method)
        self._pop_splicer('method')
        self.impl.append('')
        self.impl.append(wformat('! splicer pop class.{lower_class}.method', fmt_class))

        # Look for generics
        for key in sorted(self.f_type_generic.keys()):
            methods = self.f_type_generic[key]
            if len(methods) > 1:
                self.f_type_decl.append('generic :: %s => %s' % (
                        key, ', '.join(methods)))

        self.f_type_decl.extend([
                -1,
                 wformat('end type {F_derived_name}', fmt_class),
                 ''
                 ])
        self._pop_splicer(fmt_class.lower_class)

    def wrap_method(self, cls, node):
        """
        cls  - class node
        node - function node
        """
        if 'decl' in node:
            self.log.write("method {1[decl]}\n".format(self, node))
        else:
            self.log.write("method {1[result][name]}\n".format(self, node))

        options = node['options']
        fmt_func = node['fmt']

        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)

        if node.get('return_this', False):
            result_type = 'void'
            result_is_ptr = False

        result_typedef = self.typedef[result_type]
        is_ctor  = result['attrs'].get('constructor', False)
        is_dtor  = result['attrs'].get('destructor', False)
        is_const = result['attrs'].get('const', False)

        fmt_func.F_obj = wformat('{F_this}%{F_this}', fmt_func)
        if 'F_this' in options:
            fmt_func.F_this = options.F_this
        if 'F_result' in options:
            fmt_func.F_result = options.F_result
        F_result = fmt_func.F_result
        C_this = fmt_func.C_this

        if 'F_C_name' in options:
            fmt_func.F_C_name = options.F_C_name
        else:
            fmt_func.F_C_name = fmt_func.C_name.lower()

        util.eval_template(options, fmt_func,
                            'F_name_impl', '{lower_class}_{underscore_name}{method_suffix}')
        util.eval_template(options, fmt_func,
                            'F_name_method', '{underscore_name}{method_suffix}')
        util.eval_template(options, fmt_func,
                            'F_name_generic', '{underscore_name}')

        arg_c_names = [ ]
        arg_c_decl = [ ]
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
        else:
            subprogram = 'function'
            fmt_func.F_result_clause = ' result(%s)' % F_result
            if is_const:
                fmt_func.F_pure_clause   = 'pure '
        fmt_func.F_subprogram    = subprogram

        # Add 'this' argument
        if not is_ctor:
            arg_c_names.append(C_this)
            arg_c_decl.append('type(C_PTR), value, intent(IN) :: ' + C_this)
            arg_c_call.append(fmt_func.F_obj)
            arg_f_names.append(fmt_func.F_this)
            if is_dtor:
                arg_f_decl.append(wformat(
                        'type({F_derived_name}) :: {F_this}',
                        fmt_func))
            else:
                arg_f_decl.append(wformat(
                        'class({F_derived_name}) :: {F_this}',
                        fmt_func))

        for arg in node.get('args', []):
            # default argument's intent
            # XXX look at const, ptr
            if 'intent' not in arg['attrs']:
                arg['attrs']['intent'] = 'in'

            arg_c_names.append(arg['name'])
            arg_c_decl.append(self._c_decl(arg))

            rrr = self.typedef[arg['type']].fortran_to_c
            arg_c_call.append(rrr.format(var=arg['name']))
            arg_f_names.append(arg['name'])
            arg_f_decl.append(self._f_decl(arg))

        # declare function return value after arguments
        # since arguments may be used to compute return value
        # (for example, string lengths)
        if subprogram == 'function':
            if result_typedef.base == 'string':
                # special case returning a string
                rvlen = result['attrs'].get('len', '1')
                fmt_func.rvlen = rvlen
                arg_c_decl.append('type(C_PTR) %s' % F_result)
                arg_f_decl.append(
                    wformat('character(kind=C_CHAR, len={rvlen}) :: {F_result}',
                            fmt_func))
            else:
                # XXX - make sure ptr is set to avoid VALUE
                arg_dict = dict(name=F_result,
                                type=result_type,
                                attrs=dict(ptr=True))
#                arg_c_decl.append(self._c_decl(result, name=F_result))
                arg_c_decl.append(self._c_decl(arg_dict))
                arg_f_decl.append(self._f_decl(result, name=F_result))

        if not is_ctor and not is_dtor:
            # Add method to derived type
            F_name_method = fmt_func.F_name_method
            self.f_type_generic.setdefault(fmt_func.F_name_generic,[]).append(F_name_method)
            self.f_type_decl.append('procedure :: %s => %s' % (
                    F_name_method, fmt_func.F_name_impl))

        fmt_func.F_arg_c_call = ', '.join(arg_c_call)
        fmt_func.F_C_arguments = options.get('F_C_arguments', ', '.join(arg_c_names))
        fmt_func.F_arguments = options.get('F_arguments', ', '.join(arg_f_names))

        # body of function
        splicer_code = self.splicer_stack[-1].get(fmt_func.F_name_method, None)
        if 'F_code' in options:
            F_code = [ options.F_code ]
        elif splicer_code:
            F_code = splicer_code
        else:
            F_code = []
            if is_ctor:
                F_code.append(wformat('{F_result}%{F_this} = {F_C_name}({F_arg_c_call})', fmt_func))
            elif subprogram == 'function':
                fmt = result_typedef.f_return_code
                F_code.append(wformat(fmt, fmt_func))
            else:
                F_code.append(wformat('call {F_C_name}({F_arg_c_call})', fmt_func))
            if is_dtor:
                F_code.append(wformat('{F_this}%{F_this} = C_NULL_PTR', fmt_func))

        c_interface = self.c_interface
        c_interface.append('')
        c_interface.append(wformat(
                '{F_pure_clause}{F_subprogram} {F_C_name}({F_C_arguments}){F_result_clause} &',
                fmt_func))
        c_interface.append(2)  # extra indent for continued line
        c_interface.append(wformat(
                'bind(C, name="{C_name}")',
                fmt_func))
        c_interface.append(-1)
        c_interface.append('use iso_c_binding')
        c_interface.append('implicit none')
        c_interface.extend(arg_c_decl)
        c_interface.append(-1)
        c_interface.append(wformat('end {F_subprogram} {F_C_name}', fmt_func))

        impl = self.impl
        impl.append('')
        impl.append(wformat('{F_subprogram} {F_name_impl}({F_arguments}){F_result_clause}', fmt_func))
        impl.append(1)
        impl.extend(arg_f_use)
        impl.append('implicit none')
        impl.extend(arg_f_decl)
        impl.append(wformat('! splicer begin {F_name_method}', fmt_func))
        impl.extend(F_code)
        impl.append(wformat('! splicer end {F_name_method}', fmt_func))
        impl.append(-1)
        impl.append(wformat('end {F_subprogram} {F_name_impl}', fmt_func))

    def wrap_function(self, node):
        """
        node - function node
        """
        self.log.write("function {1[decl]}\n".format(self, node))

    def write_copyright(self, fp):
        for line in self.tree.get('copyright', []):
            if line:
                fp.write('! ' + line + '\n')
            else:
                fp.write('!\n')

    def write_module(self, node):
        options = node['options']
        fmt_class = node['fmt']
        fname = fmt_class.F_impl_filename
        module_name = fmt_class.F_module_name
        fp = open(os.path.join(self.config.binary_dir, fname), 'w')
        self.write_copyright(fp)

        output = []
        output.append('module %s' % module_name)
        output.append(1)

        output.append('use fstr_mod')
        if node['options'].F_module_per_class:
            include_option = 'F_class_module_includes'
            # XXX this will have some problems because of forward declarations
            for mname, only in node['F_module_dependencies']:
                if mname == module_name:
                    continue
                if only:
                    output.append('use %s, only : %s' % (
                            mname, ', '.join(only)))
                else:
                    output.append('use %s' % mname)
        else:
            include_option = 'F_library_module_includes'
            output.append('use, intrinsic :: iso_c_binding, only : C_PTR')
        output.append('implicit none')
        for include in options.get(include_option, '').split():
            output.append('include "{}"'.format(include))

        output.extend(self.f_type_decl)

        output.extend(self.c_interface)

        output.append(-1)
        output.append('')
        output.append('contains')
        output.append(1)

        output.extend(self.impl)

        output.append(-1)
        output.append('')
        output.append('end module %s' % module_name)

        self.indent = 0
        self.write_lines(fp, output)
        fp.close()
        self.log.write("Close %s\n" % fname)
        print("Wrote", fname)

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

