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

"""
from __future__ import print_function

import util
import fwrap_util

import os

wformat = util.wformat

class Wrapf(object):
    """Generate Fortran bindings.
    """

    def __init__(self, tree, config):
        self.tree = tree    # json tree
        self.config = config
        self.log = config.log
        self.typedef = tree['typedef']

    def _begin_file(self):
        """Start a new class for output"""
        self.f_type = []
        self.f_type_generic = {} # look for generic methods
        self.c_interface = []
        self.impl = []         # implementation

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

        typ = typedef.c_fortran
        if typedef.base == 'string':
            return (typ, True)   # is array
        else:
            #        if arg['attrs'].get('const', False):
            #            t.append('const')
            t.append(typ)
            if not is_ptr:
                t.append(', value')
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

        self.fmt_library = dict(
            lower_library = options.library.lower()
            )

        self._begin_file()
        for node in self.tree['classes']:
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
                self._begin_file()

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

        fmt_class = self.fmt_class = dict(
            cpp_class = name,
            lower_class = name.lower(),
            upper_class = name.upper(),
            F_derived_name = typedef.fortran_derived,
            F_this = node['options'].F_this,
            )
        fmt_class.update(self.fmt_library)

        unname = util.un_camel(name)
        self.f_type.extend([
                '',
                wformat('type {F_derived_name}', fmt_class),
                1,
                wformat('type(C_PTR) {F_this}', fmt_class),
                -1, 'contains', 1,
                ])

        for method in node['methods']:
            self.wrap_method(node, method)

        # Look for generics
        for key in sorted(self.f_type_generic.keys()):
            methods = self.f_type_generic[key]
            if len(methods) > 1:
                self.f_type.append('generic :: %s => %s' % (
                        key, ', '.join(methods)))

        self.f_type.extend([
                -1,
                 wformat('end type {F_derived_name}', fmt_class),
                 ''
                 ])

    def wrap_method(self, cls, node):
        """
        cls  - class node
        node - function node
        """
        if 'decl' in node:
            self.log.write("method {1[decl]}\n".format(self, node))
        else:
            self.log.write("method {1[result][name]}\n".format(self, node))

        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)
        result_typedef = self.typedef[result_type]
        C_this = node['options'].C_this
        F_result = node['options'].F_result
        is_ctor  = result['attrs'].get('constructor', False)
        is_dtor  = result['attrs'].get('destructor', False)
        is_const = result['attrs'].get('const', False)

        fmt_func = dict(
            method_name=result['name'],
            underscore_name=util.un_camel(result['name']),
            method_suffix = node.get('method_suffix', ''),

            F_this=node['options'].F_this,
            F_result=F_result,
            C_name=node['C_name'],
            )
        fmt_func['F_obj'] = wformat('{F_this}%{F_this}', fmt_func)
        fmt_func.update(self.fmt_class)

        if 'F_C_name' not in node:
            node['F_C_name'] = node['C_name'].lower()
        fmt_func['F_C_name'] = node['F_C_name']

        util.eval_templates(
            ['F_name_impl', 'F_name_method', 'F_name_generic'],
            node, fmt_func)

        arg_c_names = [ ]
        arg_c_decl = [ ]
        arg_c_call = []      # arguments to C function

        arg_f_names = [ ]
        arg_f_decl = [ ]

        # find subprogram type
        # compute first to get order of arguments correct.
        # Add 
        pure_clause = ''
        if result_type == 'void' and not result_is_ptr:
            #  void=subroutine   void *=function
            subprogram = 'subroutine'
            result_clause = ''
        else:
            subprogram = 'function'
            result_clause = ' result(%s)' % F_result
            if is_const:
                pure_clause = 'pure '
        fmt_func['subprogram'] = subprogram
        fmt_func['result_clause'] = result_clause
        fmt_func['pure_clause'] = pure_clause

        # Add 'this' argument
        if not is_ctor and (is_dtor or subprogram == 'function'):
            arg_c_names.append(C_this)
            arg_c_decl.append('type(C_PTR), value :: ' + C_this)
            arg_c_call.append(fmt_func['F_obj'])
            arg_f_names.append(fmt_func['F_this'])
            if is_dtor:
                arg_f_decl.append(wformat(
                        'type({F_derived_name}) :: {F_this}',
                        fmt_func))
            else:
                arg_f_decl.append(wformat(
                        'class({F_derived_name}) :: {F_this}',
                        fmt_func))

        for arg in node.get('args', []):
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
                fmt_func['rvlen'] = rvlen
                arg_c_decl.append('type(C_PTR) %s' % F_result)
                arg_f_decl.append(
                    wformat('character(kind=C_CHAR, len={rvlen}) :: {F_result}',
                            fmt_func))
                arg_f_decl.append('type(C_PTR) :: rv_ptr')
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
            F_name_method = node['F_name_method']
            F_name_generic = node['F_name_generic']
            F_name_impl = node['F_name_impl']
            self.f_type_generic.setdefault(F_name_generic,[]).append(F_name_method)
            self.f_type.append('procedure :: %s => %s' % (
                    F_name_method, F_name_impl))

        fmt_func['arg_c_call'] = ', '.join(arg_c_call)

        if 'F_C_arguments' not in node:
            node['F_C_arguments'] = ', '.join(arg_c_names)
        fmt_func['F_C_arguments'] = node['F_C_arguments']

        if 'F_arguments' not in node:
            node['F_arguments'] = ', '.join(arg_f_names)
        fmt_func['F_arguments'] = node['F_arguments']

        if 'F_code' not in node:
            lines = []
            if is_ctor:
                lines.append(wformat('{F_result}%{F_this} = {F_C_name}({arg_c_call})', fmt_func))
            elif subprogram == 'function':
                fmt = result_typedef.f_return_code
                lines.append(wformat(fmt, fmt_func))

            else:
                lines.append(wformat('call {F_C_name}({arg_c_call})', fmt_func))

            if is_dtor:
                lines.append(wformat('{F_this}%{F_this} = C_NULL_PTR', fmt_func))

            node['F_code'] = '\n'.join(lines)

        c_interface = self.c_interface
        c_interface.append('')
        c_interface.append(wformat(
                '{pure_clause}{subprogram} {F_C_name}({F_C_arguments}){result_clause} bind(C, name="{C_name}")',
                fmt_func))
        c_interface.append(1)
        c_interface.append('use iso_c_binding')
        c_interface.append('implicit none')
        c_interface.extend(arg_c_decl)
        c_interface.append(-1)
        c_interface.append(wformat('end {subprogram} {F_C_name}', fmt_func))

        impl = self.impl
        impl.append('')
        impl.append(wformat('{subprogram} {F_name_impl}({F_arguments}){result_clause}', fmt_func))
        impl.append(1)
        impl.append('implicit none')
        impl.extend(arg_f_decl)
        impl.append('! splicer begin')
        impl.append(node['F_code'])
        impl.append('! splicer end')
        impl.append(-1)
        impl.append(wformat('end {subprogram} {F_name_impl}', fmt_func))

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
        fname = node['fmt']['F_impl_filename']
        module_name = node['fmt'].F_module_name
        fp = open(os.path.join(self.config.binary_dir, fname), 'w')
        self.write_copyright(fp)

        output = []
        output.append('module %s' % module_name)
        output.append(1)

        output.append('use fstr_mod')
        for mname, only in node['F_module_dependencies']:
            if mname == module_name:
                continue
            if only:
                output.append('use %s, only : %s' % (
                        mname, ', '.join(only)))
            else:
                output.append('use %s' % mname)

        output.extend(self.f_type)

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

