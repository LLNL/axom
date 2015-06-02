#!/bin/env python3

"""
generate language bindings
"""
from __future__ import print_function

import copy
import os
import json
import argparse
import yaml

import util
import parse_decl
import wrapc
import wrapf


class Schema(object):
    """
    Verify that the input dictionary has the correct fields.
    Create defaults for missing fields.
    """
    def __init__(self, config):
        self.config = config

    def push_options(self, node):
        """ Push a new set of options.
        Copy current options, then update with new options.
        """
        if 'options' in node and \
                node['options'] is not None:
            if not isinstance(node['options'], dict):
                raise TypeError("options must be a dictionary")
            new = copy.deepcopy(self.options_stack[-1])
            new.update(node['options'])
        else:
            new = self.options_stack[-1]
        self.options_stack.append(new)
        node['options'] = new
        return new

    def pop_options(self):
        self.options_stack.pop()

    def update_options(self, options, d):
        """ update options with values from d,
        if key is not already in options.
        """
        for key in d:
            if key not in options:
                options[key] = d[key]

    def check_schema(self):
        node = self.config

        # default options
        def_options = dict(
            C_prefix='',
            C_this='self',   # object argument name
            F_this='obj',    # object argument name
            F_result='rv',   # function result
#            module_name='gen_module',  # Fortran module name
            C_header_filename_template='wrap{cpp_class}.h',
            C_impl_filename_template='wrap{cpp_class}.cpp',
            F_impl_filename_template='wrapf{cpp_class}.f',
            F_module_name_template='{lower_class}_mod',

            C_name_method_template='{C_prefix}{lower_class}_{underscore_name}{method_suffix}',
            F_name_method_template='{underscore_name}{method_suffix}',
            F_name_generic_template='{underscore_name}',
            F_name_impl_template  ='{lower_class}_{underscore_name}{method_suffix}',
            )
        if 'options' in node:
            def_options.update(node['options'])
        self.options_stack = [ def_options ]
        node['options'] = def_options

        def_types = dict(
            void    = dict(
                c       = 'void',
                cpp     = 'void',
#                fortran = 'subroutine',
                c_fortran = 'type(C_PTR)',
                fortran = 'type(C_PTR)',
                ),
            int    = dict(
                c       = 'int',
                cpp     = 'int',
                c_fortran = 'integer(C_INT)',
                fortran   = 'integer(C_INT)',
                ),
            long   = dict(
                c       = 'long',
                cpp     = 'long',
                c_fortran = 'integer(C_LONG)',
                fortran   = 'integer(C_LONG)',
                ),
            size_t   = dict(
                c       = 'size_t',
                cpp     = 'size_t',
                c_header = 'stdlib.h',
                c_fortran = 'integer(C_SIZE_T)',
                fortran   = 'integer(C_SIZE_T)',
                ),
            bool   = dict(
                c       = 'bool',
                cpp     = 'bool',
                c_fortran = 'logical(C_BOOL)',
                fortran   = 'logical',
                f_return_code = '{F_result} = bool2logical({F_C_name}({arg_c_call}))',
                ),
            string = dict(  # implies null terminated string
                c   = 'char',    # XXX - char *
                cpp = 'std::string',
                cpp_to_c = '{var}.c_str()',  # . or ->
                c_fortran  = 'character(kind=C_CHAR)',
                fortran = 'character(*)',
                fortran_to_c = 'trim({var}) // C_NULL_CHAR',
#                f_module = dict(iso_c_binding = [ 'C_NULL_CHAR' ]),
                f_module = dict(iso_c_binding=None),
                f_return_code = '{F_result} = fstr({F_C_name}({arg_c_call}))',
                base = 'string',
                ),
            )
        def_types['std::string'] = def_types['string']
        if 'typedef' in node:
            def_types.update(node['typedef'])

        # set some default members
        for typ in def_types.values():
            if 'base' not in typ:
                typ['base'] = 'unknown'

        node['typedef'] = def_types
        self.typedef = node['typedef']

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

        for cls in node:
            self.check_class_dependencies(cls)
        
    def check_class(self, node):
        if 'name' not in node:
            raise RuntimeError('Expeced name for class')
        name = node['name']

        options = self.push_options(node)

        fmt_dict = dict(
            cpp_class = name,
            lower_class = name.lower(),
            upper_class = name.upper(),
            C_prefix = options['C_prefix'],
            )
        util.eval_templates(
            ['C_header_filename',
             'C_impl_filename',
             'F_impl_filename',
             'F_module_name'],
            node, fmt_dict)

        # create typedef for each class before generating code
        # this allows classes to reference each other
        if name not in self.typedef:
#            unname = util.un_camel(name)
            unname = name.lower()
            cname = node['options']['C_prefix'] + unname
            self.typedef[name] = dict(
                cpp = name,
                c = cname,
                c_to_cpp = 'static_cast<%s{ptr}>({var})' % name,
                c_fortran = 'type(C_PTR)',
                fortran = 'type(%s)' % unname,
                fortran_type = unname,
                fortran_to_c = '{var}%obj',
                # XXX module name may not conflict with type name
                f_module = {node['F_module_name']:[unname]},

                # return from C function
#                f_c_return_decl = 'type(CPTR)' % unname,
                f_return_code = '{F_result}%{F_this} = {F_C_name}({arg_c_call})',

                # allow forward declarations to avoid recursive headers
                forward = name,
                base = 'wrapped',
                )

        method_dict = {}
        methods = node.setdefault('methods', [])
        for method in methods:
            if not isinstance(method, dict):
                raise TypeError("classes[n]['methods'] must be a dictionary")
            self.check_function(method)
            method_dict.setdefault(method['result']['name'], []).append(method)

        # look for function overload and compute method_suffix
        for mname, methods in method_dict.items():
            if len(methods) > 1:
                for i, method in enumerate(methods):
                    if 'method_suffix' not in method:
                        method['method_suffix'] = '_%d' % i

        self.pop_options()

    def check_function(self, node):
        self.push_options(node)

        if 'decl' in node:
            # parse decl and add to dictionary
            values = parse_decl.check_decl(node['decl'])
            util.update(node, values)  # recursive update
        if 'result' not in node:
            raise RuntimeError("Missing result")
        result = node['result']
        if 'name' not in result:
            raise RuntimeError("Missing result.name")
        if 'type' not in result:
            raise RuntimeError("Missing result.type")

        # docs
        self.pop_options()

    def check_functions(self, node):
        if not isinstance(node, list):
            raise TypeError("functions must be a list")
        for func in node:
            self.check_function(func)

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
            if 'f_module' in typ:
                for mname, only in typ['f_module'].items():
                    module = modules.setdefault(mname, {})
                    if only:  # Empty list means no ONLY clause
                        for oname in only:
                            module[oname] = True
        F_modules = []  # array of tuples ( name, (only1, only2) )
        for mname in sorted(modules):
            F_modules.append((mname, sorted(modules[mname])))
        node['F_module_dependencies'] = F_modules

    def check_function_dependencies(self, node, used_types):
        result = node['result']
        rv_type = result['type']
        if rv_type not in self.typedef:
            raise RuntimeError("Unknown type {} for function decl: {}".format(rv_type, node['decl']))
        result_typedef = self.typedef[rv_type]
        # XXX - make sure it exists
        used_types[result['type']] = result_typedef
        for arg in node.get('args', []):
            used_types[arg['type']] = self.typedef[arg['type']]


if __name__ == '__main__':

    appname = 'modulator3'
    appver = '0.1'

    parser = argparse.ArgumentParser(prog=appname)
    parser.add_argument('--version', action='version', version=appver)
    parser.add_argument('filename', nargs='*',
                        help='Input file to process')

    args = parser.parse_args()

    if len(args.filename) == 0:
        raise RuntimeError("Must give at least one input file")

    basename = os.path.splitext(os.path.basename(args.filename[0]))[0]
    log = open(basename + '.log', 'w')

    # accumulated input
    all = {}

    for filename in args.filename:
        root, ext = os.path.splitext(filename)
        if ext == '.yaml':
            fp = open(filename, 'r')
            d = yaml.load(fp.read())
            fp.close()
            all.update(d)
        else:
            raise RuntimeError("File must end in .yaml for now")

#    print(all)



    Schema(all).check_schema()

    wrapc.Wrapc(all, log).wrap()

    wrapf.Wrapf(all, log).wrap()

    fp = open(basename + '.json', 'w')
    json.dump(all, fp, sort_keys=True, indent=4)
    fp.close()

    log.close()
