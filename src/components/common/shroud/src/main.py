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
        node = self.config

        # default options
        def_options = util.Options(
            parent=None,

#            library='default_library',
            namespace='',
            cpp_header='',

            F_module_per_class=True,
            )
        if 'options' in node:
            def_options.update(node['options'])
        self.options_stack = [ def_options ]

        fmt_library = node['fmt'] = util.Options(None)
        fmt_library.library       = def_options.get('library', 'default_library')
        fmt_library.lower_library = fmt_library.library.lower()
        fmt_library.method_suffix = ''   # assume no suffix
        fmt_library.overloaded    = False
        fmt_library.C_prefix      = def_options.get('C_prefix', '')
        util.eval_template(def_options, fmt_library,
                           'C_header_filename', 'wrap{library}.h')
        util.eval_template(def_options, fmt_library,
                           'C_impl_filename', 'wrap{library}.cpp')
        self.fmt_stack.append(fmt_library)

        # default some options based on other options
        # All class/methods and functions may go into this file or just functions.
        util.eval_template(def_options, fmt_library,
                           'F_module_name', '{lower_library}_mod')
        util.eval_template(def_options, fmt_library,
                           'F_impl_filename', 'wrapf{lower_library}.f')

        node['options'] = def_options

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
                c_fortran = 'integer(C_INT)',
                f_type    = 'integer(C_INT)',
                PY_format = 'i',
                ),
            long   = util.Typedef('long',
                c_type    = 'long',
                cpp_type  = 'long',
                c_fortran = 'integer(C_LONG)',
                f_type    = 'integer(C_LONG)',
                PY_format = 'l',
                ),
            size_t   = util.Typedef('size_t',
                c_type    = 'size_t',
                cpp_type  = 'size_t',
                c_header  = 'stdlib.h',
                c_fortran = 'integer(C_SIZE_T)',
                f_type    = 'integer(C_SIZE_T)',
                ),
            bool   = util.Typedef('bool',
                c_type    = 'bool',
                cpp_type  = 'bool',
                c_fortran = 'logical(C_BOOL)',
                fortran_to_c  = 'logicaltobool({var})',
                f_type    = 'logical',
                f_return_code = '{F_result} = booltological({F_C_name}({F_arg_c_call}))',
                ),
            string = util.Typedef('string',  # implies null terminated string
                c_type   = 'char',    # XXX - char *
                cpp_type = 'std::string',
                cpp_to_c = '{var}.c_str()',  # . or ->
                c_fortran  = 'character(kind=C_CHAR)',
                f_type     = 'character(*)',
                fortran_to_c = 'trim({var}) // C_NULL_CHAR',
#                f_module = dict(iso_c_binding = [ 'C_NULL_CHAR' ]),
                f_module = dict(iso_c_binding=None),
                f_return_code = '{F_result} = fstr({F_C_name}({F_arg_c_call}))',
                PY_format = 's',
                base = 'string',
                ),
            )
        def_types['std::string'] = def_types['string']
        if 'typedef' in node and \
                node['typedef'] is not None:
            if not isinstance(node['typedef'], dict):
                raise TypeError("typedef must be a dictionary")
            for key, value in node['typedef'].items():
                if not isinstance(value, dict):
                    raise TypeError("typedef '%s' must be a dictionary" % key)
                if key in def_types:
                    def_types[key].update(value)
                else:
                    def_types[key] = util.Typedef(key, **value)

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
        fmt_class = self.push_fmt(node)
        fmt_class.cpp_class = name
        fmt_class.lower_class = name.lower()
        fmt_class.upper_class = name.upper()
        if 'C_prefix' in options:
            fmt_class.C_prefix = options.C_prefix

        if options.F_module_per_class:
            util.eval_template(options, fmt_class,
                               'F_module_name', '{lower_class}_mod')
            util.eval_template(options, fmt_class,
                               'F_impl_filename', 'wrapf{cpp_class}.f')   # XXX lower_class

        # Only one file per class for C.
        util.eval_template(options, fmt_class,
                           'C_header_filename', 'wrap{cpp_class}.h')
        util.eval_template(options, fmt_class,
                           'C_impl_filename', 'wrap{cpp_class}.cpp')

        # create typedef for each class before generating code
        # this allows classes to reference each other
        if name not in self.typedef:
#            unname = util.un_camel(name)
            unname = name.lower()
            cname = node['options'].C_prefix + unname
            self.typedef[name] = util.Typedef(
                name,
                cpp_type = name,
#                cpp_to_c = 'static_cast<void *>({var})',
                c_type = cname,
                c_to_cpp = 'static_cast<%s{ptr}>({var})' % name,
                c_fortran = 'type(C_PTR)',
                f_type = 'type(%s)' % unname,
                fortran_derived = unname,
                fortran_to_c = '{var}%{F_derived_member}',
                # XXX module name may not conflict with type name
                f_module = {fmt_class.F_module_name:[unname]},

                # return from C function
#                f_c_return_decl = 'type(CPTR)' % unname,
                f_return_code = '{F_result}%{F_derived_member} = {F_C_name}({F_arg_c_call})',

                # allow forward declarations to avoid recursive headers
                forward = name,
                base = 'wrapped',
                )

        typedef = self.typedef[name]
        fmt_class.C_type_name = typedef.c_type

        overloaded_methods = {}
        methods = node.setdefault('methods', [])
        for method in methods:
            if not isinstance(method, dict):
                raise TypeError("classes[n]['methods'] must be a dictionary")
            self.check_function(method)
            overloaded_methods.setdefault(method['result']['name'], []).append(method)

        # look for function overload and compute method_suffix
        for mname, methods in overloaded_methods.items():
            if len(methods) > 1:
                for i, method in enumerate(methods):
#                    method['fmt'].overloaded = True
                    if 'method_suffix' not in method:
                        method['fmt'].method_suffix =  '_%d' % i

        self.pop_fmt()
        self.pop_options()

    def check_function(self, node):
        self.push_options(node)
        fmt_func = self.push_fmt(node)

#        func = util.FunctionNode()
#        func.update(node)
#        func.dump()

        node['qualifiers'] = {}

        if 'decl' in node:
            # parse decl and add to dictionary
            values = parse_decl.check_decl(node['decl'])
            util.update(node, values)  # recursive update
        if 'method_suffix' in node and node['method_suffix'] is None:
            # YAML turns blanks strings to None
            node['method_suffix'] = ''
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
        if 'method_suffix' in node:
            fmt_func.method_suffix =   node['method_suffix']

        # docs
        self.pop_fmt()
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
    parser.add_argument('--indir', default='',
                        help='directory for input source files')
    parser.add_argument('--outdir', default='',
                        help='directory for output files')
    parser.add_argument('--logdir', default='',
                        help='directory for log files')
    parser.add_argument('filename', nargs='*',
                        help='Input file to process')

    args = parser.parse_args()

    # check command line options
    if len(args.filename) == 0:
        raise SystemExit("Must give at least one input file")
    if args.indir and not os.path.isdir(args.indir):
        raise SystemExit("indir %s does not exist" % args.indir)
    if args.outdir and not os.path.isdir(args.outdir):
        raise SystemExit("outdir %s does not exist" % args.outdir)
    if args.logdir and not os.path.isdir(args.logdir):
        raise SystemExit("logdir %s does not exist" % args.logdir)

    basename = os.path.splitext(os.path.basename(args.filename[0]))[0]
    logpath = os.path.join(args.logdir, basename + '.log')
    log = open(logpath, 'w')

    # pass around a configuration object
    config = Config()
    config.binary_dir = args.outdir
    config.source_dir = args.indir
    config.log = log

    # accumulated input
    all = {}
    splicers = dict(c={}, f={}, py={})

    for filename in args.filename:
        root, ext = os.path.splitext(filename)
        if ext == '.yaml':
            fp = open(filename, 'r')
            d = yaml.load(fp.read())
            fp.close()
            all.update(d)
#            util.update(all, d)  # recursive update
        elif ext == '.json':
            raise NotImplemented("Can not deal with json input for now")
        else:
            # process splicer file on command line
            splicer.get_splicer_based_on_suffix(filename, splicers)

#    print(all)


    Schema(all).check_schema()

    if 'splicer' in all:
        # read splicer files defined in input yaml file
        for suffix, names in all['splicer'].items():
            subsplicer = splicers.setdefault(suffix, {})
            for name in names:
                fullname = os.path.join(config.source_dir, name)
                splicer.get_splicers(fullname, subsplicer)


    wrapc.Wrapc(all, config, splicers['c']).wrap_library()

    wrapf.Wrapf(all, config, splicers['f']).wrap_library()

    wrapp.Wrapp(all, config, splicers['py']).wrap_library()

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

