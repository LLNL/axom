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
import wrapc
import wrapf

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

            C_this='self',   # object argument name
            F_this='obj',    # object argument name
            F_result='rv',   # function result
            F_module_per_class=True,

            C_name_method_template='{C_prefix}{lower_class}_{underscore_name}{method_suffix}',
#            F_name_method_template='{underscore_name}{method_suffix}',
#            F_name_generic_template='{underscore_name}',
#            F_name_impl_template  ='{lower_class}_{underscore_name}{method_suffix}',
            )
        if 'options' in node:
            def_options.update(node['options'])
        self.options_stack = [ def_options ]

        fmt2_library = node['fmt'] = util.Options(None)
        fmt2_library.library       = def_options.get('library', 'default_library')
        fmt2_library.lower_library = fmt2_library.library.lower()
        fmt2_library.method_suffix = ''   # assume no suffix
        self.fmt_stack.append(fmt2_library)

        # default some options based on other options
        if not def_options.F_module_per_class:
            util.eval_template(def_options, fmt2_library,
                               'F_module_name', '{lower_library}_mod')
            util.eval_template(def_options, fmt2_library,
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
                ),
            long   = util.Typedef('long',
                c_type    = 'long',
                cpp_type  = 'long',
                c_fortran = 'integer(C_LONG)',
                f_type    = 'integer(C_LONG)',
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
                fortran_to_c  = 'logical2bool({var})',
                f_type    = 'logical',
                f_return_code = '{F_result} = bool2logical({F_C_name}({arg_c_call}))',
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
                f_return_code = '{F_result} = fstr({F_C_name}({arg_c_call}))',
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
                    def_types.update(value)
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
        fmt2_class = self.push_fmt(node)
        fmt_class = dict(
            cpp_class = name,
            lower_class = name.lower(),
            upper_class = name.upper(),
            C_prefix = options.get('C_prefix', ''),
            )
        fmt2_class.update(fmt_class)

        if options.F_module_per_class:
            util.eval_template(options, fmt2_class,
                               'F_module_name', '{lower_class}_mod')
            util.eval_template(options, fmt2_class,
                               'F_impl_filename', 'wrapf{cpp_class}.f')   # XXX lower_class

        # Only one file per class for C.
        util.eval_template(options, fmt2_class,
                           'C_header_filename', 'wrap{cpp_class}.h')
        util.eval_template(options, fmt2_class,
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
                c_type = cname,
                c_to_cpp = 'static_cast<%s{ptr}>({var})' % name,
                c_fortran = 'type(C_PTR)',
                f_type = 'type(%s)' % unname,
                fortran_derived = unname,
                fortran_to_c = '{var}%obj',
                # XXX module name may not conflict with type name
                f_module = {fmt2_class.F_module_name:[unname]},

                # return from C function
#                f_c_return_decl = 'type(CPTR)' % unname,
                f_return_code = '{F_result}%{F_this} = {F_C_name}({arg_c_call})',

                # allow forward declarations to avoid recursive headers
                forward = name,
                base = 'wrapped',
                )

        typedef = self.typedef[name]
        fmt2_class.C_type_name = typedef.c_type

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
                    if 'method_suffix' not in method:
                        method['method_suffix'] = '_%d' % i
                        method['fmt'].method_suffix =  '_%d' % i

        self.pop_fmt()
        self.pop_options()

    def check_function(self, node):
        self.push_options(node)
        fmt2_func = self.push_fmt(node)

#        func = util.FunctionNode()
#        func.update(node)
#        func.dump()

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

        result = node['result']

        fmt2_func.method_name =     result['name']
        fmt2_func.underscore_name = util.un_camel(result['name'])
        fmt2_func.method_suffix =   node.get('method_suffix', '')

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

    for filename in args.filename:
        root, ext = os.path.splitext(filename)
        if ext == '.yaml':
            fp = open(filename, 'r')
            d = yaml.load(fp.read())
            fp.close()
            all.update(d)
        else:
            raise SystemExit("File must end in .yaml for now")

#    print(all)



    Schema(all).check_schema()

    wrapc.Wrapc(all, config).wrap_library()

    wrapf.Wrapf(all, config).wrap_library()

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

