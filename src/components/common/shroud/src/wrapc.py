#!/bin/env python3
"""
Generate C bindings for C++ classes

#ifdef __cplusplus
typedef void {C_type_name}
#else
struct s_{C_type_name};
typedef struct s_{C_type_name} {C_type_name};
#endif


"""
from __future__ import print_function

import util

import os

wformat = util.wformat

class Wrapc(object):
    """Generate C bindings for C++ classes

    """
    def __init__(self, tree, config, splicers):
        self.tree = tree    # json tree
        self.config = config
        self.log = config.log
        self.typedef = tree['typedef']
        self._init_splicer(splicers)
        self.comment = '//'

    def _begin_output_file(self):
        """Start a new class for output"""
        self.header_forward = {}          # forward declarations of C++ class as opaque C struct.
        self.header_typedef_include = {}  # include files required by typedefs
        self.header_impl_include = {}     # headers needed by implementation, i.e. helper functions
        self.header_proto_c = []
        self.impl = []

#####

    def _init_splicer(self, splicers):
        self.splicers = splicers
        self.splicer_stack = [ splicers ]
        self.splicer_names = [ ]
        self.splicer_path = ''

    def _push_splicer(self, name, out):
        level = self.splicer_stack[-1].setdefault(name, {})
        self.splicer_stack.append(level)
        self.splicer_names.append(name)
        self.splicer_path = '.'.join(self.splicer_names) + '.'
#        out.append('// splicer push %s' % name)

#X changes for push/pop instead of full paths
    def _pop_splicer(self, name, out):
        # XXX maybe use name for error checking, must pop in reverse order
        self.splicer_stack.pop()
        self.splicer_names.pop()
        if self.splicer_names:
            self.splicer_path = '.'.join(self.splicer_names) + '.'
        else:
            self.splicer_path = ''
#X        out.append('// splicer pop %s' % name)

    def _create_splicer(self, name, out, default=None):
        # The prefix is needed when two different sets of output are being create
        # and they are not in sync.
        # Creating methods and derived types together.
#X        out.append('// splicer begin %s' % name)
        out.append('// splicer begin %s%s' % (self.splicer_path, name))
        if default:
            out.extend(default)
        else:
            out.extend(self.splicer_stack[-1].get(name, []))
#X        out.append('// splicer end %s' % name)
        out.append('// splicer end %s%s' % (self.splicer_path, name))

#####

    def _c_type(self, lang, arg):
        """
        Return the C type.
        pass-by-value default

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
        if arg['attrs'].get('const', False):
            t.append('const')
        t.append(getattr(typedef, lang))
        if arg['attrs'].get('ptr', False):
            t.append('*')
        elif arg['attrs'].get('reference', False):
            if lang == 'cpp_type':
                t.append('&')
            else:
                t.append('*')
        return ' '.join(t)

    def _c_decl(self, lang, arg, name=None):
        """
        Return the C declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
#        if lang not in [ 'c_type', 'cpp_type' ]:
#            raise RuntimeError
        typ = self._c_type(lang, arg)
        return typ + ' ' + ( name or arg['name'] )

    def wrap_library(self):
        options = self.tree['options']
        fmt_library = self.tree['fmt']
        fmt_library.C_this = options.get('C_this', 'self')
        fmt_library.C_const = ''

        self._push_splicer('class', [])
        for node in self.tree['classes']:
            self.write_file(node, self.wrap_class, True)
        self._pop_splicer('class', [])

        if self.tree['functions']:
            self.write_file(self.tree, self.wrap_functions, False)

    def write_file(self, node, worker, cls):
        """Write a file for the library and it's functions or
        a class and it's methods.
        """
        fmt = node['fmt']
        self._begin_output_file()
        worker(node)
        c_header = fmt.C_header_filename
        c_impl   = fmt.C_impl_filename
        self.write_header(node, c_header, cls)
        self.write_impl(node, c_header, c_impl, cls)

    def wrap_functions(self, tree):
        # worker function for write_file
        self._push_splicer('function', [])
        for node in tree['functions']:
            self.wrap_method(None, node)
        self._pop_splicer('function', [])

    def write_header(self, node, fname, cls=False):
        guard = fname.replace(".", "_").upper()

        output = []
        output.extend([
                '// %s' % fname,
                '// For C users and C++ implementation',
                '',
                '#ifndef %s' % guard,
                '#define %s' % guard,
                ])
        if cls:
            self._push_splicer('class', output)

        # headers required by typedefs
        if self.header_typedef_include:
#            output.append('// header_typedef_include')
            output.append('')
            headers = self.header_typedef_include.keys()
            headers.sort()
            for header in headers:
                output.append('#include "%s"' % header)

        output.extend([
                '',
                '#ifdef __cplusplus',
                'extern "C" {',
                '#endif',
                '',
                '// declaration of wrapped types',
                '#ifdef EXAMPLE_WRAPPER_IMPL',
                ])
        names = self.header_forward.keys()
        names.sort()
        for name in names:
            output.append('typedef void {C_type_name};'.format(C_type_name=name))
        output.append('#else')
        for name in names:
            output.append('struct s_{C_type_name};\ntypedef struct s_{C_type_name} {C_type_name};'.
                     format(C_type_name=name))
        output.append('#endif')
        output.extend(self.header_proto_c);
        if cls:
            self._pop_splicer('class', output)
        output.extend([
                '',
                '#ifdef __cplusplus',
                '}',
                '#endif',
                '',
                '#endif  // %s' % guard
                ])

        self.write_output_file(fname, self.config.binary_dir, output)

    def write_impl(self, node, hname, fname, cls=False):
        # node = class node
        options = node['options']
        namespace = options.namespace

        output = []
        output.append('// ' + fname)
        output.append('#define EXAMPLE_WRAPPER_IMPL')

        output.append('#include "%s"' % hname)
        if node['options'].cpp_header:
            for include in node['options'].cpp_header.split():
                self.header_impl_include[include] = True
        # headers required by implementation
        if self.header_impl_include:
            headers = self.header_impl_include.keys()
            headers.sort()
            for header in headers:
                output.append('#include "%s"' % header)

        output.append('\nextern "C" {')
        for name in namespace.split():
            output.append('namespace %s {' % name)
        output.extend(self.impl)
        output.append('')
        for name in namespace.split():
            output.append('}  // namespace %s' % name)

        output.append('}  // extern "C"')

        self.write_output_file(fname, self.config.binary_dir, output)

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        typedef = self.typedef[name]
        cname = typedef.c_type

        self._push_splicer(name, [])
#        fmt_class = node['fmt']
#        fmt_class.update(dict(
#                ))

        # create a forward declaration for this type
        self.header_forward[cname] = True

        self._push_splicer('method', [])
        for method in node['methods']:
            self.wrap_method(node, method)
        self._pop_splicer('method', [])

        self._pop_splicer(name, [])

    def wrap_method(self, cls, node):
        """
        cls  - class node or None for functions
        node - function/method node
        """
        if cls:
            cls_function = 'method'
        else:
            cls_function = 'function'
        if 'decl' in node:
            self.log.write("{0} {1[decl]}\n".format(cls_function, node))
        else:
            self.log.write("{0} {1[result][name]}\n".format(cls_function, node))

        fmt_func = node['fmt']

        # return type
        options = node['options']
        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)

        if node.get('return_this', False):
            result_type = 'void'
            result_is_ptr = False

        result_typedef = self.typedef[result_type]
        is_const = result['attrs'].get('const', False)
        is_ctor  = result['attrs'].get('constructor', False)
        is_dtor  = result['attrs'].get('destructor', False)

        if 'C_this' in options:
            fmt_func.C_this = options.C_this
        C_this = fmt_func.C_this

        if result_typedef.c_header:
            # include any dependent header in generated header
            self.header_typedef_include[result_typedef.c_header] = True
        if result_typedef.cpp_header:
            # include any dependent header in generated source
            self.header_impl_include[result_typedef.cpp_header] = True
        if result_typedef.forward:
            # create forward references for other types being wrapped
            # i.e. This method returns a wrapped type
            self.header_forward[result_typedef.c_type] = True

        if is_const:
            fmt_func.C_const = 'const '
        fmt_func.CPP_this = C_this + 'obj'
        fmt_func.CPP_name = result['name']
        fmt_func.rv_decl = self._c_decl('cpp_type', result, name='rv')  # return value

        arguments = []
        anames = []
        if cls:
            # object pointer
            fmt_func.CPP_this_call = fmt_func.CPP_this + '->'  # call method syntax
            arg_dict = dict(name=C_this,
                            type=cls['name'], 
                            attrs=dict(ptr=True,
                                       const=is_const))
            C_this_type = self._c_type('c_type', arg_dict)
            if not is_ctor:
                arg = self._c_decl('c_type', arg_dict)
                arguments.append(arg)
        else:
            fmt_func.CPP_this_call = ''  # call method syntax


        for arg in node.get('args', []):
            arguments.append(self._c_decl('c_type', arg))
            arg_typedef = self.typedef[arg['type']]
            # convert C argument to C++
            anames.append(arg_typedef.c_to_cpp.format(
                    var=arg['name'],
                    ptr=' *' if arg['attrs'].get('ptr', False) else ''))
            if arg_typedef.c_header:
                # include any dependent header in generated header
                self.header_typedef_include[arg_typedef.c_header] = True
            if arg_typedef.cpp_header:
                # include any dependent header in generated source
                self.header_impl_include[arg_typedef.cpp_header] = True
            if arg_typedef.forward:
                # create forward references for other types being wrapped
                # i.e. This argument is another wrapped type
                self.header_forward[arg_typedef.c_type] = True
        fmt_func.C_call_list = ', '.join(anames)

        fmt_func.C_arguments = options.get('C_arguments', ', '.join(arguments))

        if node.get('return_this', False):
            fmt_func.C_return_type = 'void'
        else:
            fmt_func.C_return_type = options.get('C_return_type', self._c_type('c_type', result))

        if cls:
            util.eval_template(options, fmt_func, 'C_name',
                               '{C_prefix}{lower_class}_{underscore_name}{method_suffix}')
            if 'C_object' in options:
                fmt_func.C_object = options.C_object
            else:
                if is_ctor:
                    template = '{C_const}{cpp_class} *{C_this}obj = new {cpp_class}({C_call_list});'
                else:
                    template = '{C_const}{cpp_class} *{C_this}obj = static_cast<{C_const}{cpp_class} *>({C_this});'
                fmt_func.C_object = wformat(template, fmt_func)
        else:
            util.eval_template(options, fmt_func, 'C_name',
                               '{C_prefix}{underscore_name}{method_suffix}')
            fmt_func.C_object = ''

        # body of function
        splicer_code = self.splicer_stack[-1].get(fmt_func.method_name, None)
        if 'C_code' in options:
            C_code = [ options.C_code ]
        elif splicer_code:
            C_code = splicer_code
        else:
            # generate the C body
            C_code = []
            if is_ctor:
                C_code.append('return (%s) %sobj;' % (C_this_type, C_this))
            elif is_dtor:
                C_code.append('delete %sobj;' % C_this)
            elif result_type == 'void' and not result_is_ptr:
                line = wformat('{CPP_this_call}{CPP_name}({C_call_list});',
                               fmt_func)
                C_code.append(line)
                C_code.append('return;')
            else:
                line = wformat('{rv_decl} = {CPP_this_call}{CPP_name}({C_call_list});',
                               fmt_func)
                C_code.append(line)

                ret = result_typedef.cpp_to_c
                line = 'return ' + ret.format(var='rv') + ';'
                C_code.append(line)

        self.header_proto_c.append('')
        self.header_proto_c.append(wformat('{C_return_type} {C_name}({C_arguments});',
                                           fmt_func))

        impl = self.impl
        impl.append('')
        impl.append(wformat('{C_return_type} {C_name}({C_arguments})', fmt_func))
        impl.append('{')
        if cls:
            impl.append(fmt_func.C_object )
        self._create_splicer(fmt_func.method_name, impl, C_code)
        impl.append('}')

#####

    def write_output_file(self, fname, directory, output):
        fp = open(os.path.join(directory, fname), 'w')
        fp.write('%s %s\n' % (self.comment, fname))
        self.write_copyright(fp)
        self.indent = 0
        self.write_lines(fp, output)
        fp.close()
        self.log.write("Close %s\n" % fname)
        print("Wrote", fname)

    def write_copyright(self, fp):
        for line in self.tree.get('copyright', []):
            if line:
                fp.write(self.comment + ' ' + line + '\n')
            else:
                fp.write(self.comment + '\n')

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
