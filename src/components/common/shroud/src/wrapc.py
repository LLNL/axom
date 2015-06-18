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
        self.splicers = splicers
        self.log = config.log
        self.typedef = tree['typedef']

    def _begin_output_file(self):
        """Start a new class for output"""
        self.header_forward = {}          # forward declarations of C++ class as opaque C struct.
        self.header_typedef_include = {}  # include files required by typedefs
        self.header_impl_include = {}     # headers needed by implementation, i.e. helper functions
        self.header_proto_c = []
        self.impl = []

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

        for node in self.tree['classes']:
            self.write_file(node, self.wrap_class, True)

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
        for node in tree['functions']:
            self.wrap_method(None, node)

    def write_copyright(self, fp):
        for line in self.tree.get('copyright', []):
            if line:
                fp.write('// ' + line + '\n')
            else:
                fp.write('//\n')

    def write_header(self, node, fname, cls=False):
        fp = open(os.path.join(self.config.binary_dir, fname), 'w')
        self.write_copyright(fp)

        guard = fname.replace(".", "_").upper()

        fp.writelines([
                '// %s\n' % fname,
                '// For C users and C++ implementation\n',
                '\n',
                '#ifndef %s\n' % guard,
                '#define %s\n' % guard,
                ])
        if cls:
            fp.write('// splicer push.class.%s\n' % node['name'])

        # headers required by typedefs
        if self.header_typedef_include:
#            fp.write('// header_typedef_include\n')
            fp.write('\n')
            headers = self.header_typedef_include.keys()
            headers.sort()
            for header in headers:
                fp.write('#include "%s"\n' % header)

        fp.writelines([
                '\n',
                '#ifdef __cplusplus\n',
                'extern "C" {\n',
                '#endif\n',
                '\n',
                '// declaration of wrapped types\n',
                '#ifdef EXAMPLE_WRAPPER_IMPL\n',
                ])
        names = self.header_forward.keys()
        names.sort()
        for name in names:
            fp.write('typedef void {C_type_name};\n'.format(C_type_name=name))
        fp.write('#else\n')
        for name in names:
            fp.write('struct s_{C_type_name};\ntypedef struct s_{C_type_name} {C_type_name};\n'.
                     format(C_type_name=name))
        fp.write('#endif\n')
        fp.writelines(self.header_proto_c);
        if cls:
            fp.write('\n// splicer pop.class.%s\n' % node['name'])
        fp.writelines([
                '\n',
                '#ifdef __cplusplus\n',
                '}\n',
                '#endif\n',
                '\n',
                '#endif  // %s\n' % guard
                ])
        fp.close()
        self.log.write("Close %s\n" % fname)
        print("Wrote", fname)

    def write_impl(self, node, hname, fname, cls=False):
        # node = class node
        options = node['options']
        namespace = options.namespace
        fp = open(os.path.join(self.config.binary_dir, fname), "w")
        self.write_copyright(fp)
        fp.write('// {}\n'.format(fname))
        fp.write('#define EXAMPLE_WRAPPER_IMPL\n')

        fp.write('#include "%s"\n' % hname)
        if node['options'].cpp_header:
            for include in node['options'].cpp_header.split():
                self.header_impl_include[include] = True
        # headers required by implementation
        if self.header_impl_include:
            headers = self.header_impl_include.keys()
            headers.sort()
            for header in headers:
                fp.write('#include "%s"\n' % header)

        fp.write('\nextern "C" {\n')
        for name in namespace.split():
            fp.write('namespace %s {\n' % name)
        if cls:
            fp.write('// splicer push class.%s.method\n' % node['name'])
        fp.writelines(self.impl)
        if cls:
            fp.write('\n// splicer pop.class.%s method\n' % node['name'])
        fp.write('\n')
        for name in namespace.split():
            fp.write('}  // namespace %s\n' % name)
        fp.write('}  // extern "C"\n')
        fp.close()
        self.log.write("Close %s\n" % fname)
        print("Wrote", fname)

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        typedef = self.typedef[name]
        cname = typedef.c_type

#        fmt_class = node['fmt']
#        fmt_class.update(dict(
#                ))

        # create a forward declaration for this type
        self.header_forward[cname] = True

        for method in node['methods']:
            self.wrap_method(node, method)

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

        if 'C_code' in options:
            fmt_func.C_code = options.C_code
        else:
            # generate the C body
            lines = []
            if is_ctor:
                lines.append('return (%s) %sobj;' % (C_this_type, C_this))
            elif is_dtor:
                lines.append('delete %sobj;' % C_this)
            elif result_type == 'void' and not result_is_ptr:
                line = wformat('{CPP_this_call}{CPP_name}({C_call_list});',
                               fmt_func)
                lines.append(line)
                lines.append('return;')
            else:
                line = wformat('{rv_decl} = {CPP_this_call}{CPP_name}({C_call_list});',
                               fmt_func)
                lines.append(line)

                ret = result_typedef.cpp_to_c
                line = 'return ' + ret.format(var='rv') + ';'
                lines.append(line)

            fmt_func.C_code = "\n".join(lines)

        self.header_proto_c.append(wformat('\n{C_return_type} {C_name}({C_arguments});\n',
                                           fmt_func))

        self.impl.append(wformat("""
{C_return_type} {C_name}({C_arguments})
{{
{C_object}
// splicer begin {C_name}
{C_code}
// splicer end {C_name}
}}
""", fmt_func))
