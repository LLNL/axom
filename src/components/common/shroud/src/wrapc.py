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

wformat = util.wformat

class Wrapc(object):
    """Generate C bindings for C++ classes

    """
    def __init__(self, tree, log):
        self.tree = tree    # json tree
        self.log = log
        self.typedef = tree['typedef']

    def _clear_class(self):
        """Start a new class for output"""
        self.header_forward = {}          # forward declarations of C++ class as opaque C struct.
        self.header_typedef_include = {}  # include files required by typedefs
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
        t = []
        typedef = self.typedef.get(arg['type'], None)
        if typedef is None:
            raise RuntimeError("No such type %s" % arg['type'])
        if arg['attrs'].get('const', False):
            t.append('const')
        t.append(typedef[lang])
        if arg['attrs'].get('ptr', False):
            t.append('*')
        elif arg['attrs'].get('reference', False):
            if lang == 'cpp':
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
        typ = self._c_type(lang, arg)
        return typ + ' ' + ( name or arg['name'] )

    def wrap(self):
        for node in self.tree['classes']:
            self._clear_class()
            name = node['name']
            self.wrap_class(node)
            c_header = node['C_header_filename']
            c_impl   = node['C_impl_filename']
            self.write_header(node, c_header)
            self.write_impl(node, c_header, c_impl)

        for node in self.tree['functions']:
            self.wrap_function(node)

    def write_copyright(self, fp):
        for line in self.tree.get('copyright', []):
            if line:
                fp.write('// ' + line + '\n')
            else:
                fp.write('//\n')

    def write_header(self, node, fname):
        fp = open(fname, 'w')
        self.write_copyright(fp)

        guard = fname.replace(".", "_").upper()

        fp.writelines([
                '// %s\n' % fname,
                '// For C users and C++ implementation\n',
                '\n',
                '#ifndef %s\n' % guard,
                '#define %s\n' % guard,
                ])

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

    def write_impl(self, node, hname, fname):
        # node = class node
        namespace = node['options'].get('namespace',None)
        fp = open(fname, "w")
        self.write_copyright(fp)
        fp.write('// {}\n'.format(fname))
        fp.write('#define EXAMPLE_WRAPPER_IMPL\n')
        fp.write('#include "%s"\n' % hname)
        if 'cpp_header' in node['options']:
            fp.write('#include "%s"\n' % node['options']['cpp_header'])
        fp.write('\nextern "C" {\n')
        for name in namespace.split():
            fp.write('namespace %s {\n' % name)
        fp.writelines(self.impl)
        fp.write('\n')
        for name in namespace.split():
            fp.write('}  // namespace %s\n' % name)
        fp.write('}  // extern "C"\n')
        fp.close()
        self.log.write("Close %s\n" % fname)
        print("Wrote", fname)

    def wrap_function(self, node):
        """
        node - function node
        """
        self.log.write("function {1[decl]}\n".format(self, node))

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        typedef = self.typedef[name]
        cname = typedef['c']

        self.fmt_dict = dict(
            cpp_class = name,
            lower_class = name.lower(),
            upper_class = name.upper(),
            C_prefix = node['options']['C_prefix'],
            C_type_name = cname,
            )
        fmt_dict = self.fmt_dict

        # create a forward declaration for this type
        self.header_forward[cname] = True

        for method in node['methods']:
            self.wrap_method(node, method)

    def wrap_method(self, cls, node):
        """
        cls  - class node
        node - function node
        """
        if 'decl' in node:
            self.log.write("method {1[decl]}\n".format(self, node))
        else:
            self.log.write("method {1[result][name]}\n".format(self, node))
        # assume a C++ method

        # return type
        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)
        result_typedef = self.typedef[result_type]
        C_this = node['options']['C_this']
        is_const = result['attrs'].get('const', False)
        is_ctor  = result['attrs'].get('constructor', False)
        is_dtor  = result['attrs'].get('destructor', False)

        if 'c_header' in result_typedef:
            # include any dependent header in generated header
            self.header_typedef_include[result_typedef['c_header']] = True
        if 'forward' in result_typedef:
            # create forward references for other types being wrapped
            # i.e. This method returns a wrapped type
            self.header_forward[result_typedef['c']] = True

        fmt_dict = dict(
            method_name=result['name'],
            underscore_name=util.un_camel(result['name']),
            method_suffix = node.get('method_suffix', ''),

            const='const ' if is_const else '',
            this=C_this,
            cpp_this = C_this + 'obj',
            cpp_name = result['name'],
            rv_decl = self._c_decl('cpp', result, name='rv'),  # return value
            )
        fmt_dict.update(self.fmt_dict)

        arguments = []
        anames = []
        # object pointer
        arg_dict = dict(name=C_this,
                        type=cls['name'], 
                        attrs=dict(ptr=True,
                                   const=is_const))
        C_this_type = self._c_type('c', arg_dict)
        if not is_ctor:
            arg = self._c_decl('c', arg_dict)
            arguments.append(arg)

        for arg in node.get('args', []):
            arguments.append(self._c_decl('c', arg))
            arg_typedef = self.typedef[arg['type']]
            # convert C argument to C++
            anames.append(arg_typedef.get('c_to_cpp', '{var}').format(
                    var=arg['name'],
                    ptr=' *' if arg['attrs'].get('ptr', False) else ''))
            if 'c_header' in arg_typedef:
                # include any dependent header in generated header
                self.header_typedef_include[arg_typedef['c_header']] = True
            if 'forward' in arg_typedef:
                # create forward references for other types being wrapped
                # i.e. This argument is another wrapped type
                self.header_forward[arg_typedef['c']] = True
        fmt_dict['call_list'] = ', '.join(anames)

        if 'C_name' not in node:
            node['C_name'] = wformat(
                node['options']['C_name_method_template'],
                fmt_dict)

        if 'C_return_type' not in node:
            node['C_return_type'] = self._c_type('c', result)

        if 'C_arguments' not in node:
            node['C_arguments'] = ', '.join(arguments)

        if 'C_object' not in node:
            if is_ctor:
                template = '{const}{cpp_class} *{this}obj = new {cpp_class}({call_list});'
            else:
                template = '{const}{cpp_class} *{this}obj = static_cast<{const}{cpp_class} *>({this});'
            node['C_object'] = wformat(template, fmt_dict)

        if 'C_code' not in node:
            # generate the C body
            lines = []
            if is_ctor:
                lines.append('return (%s) %sobj;' % (C_this_type, C_this))
            elif is_dtor:
                lines.append('delete %sobj;' % C_this)
            elif result_type == 'void' and not result_is_ptr:
                line = wformat('{cpp_this}->{cpp_name}({call_list});',
                               fmt_dict)
                lines.append(line)
                lines.append('return;')
            else:
                line = wformat('{rv_decl} = {cpp_this}->{cpp_name}({call_list});',
                               fmt_dict)
                lines.append(line)

                ret = result_typedef.get('cpp_to_c', '{var}')
                line = 'return ' + ret.format(var='rv') + ';'
                lines.append(line)

            node['C_code'] = "\n".join(lines)

        self.header_proto_c.append(wformat('\n{C_return_type} {C_name}({C_arguments});\n',
                                           node))

        self.impl.append(wformat("""
{C_return_type} {C_name}({C_arguments})
{{
{C_object}
// splicer begin
{C_code}
// splicer end
}}
""", node))
