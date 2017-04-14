#!/bin/env python3
"""
Generate C bindings for C++ classes

"""
from __future__ import print_function
from __future__ import absolute_import

import os

from . import whelpers
from . import util
from .util import append_format

wformat = util.wformat


class Wrapc(util.WrapperMixin):
    """Generate C bindings for C++ classes

    """
    def __init__(self, tree, config, splicers):
        self.tree = tree    # json tree
        self.patterns = tree['patterns']
        self.config = config
        self.log = config.log
        self.typedef = tree['types']
        self._init_splicer(splicers)
        self.comment = '//'
        self.doxygen_begin = '/**'
        self.doxygen_cont = ' *'
        self.doxygen_end = ' */'

    def _begin_output_file(self):
        """Start a new class for output"""
        # forward declarations of C++ class as opaque C struct.
        self.header_forward = {}
        # include files required by typedefs
        self.header_typedef_include = {}
        # headers needed by implementation, i.e. helper functions
        self.header_impl_include = {}
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
        typ = getattr(typedef, lang)
        if typ is None:
            raise RuntimeError(
                "Type {} has no value for {}".format(arg['type'], lang))
        t.append(typ)
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
        return typ + ' ' + (name or arg['name'])

    def wrap_library(self):
        fmt_library = self.tree['fmt']

        self._push_splicer('class')
        for node in self.tree['classes']:
            self._push_splicer(node['name'])
            self.write_file(self.tree, node)
            self._pop_splicer(node['name'])
        self._pop_splicer('class')

        if self.tree['functions']:
            self.write_file(self.tree, None)

        self.write_helper_files()

    def write_file(self, library, cls):
        """Write a file for the library and its functions or
        a class and its methods.
        """
        node = cls or library
        fmt = node['fmt']
        self._begin_output_file()
        if cls:
            self.wrap_class(cls)
        else:
            self.wrap_functions(library)
        c_header = fmt.C_header_filename
        c_impl = fmt.C_impl_filename
        self.write_header(library, cls, c_header)
        self.write_impl(library, cls, c_header, c_impl)

    def wrap_functions(self, tree):
        # worker function for write_file
        self._push_splicer('function')
        for node in tree['functions']:
            self.wrap_function(None, node)
        self._pop_splicer('function')

    def write_header(self, library, cls, fname):
        """ Write header file for a library node or a class node.
        """
        guard = fname.replace(".", "_").upper()
        node = cls or library
        options = node['options']

        output = []

        if options.doxygen:
            self.write_doxygen_file(output, fname, library, cls)

        output.extend([
                '// For C users and C++ implementation',
                '',
                '#ifndef %s' % guard,
                '#define %s' % guard,
                ])
        # headers required by typedefs
        if self.header_typedef_include:
            # output.append('// header_typedef_include')
            output.append('')
            headers = self.header_typedef_include.keys()
            self.write_headers(headers, output)

        output.extend([
                '',
                '#ifdef __cplusplus',
                'extern "C" {',
                '#endif',
                '',
                '// declaration of wrapped types',
                ])
        names = sorted(self.header_forward.keys())
        for name in names:
            output.append(
                'struct s_{C_type_name};\n'
                'typedef struct s_{C_type_name} {C_type_name};'.
                format(C_type_name=name))
        output.append('')
        self._create_splicer('C_definition', output)
        output.extend(self.header_proto_c)
        output.extend([
                '',
                '#ifdef __cplusplus',
                '}',
                '#endif',
                '',
                '#endif  // %s' % guard
                ])

        self.config.cfiles.append(
            os.path.join(self.config.c_fortran_dir, fname))
        self.write_output_file(fname, self.config.c_fortran_dir, output)

    def write_impl(self, library, cls, hname, fname):
        """Write implementation
        """
        node = cls or library
        options = node['options']

        output = []
        output.append('// ' + fname)

        output.append('#include "%s"' % hname)

        # Use headers from class if they exist or else library
        if cls and cls['cpp_header']:
            for include in cls['cpp_header'].split():
                self.header_impl_include[include] = True
        else:
            for include in library['cpp_header'].split():
                self.header_impl_include[include] = True

        # headers required by implementation
        if self.header_impl_include:
            headers = self.header_impl_include.keys()
            self.write_headers(headers, output)

        output.append('\nextern "C" {')
        self.namespace(library, cls, 'begin', output)
        output.extend(self.impl)
        output.append('')
        self._create_splicer('additional_functions', output)
        output.append('')
        self.namespace(library, cls, 'end', output)

        output.append('}  // extern "C"')

        self.config.cfiles.append(
            os.path.join(self.config.c_fortran_dir, fname))
        self.write_output_file(fname, self.config.c_fortran_dir, output)

    def write_headers(self, headers, output):
        for header in sorted(headers):
            if header[0] == '<':
                output.append('#include %s' % header)
            else:
                output.append('#include "%s"' % header)

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        typedef = self.typedef[name]
        cname = typedef.c_type

        fmt_class = node['fmt']
        # call method syntax
        fmt_class.CPP_this_call = fmt_class.CPP_this + '->'
#        fmt_class.update(dict(
#                ))

        # create a forward declaration for this type
        self.header_forward[cname] = True

        self._push_splicer('method')
        for method in node['methods']:
            self.wrap_function(node, method)
        self._pop_splicer('method')

    def wrap_function(self, cls, node):
        """
        Wrap a C++ function with C
        cls  - class node or None for functions
        node - function/method node
        """
        options = node['options']
        if not options.wrap_c:
            return

        if cls:
            cls_function = 'method'
        else:
            cls_function = 'function'
        self.log.write("C {0} {1[_decl]}\n".format(cls_function, node))

        fmt_func = node['fmt']

        # Look for C++ routine to wrap
        # Usually the same node unless it is bufferified
        CPP_node = node
        generated = []
        if '_generated' in CPP_node:
            generated.append(CPP_node['_generated'])
        while '_PTR_C_CPP_index' in CPP_node:
            CPP_node = self.tree['function_index'][
                CPP_node['_PTR_C_CPP_index']]
            if '_generated' in CPP_node:
                generated.append(CPP_node['_generated'])
        CPP_result = CPP_node['result']
        CPP_result_type = CPP_result['type']
        CPP_subprogram = CPP_node['_subprogram']

        # C return type
        result = node['result']
        result_type = result['type']
        subprogram = node['_subprogram']

        # C++ functions which return 'this',
        # are easier to call from Fortran if they are subroutines.
        # There is no way to chain in Fortran:  obj->doA()->doB();
        if node.get('return_this', False):
            CPP_result_type = 'void'
            CPP_subprogram = 'subroutine'

        result_typedef = self.typedef[result_type]
        result_is_const = result['attrs'].get('const', False)
        is_ctor = node['attrs'].get('constructor', False)
        is_dtor = node['attrs'].get('destructor', False)
        is_const = node['attrs'].get('const', False)

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

        if result_is_const:
            fmt_func.c_const = 'const '
        else:
            fmt_func.c_const = ''
        # return value
        fmt_func.C_rv_decl = self._c_decl('cpp_type', CPP_result, name=fmt_func.C_result)

        proto_list = []
        call_list = []
        if cls:
            # object pointer
            arg_dict = dict(name=fmt_func.C_this,
                            type=cls['name'],
                            attrs=dict(ptr=True,
                                       const=is_const))
            C_this_type = self._c_type('c_type', arg_dict)
            if not is_ctor:
                arg = self._c_decl('c_type', arg_dict)
                proto_list.append(arg)

        # indicate which argument contains function result, usually none
        result_arg = None
        pre_call = []      # list of temporary variable declarations
        post_call = []

        if cls and not is_ctor:
            if is_const:
                fmt_func.c_const = 'const '
            else:
                fmt_func.c_const = ''
            fmt_func.c_ptr = ' *'
            fmt_func.c_var = fmt_func.C_this
            # LHS is class' cpp_to_c
            cls_typedef = self.typedef[cls['name']]
            append_format(pre_call, 
                          '{c_const}{cpp_class} *{CPP_this} = ' +
                          cls_typedef.c_to_cpp + ';', fmt_func)

#    c_var     - argument to C function  (wrapper function)
#    c_var_len - variable with trimmed length of c_var
#    c_var_num - variable with length of c_var
#    cpp_var   - argument to C++ function  (wrapped function).
#                Usually same as c_var but may be a new local variable
#                or the funtion result variable.

        for arg in node['args']:
            fmt_arg = arg.setdefault('fmtc', util.Options(fmt_func))
            c_attrs = arg['attrs']
            arg_typedef = self.typedef[arg['type']]
            c_statements = arg_typedef.c_statements
            fmt_arg.c_var = arg['name']
            fmt_arg.c_var_len = 'L' + fmt_arg.c_var
            fmt_arg.c_var_num = 'N' + fmt_arg.c_var
            if arg['attrs'].get('const', False):
                fmt_arg.c_const = 'const '
            else:
                fmt_arg.c_const = ''
            if arg['attrs'].get('ptr', False):
                fmt_arg.c_ptr = ' *'
            else:
                fmt_arg.c_ptr = ''
            fmt_arg.cpp_type = arg_typedef.cpp_type

            if c_attrs.get('_is_result', False):
                arg_call = False
                result_arg = arg
                slist = ['result']
            else:
                arg_call = arg
                fmt_arg.cpp_var = fmt_arg.c_var      # name in c++ call.
                slist = []
                if c_attrs['intent'] in ['inout', 'in']:
                    slist.append('intent_in')
                if c_attrs['intent'] in ['inout', 'out']:
                    slist.append('intent_out')

            proto_list.append(self._c_decl('c_type', arg))

            len_trim = arg['attrs'].get('len_trim', False)
            if len_trim:
                fmt_arg.c_var_trim = len_trim
                append_format(proto_list, 'int {c_var_trim}', fmt_arg)
            len_arg = arg['attrs'].get('len', False)
            if len_arg:
                fmt_arg.len_arg = len_arg
                fmt_arg.c_var_len = len_arg
                append_format(proto_list, 'int {c_var_len}', fmt_arg)
#            else:
#                fmt_arg.c_var_len = 'NNNN' + fmt_arg.c_var

            # Add any code needed for intent(IN).
            # Usually to convert types.
            # For example, convert char * to std::string
            # Skip input arguments generated by F_string_result_as_arg
            have_cpp_local_var = arg_typedef.cpp_local_var
            for intent in slist:
                have_cpp_local_var = have_cpp_local_var or \
                    c_statements.get(intent, {}).get('cpp_local_var', False)
            if have_cpp_local_var:
                fmt_arg.cpp_var = 'SH_' + fmt_arg.c_var

            for intent in slist:
                # pre_call.append('// intent=%s' % intent)
                if len_trim:
                    cmd_list = c_statements.get(intent, {}) \
                                           .get('pre_call_trim', [])
                else:
                    cmd_list = c_statements.get(intent, {}) \
                                           .get('pre_call', [])
                if cmd_list:
                    for cmd in cmd_list:
                        append_format(pre_call, cmd, fmt_arg)

                cmd_list = c_statements.get(intent, {}).get('post_call', [])
                if cmd_list:
                    # pick up c_str() from cpp_to_c
                    if intent == 'result':
                        fmt_arg.cpp_var = fmt_arg.C_result
                    fmt_arg.cpp_val = wformat(arg_typedef.cpp_to_c, fmt_arg)
                    # append_format(post_call, '// c_var={c_var}  cpp_var={cpp_var}  cpp_val={cpp_val}', fmt_arg)
                    for cmd in cmd_list:
                        append_format(post_call, cmd, fmt_arg)

                cpp_header = c_statements.get(intent, {}) \
                                         .get('cpp_header', None)
                if cpp_header:
                    # include any dependent header in generated source
                    self.header_impl_include[cpp_header] = True

            if arg_typedef.cpp_local_var:
                # cpp_local_var should only be set if c_statements are not used
                if arg_typedef.c_statements:
                    raise RuntimeError(
                        'c_statements and cpp_local_var are both '
                        'defined for {}'
                        .format(arg_typedef.name))
                append_format(pre_call,
                              '{c_const}{cpp_type}{c_ptr} {cpp_var} = ' +
                              arg_typedef.c_to_cpp + ';', fmt_arg)

            if arg_call:
                if have_cpp_local_var:
                    call_list.append(fmt_arg.cpp_var)
                else:
                    # convert C argument to C++
                    append_format(call_list, arg_typedef.c_to_cpp, fmt_arg)

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
        fmt_func.C_call_list = ', '.join(call_list)

        fmt_func.C_prototype = options.get('C_prototype', ', '.join(proto_list))

        fmt_func.cpp_var = fmt_func.C_result
        if node.get('return_this', False):
            fmt_func.C_return_type = 'void'
        else:
            fmt_func.C_return_type = options.get(
                'C_return_type', self._c_type('c_type', result))

        # body of function
        splicer_code = self.splicer_stack[-1].get(fmt_func.function_name, None)
        if 'C_code' in node:
            C_code = [wformat(node['C_code'], fmt_func)]
        elif splicer_code:
            C_code = splicer_code
        else:
            # generate the C body
            C_code = [1]
            C_code.extend(pre_call)
            return_line = ''
            if is_ctor:
                line = wformat('{C_rv_decl} = new {cpp_class}'
                               '({C_call_list});', fmt_func)
                C_code.append(line)
                C_code.append('return '
                              + wformat(result_typedef.cpp_to_c, fmt_func)
                              + ';')
            elif is_dtor:
                C_code.append('delete %s;' % fmt_func.CPP_this)
            elif CPP_subprogram == 'subroutine':
                line = wformat(
                    '{CPP_this_call}{function_name}'
                    '{CPP_template}({C_call_list});',
                    fmt_func)
                C_code.append(line)
                return_line = 'return;'
            else:
                line = wformat(
                    '{C_rv_decl} = {CPP_this_call}{function_name}'
                    '{CPP_template}({C_call_list});',
                    fmt_func)
                C_code.append(line)

                if result_arg:
                    # may be used in C_error_pattern
                    fmt_func.f_string = result_arg['name']
                    fmt_func.f_string_len = result_arg['attrs'].get('len', '')
                    # pick up rv.c_str() from cpp_to_c
                    fmt_func.c_string = wformat(arg_typedef.cpp_to_c, fmt_func)
                if 'C_error_pattern' in node:
                    C_error_pattern = node['C_error_pattern'] + \
                        node.get('_error_pattern_suffix', '')
                    if C_error_pattern in self.patterns:
                        lfmt = util.Options(fmt_func)
                        C_code.append('// check for error')
                        append_format(
                            C_code, self.patterns[C_error_pattern], lfmt)

#                if result_arg:
#                    c_post_call = self.typedef[ result_arg['type'] ].c_post_call
#                    if c_post_call:
#                        append_format(C_code, c_post_call, fmt_func)
                    # XXX release rv is necessary

                if subprogram == 'subroutine':
                    # function result is returned as an argument
                    return_line = 'return;'
                else:
                    return_line = ('return '
                                   + wformat(result_typedef.cpp_to_c, fmt_func)
                                   + ';')

            # copy-out values, clean up
            C_code.extend(post_call)
            if return_line:
                C_code.append(return_line)
            C_code.append(-1)

        self.header_proto_c.append('')
        self.header_proto_c.append(
            wformat('{C_return_type} {C_name}({C_prototype});',
                    fmt_func))

        impl = self.impl
        impl.append('')
        if options.debug:
            impl.append('// %s' % node['_decl'])
            impl.append('// function_index=%d' % node['_function_index'])
        if options.doxygen and 'doxygen' in node:
            self.write_doxygen(impl, node['doxygen'])
        impl.append(wformat('{C_return_type} {C_name}({C_prototype})', fmt_func))
        impl.append('{')
        self._create_splicer(fmt_func.underscore_name +
                             fmt_func.function_suffix, impl, C_code)
        impl.append('}')

    def write_helper_files(self):
        output = [whelpers.FccHeaders]
        self.write_output_file('shroudrt.hpp',
                               self.config.c_fortran_dir, output)

        output = [whelpers.FccCSource]
        self.write_output_file('shroudrt.cpp',
                               self.config.c_fortran_dir, output)
