#!/bin/env python3
"""
Generate C bindings for C++ classes

"""
from __future__ import print_function

import whelpers
import util
from util import append_format

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
        fmt_library = self.tree['fmt']
        fmt_library.C_const = ''

        self._push_splicer('class')
        for node in self.tree['classes']:
            self._push_splicer(node['name'])
            self.write_file(node, self.wrap_class, True)
            self._pop_splicer(node['name'])
        self._pop_splicer('class')

        if self.tree['functions']:
            self.write_file(self.tree, self.wrap_functions, False)

        self.write_helper_files()

    def write_file(self, node, worker, cls):
        """Write a file for the library and its functions or
        a class and its methods.
        """
        fmt = node['fmt']
        self._begin_output_file()
        worker(node)    # self.wrap_class or self.wrap_functions
        c_header = fmt.C_header_filename
        c_impl   = fmt.C_impl_filename
        self.write_header(node, c_header, cls)
        self.write_impl(node, c_header, c_impl, cls)

    def wrap_functions(self, tree):
        # worker function for write_file
        self._push_splicer('function')
        for node in tree['functions']:
            self.wrap_function(None, node)
        self._pop_splicer('function')

    def write_header(self, node, fname, cls=False):
        guard = fname.replace(".", "_").upper()
        options = node['options']

        output = []

        if options.doxygen:
            self.write_doxygen_file(output, fname, node, cls)

        output.extend([
                '// For C users and C++ implementation',
                '',
                '#ifndef %s' % guard,
                '#define %s' % guard,
                ])
        # headers required by typedefs
        if self.header_typedef_include:
#            output.append('// header_typedef_include')
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
        names = self.header_forward.keys()
        names.sort()
        for name in names:
            output.append('struct s_{C_type_name};\ntypedef struct s_{C_type_name} {C_type_name};'.
                     format(C_type_name=name))
        output.append('')
        self._create_splicer('C_definition', output)
        output.extend(self.header_proto_c);
        output.extend([
                '',
                '#ifdef __cplusplus',
                '}',
                '#endif',
                '',
                '#endif  // %s' % guard
                ])

        self.config.cfiles.append(fname)
        self.write_output_file(fname, self.config.binary_dir, output)

    def write_impl(self, node, hname, fname, cls=False):
        # node = class node
        options = node['options']
        namespace = options.namespace

        output = []
        output.append('// ' + fname)

        output.append('#include "%s"' % hname)
        if options.cpp_header:
            for include in options.cpp_header.split():
                self.header_impl_include[include] = True
        # headers required by implementation
        if self.header_impl_include:
            headers = self.header_impl_include.keys()
            self.write_headers(headers, output)

        output.append('\nextern "C" {')
        self.namespace(node, 'begin', output)
        output.extend(self.impl)
        output.append('')
        self._create_splicer('additional_functions', output)
        output.append('')
        self.namespace(node, 'end', output)

        output.append('}  // extern "C"')

        self.config.cfiles.append(fname)
        self.write_output_file(fname, self.config.binary_dir, output)

    def write_headers(self, headers, output):
        headers.sort()
        for header in headers:
            if header[0] == '<':
                output.append('#include %s' % header)
            else:
                output.append('#include "%s"' % header)
        
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
        fmt = util.Options(fmt_func)
        if 'CPP_template' not in fmt:
            fmt.CPP_template = ''

        # Look for C++ routine to wrap
        # Usually the same node unless it is bufferified
        CPP_node = node
        generated = []
        if '_generated' in CPP_node:
            generated.append(CPP_node['_generated'])
        while '_PTR_C_CPP_index' in CPP_node:
            CPP_node = self.tree['function_index'][CPP_node['_PTR_C_CPP_index']]
            if '_generated' in CPP_node:
                generated.append(CPP_node['_generated'])
        CPP_result = CPP_node['result']
        CPP_result_type = CPP_result['type']
        CPP_subprogram = CPP_node['_subprogram']

        # C return type
        result = node['result']
        result_type = result['type']
        subprogram = node['_subprogram']

        # C++ functions which return 'this', are easier to call from Fortran if they are subroutines.
        # There is no way to chain in Fortran:  obj->doA()->doB();
        if node.get('return_this', False):
            CPP_result_type = 'void'
            CPP_subprogram = 'subroutine'

        result_typedef = self.typedef[result_type]
        is_const = result['attrs'].get('const', False)
        is_ctor  = node['attrs'].get('constructor', False)
        is_dtor  = node['attrs'].get('destructor', False)

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
            fmt.C_const = 'const '
        fmt.CPP_this = fmt_func.C_this + 'obj'
        fmt.rv_decl = self._c_decl('cpp_type', CPP_result, name='rv')  # return value

        proto_list = []
        call_list = []
        if cls:
            # object pointer
            fmt.CPP_this_call = fmt.CPP_this + '->'  # call method syntax
            arg_dict = dict(name=fmt_func.C_this,
                            type=cls['name'], 
                            attrs=dict(ptr=True,
                                       const=is_const))
            C_this_type = self._c_type('c_type', arg_dict)
            if not is_ctor:
                arg = self._c_decl('c_type', arg_dict)
                proto_list.append(arg)
        else:
            fmt.CPP_this_call = ''  # call function syntax

        result_arg = None  # indicate which argument contains function result, usually none
        pre_call   = []    # list of temporary variable declarations
        post_call  = []

#
#    c_var     - argument to C function  (wrapper function)
#    c_var_len - variable with trimmed length of c_var
#    c_var_num - variable with length of c_var
#    cpp_var   - argument to C++ function  (wrapped function).
#                Usually same as c_var but may be a new local variable
#                or the funtion result variable.
#
#
#
        for arg in node['args']:
            c_attrs = arg['attrs']
            arg_typedef = self.typedef[arg['type']]
            c_statements = arg_typedef.c_statements
            fmt.c_var = arg['name']
            fmt.c_var_len = 'L' + fmt.c_var
            fmt.c_var_num = 'N' + fmt.c_var
            fmt.ptr = ' *' if arg['attrs'].get('ptr', False) else ''

            if c_attrs.get('_is_result', False):
                arg_call = False
                result_arg = arg
#                fmt.var = 'rv'           # name of result variable
#                fmt.c_var = 'rv'           # name of result variable
#                fmt.cpp_var = 'rv'           # name of result variable
                slist = [ 'result' ]
            else:
                arg_call = arg
                fmt.cpp_var = fmt.c_var      # name in c++ call.
                slist = [ ]
                if c_attrs['intent'] in [ 'inout', 'in']:
                    slist.append('intent_in')
                if c_attrs['intent'] in [ 'inout', 'out']:
                    slist.append('intent_out')

            proto_list.append(self._c_decl('c_type', arg))

            len_trim = arg['attrs'].get('len_trim', False)
            if len_trim:
                fmt.c_var_trim = len_trim
                append_format(proto_list, 'int {c_var_trim}', fmt)
            len_arg = arg['attrs'].get('len', False)
            if len_arg:
                fmt.len_arg = len_arg
                fmt.c_var_len = len_arg
                append_format(proto_list, 'int {c_var_len}', fmt)
#            else:
#                fmt.c_var_len = 'NNNN' + fmt.c_var

            # Add any code needed for intent(IN).
            # Usually to convert types. For example, convert char * to std::string
            # Skip input arguments generated by F_string_result_as_arg
            for intent in slist:
#                pre_call.append('// intent=%s' % intent)
                if len_trim:
                    cmd_list = c_statements.get(intent,{}).get('pre_call_trim',[])
                else:
                    cmd_list = c_statements.get(intent,{}).get('pre_call',[])
                if cmd_list:
                    fmt.cpp_var = 'SH_' + fmt.c_var
                    for cmd in cmd_list:
                        append_format(pre_call, cmd, fmt)

                cmd_list = c_statements.get(intent,{}).get('post_call',[])
                if cmd_list:
                    # pick up c_str() from cpp_to_c
                    if intent == 'result':
                        fmt.cpp_var = 'rv'
                    fmt.cpp_val = wformat(arg_typedef.cpp_to_c, fmt)
#                    append_format(post_call, '// c_var={c_var}  cpp_var={cpp_var}  cpp_val={cpp_val}', fmt)
                    for cmd in cmd_list:
                        append_format(post_call, cmd, fmt)

                cpp_header = c_statements.get(intent,{}).get('cpp_header',None)
                if cpp_header:
                    # include any dependent header in generated source
                    self.header_impl_include[cpp_header] = True

            if arg_call:
                # convert C argument to C++
                append_format(call_list, arg_typedef.c_to_cpp, fmt)

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
        fmt.C_call_list = ', '.join(call_list)

        fmt.C_prototype = options.get('C_prototype', ', '.join(proto_list))

        fmt.cpp_var = 'rv'
        if node.get('return_this', False):
            fmt.C_return_type = 'void'
        else:
            fmt.C_return_type = options.get('C_return_type', self._c_type('c_type', result))

        fmt.C_object = ''
        if cls:
            if 'C_object' in options:
                fmt.C_object = options.C_object
            elif not is_ctor:
                template = '{C_const}{cpp_class} *{C_this}obj = static_cast<{C_const}{cpp_class} *>(static_cast<{C_const}void *>({C_this}));'
                fmt.C_object = wformat(template, fmt)

        # body of function
        splicer_code = self.splicer_stack[-1].get(fmt_func.function_name, None)
        if 'C_code' in options:
            C_code = [   wformat(options.C_code, fmt) ]
        elif splicer_code:
            C_code = splicer_code
        else:
            # generate the C body
            C_code = []
            C_code.extend(pre_call)
            return_line = ''
            if is_ctor:
                line = wformat('{rv_decl} = new {cpp_class}({C_call_list});', fmt)
                C_code.append(line)
                C_code.append('return ' + wformat(result_typedef.cpp_to_c, fmt) + ';')
            elif is_dtor:
                C_code.append('delete %sobj;' % fmt_func.C_this)
            elif CPP_subprogram == 'subroutine':
                line = wformat('{CPP_this_call}{function_name}{CPP_template}({C_call_list});',
                               fmt)
                C_code.append(line)
                return_line = 'return;'
            else:
                line = wformat('{rv_decl} = {CPP_this_call}{function_name}{CPP_template}({C_call_list});',
                               fmt)
                C_code.append(line)

                if result_arg:
                    # may be used in C_error_pattern
                    fmt.f_string = result_arg['name']
                    fmt.f_string_len = result_arg['attrs'].get('len', '')
                    fmt.c_string = wformat(arg_typedef.cpp_to_c, fmt)  # pick up rv.c_str() from cpp_to_c
                if 'C_error_pattern' in node:
                    C_error_pattern = node['C_error_pattern'] + \
                        node.get('_error_pattern_suffix', '')
                    if C_error_pattern in self.patterns:
                        lfmt = util.Options(fmt)
#                        lfmt.var = fmt.rv
                        C_code.append('// check for error')
                        append_format(C_code, self.patterns[C_error_pattern], lfmt)

#                if result_arg:
#                    c_post_call = self.typedef[ result_arg['type'] ].c_post_call
#                    if c_post_call:
#                        append_format(C_code, c_post_call, fmt)
                    # XXX release rv is necessary

                if subprogram == 'subroutine':
                    # function result is returned as an argument
                    return_line = 'return;'
                else:
                    return_line = 'return ' + wformat(result_typedef.cpp_to_c, fmt) + ';'

            # copy-out values, clean up
            C_code.extend(post_call)

            if return_line:
                C_code.append(return_line)


        self.header_proto_c.append('')
        self.header_proto_c.append(wformat('{C_return_type} {C_name}({C_prototype});',
                                           fmt))

        impl = self.impl
        impl.append('')
        if options.debug:
            impl.append('// %s' % node['_decl'])
            impl.append('// function_index=%d' % node['_function_index'])
        if options.doxygen and 'doxygen' in node:
            self.write_doxygen(impl, node['doxygen'])
        impl.append(wformat('{C_return_type} {C_name}({C_prototype})', fmt))
        impl.append('{')
        if cls:
            impl.append(fmt.C_object )
        self._create_splicer(fmt_func.underscore_name + 
                             fmt_func.function_suffix, impl, C_code)
        impl.append('}')

    def write_helper_files(self):
        output = [ whelpers.FccHeaders ]
        self.write_output_file('shroudrt.hpp', self.config.binary_dir, output)
