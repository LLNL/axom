"""
Generate Lua module for C++ code.
"""
from __future__ import print_function

import util
from util import wformat, append_format

def add_templates(options):
    options.update(dict(
        LUA_package_filename_template = 'lua{library}package.cpp',
        LUA_header_filename_template = 'lua{library}package.hpp',
#        LUA_helper_filename_template = 'lua{library}helper.cpp',
#        LUA_PyTypeObject_template    = '{LUA_prefix}{cpp_class}_Type',
#        LUA_PyObject_template        = '{LUA_prefix}{cpp_class}',
#        LUA_type_filename_template   = 'lua{cpp_class}type.cpp',
        LUA_name_impl_template       = '{LUA_prefix}{class_name}{underscore_name}{function_suffix}',
        ))


class Wrapl(util.WrapperMixin):
    """Generate Lua bindings.
    """

    def __init__(self, tree, config, splicers):
        self.tree = tree    # json tree
        self.patterns = tree['patterns']
        self.config = config
        self.log = config.log
        self.typedef = tree['types']
        self._init_splicer(splicers)
        self.comment = '//'

    def reset_file(self):
        pass

    def wrap_library(self):
        top = self.tree
        options = self.tree['options']
        fmt_library = self.tree['fmt']

        # Format variables
        fmt_library.LUA_prefix        = options.get('LUA_prefix', 'l_')
        fmt_library.LUA_package_name  = fmt_library.library_lower
        fmt_library.LUA_state = 'L'
        fmt_library.LUA_package_reg = 'XXX1'
        util.eval_template(top, 'LUA_package_filename')
        util.eval_template(top, 'LUA_header_filename')

        # Variables to accumulate output lines
        self.luaL_Reg_package = []
        self.body_lines = []

        self._push_splicer('class')
        for node in self.tree['classes']:
            name = node['name']
            self.reset_file()
            self._push_splicer(name)
            self.wrap_class(node)
#            self.write_extension_type(node)
            self._pop_splicer(name)
        self._pop_splicer('class')

        self.reset_file()
        if self.tree['functions']:
            self._push_splicer('function')
            self.wrap_functions(None, self.tree['functions'])
            self._pop_splicer('function')

        self.write_header(self.tree)
        self.write_package(self.tree)
#        self.write_helper()

    def wrap_class(self, node):
        options = node['options']
        fmt_class = node['fmt']

        # wrap methods
        self._push_splicer('method')
        self.wrap_functions(node, node['methods'])
        self._pop_splicer('method')

    def wrap_functions(self, cls, functions):
        """Wrap functions for a library or class.
        """
        for function in functions:
            self.wrap_function(cls, function)

    def wrap_function(self, cls, node):
        """Write a Lua wrapper for a C++ function.

        cls  - class node or None for functions
        node - function/method node

        fmt.c_var   - name of variable in PyArg_ParseTupleAndKeywords
        fmt.cpp_var - name of variable in c++ call.
        fmt.py_var  - name of PyObject variable
        """
        options = node['options']
        if not options.wrap_lua:
            return

        if cls:
            cls_function = 'method'
        else:
            cls_function = 'function'
        self.log.write("Lua {0} {1[_decl]}\n".format(cls_function, node))

        fmt_func = node['fmt']
        fmt = util.Options(fmt_func)
        fmt.doc_string = 'documentation'
        util.eval_template(node, 'LUA_name_impl')

        CPP_subprogram = node['_subprogram']

        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)
        result_is_ref = result['attrs'].get('reference', False)

        if node.get('return_this', False):
            result_type = 'void'
            result_is_ptr = False
            CPP_subprogram = 'subroutine'

        result_typedef = self.typedef[result_type]
        is_ctor  = node['attrs'].get('constructor', False)
        is_dtor  = node['attrs'].get('destructor', False)
#        is_const = result['attrs'].get('const', False)
##-        if is_ctor:   # or is_dtor:
##-            # XXX - have explicit delete
##-            # need code in __init__ and __del__
##-            return

        # XXX if a class, then knock off const since the PyObject
        # is not const, otherwise, use const from result.
        if result_typedef.base == 'wrapped':
            is_const = False
        else:
            is_const = None
        fmt.rv_decl = self.std_c_decl('cpp_type', result, name=fmt.rv, const=is_const)  # return value

        LUA_decl = []  # declare variables and pop values
        LUA_code = []  # call C++ function
        LUA_push = []  # push results

        post_parse = []

        cpp_call_list = []

        # parse arguments
        # call function based on number of default arguments provided
        default_calls = []   # each possible default call
        found_default = False
        if '_has_default_arg' in node:
            append_format(LUA_decl, 'int shroud_nargs = lua_gettop({LUA_state});', fmt)

        if True:
            fmt.LUA_index = 1
            for arg in node['args']:
                arg_name = arg['name']
                fmt.c_var = arg['name']
                fmt.cpp_var = fmt.c_var
                fmt.lua_var = 'SH_Lua_' + fmt.c_var
                fmt.c_var_len = 'L' + fmt.c_var
                attrs = arg['attrs']

                lua_pop = None

                arg_typedef = self.typedef[arg['type']]
                LUA_statements = arg_typedef.LUA_statements
                if attrs['intent'] in [ 'inout', 'in']:
                    lua_pop = wformat(arg_typedef.LUA_pop, fmt)
#x                    # names to PyArg_ParseTupleAndKeywords
#X                    arg_names.append(arg_name)
#X                    arg_offsets.append( '(char *) kwcpp+%d' % offset)
#X                    offset += len(arg_name) + 1

                    # XXX default should be handled differently
                    if 'default' in attrs:
                        if not found_default:
#                            parse_format.append('|')  # add once
                            found_default = True
                        # call for default arguments  (num args, arg string)
                        default_calls.append(
                            (len(cpp_call_list), len(post_parse), ', '.join(cpp_call_list)))
                        LUA_code.extend([
                                'if (shroud_nargs > {}) {{'.format(fmt.LUA_index-1),
                                1,
                                '{} = {}'.format(fmt.c_var, lua_pop),
                                -1,
                                '}'
                                ])
                        lua_pop = False;

#                    parse_format.append(arg_typedef.PY_format)
#                    if arg_typedef.PY_PyTypeObject:
#                        # Expect object of given type
#                        parse_format.append('!')
#                        parse_vargs.append('&' + arg_typedef.PY_PyTypeObject)
#                        arg_name = fmt.py_var
#                    elif arg_typedef.PY_from_object:
#                        # Use function to convert object
#                        parse_format.append('&')
#                        parse_vargs.append(arg_typedef.PY_from_object)

                    # add argument to call to PyArg_ParseTypleAndKeywords
#                    parse_vargs.append('&' + arg_name)

#                    cmd_list = py_statements.get('intent_in',{}).get('post_parse',[])
#                    if cmd_list:
#                        fmt.cpp_var = 'SH_' + fmt.c_var
#                        for cmd in cmd_list:
#                            append_format(post_parse, cmd, fmt)
                    fmt.LUA_index += 1 

                if attrs['intent'] in [ 'inout', 'out']:
                    # output variable must be a pointer
                    # XXX - fix up for strings
#                    format, vargs = self.intent_out(arg_typedef, fmt, post_call)
#                    build_format.append(format)
#                    build_vargs.append('*' + vargs)
                    append_format(LUA_push, arg_typedef.LUA_push, fmt)
                   

                # argument for C++ function
                lang = 'cpp_type'
                if arg_typedef.base == 'string':
                    # C++ will coerce char * to std::string
                    lang = 'c_type'
                if attrs.get('reference', False):
                    # convert a reference to a pointer
                    ptr = True
                else:
                    ptr = False

                if lua_pop:
                    decl_suffix = ' = {};'.format(lua_pop)
                else:
                    decl_suffix = ';'
                LUA_decl.append(self.std_c_decl(lang, arg, const=False, ptr=ptr) + decl_suffix)
                
#                if arg_typedef.PY_PyTypeObject:
#                    # A Python Object which must be converted to C++ type.
#                    objtype = arg_typedef.PY_PyObject or 'PyObject'
#                    LUA_decl.append(objtype + ' * ' + fmt.py_var + ';')
#                    cpp_call_list.append(fmt.cpp_var)
#                elif arg_typedef.PY_from_object:
#                    # already a C++ type
#                    cpp_call_list.append(fmt.cpp_var)
#                else:
                # convert to C++ type
                fmt.ptr=' *' if arg['attrs'].get('ptr', False) else ''
                append_format(cpp_call_list, arg_typedef.c_to_cpp, fmt)

        if cls:
#                    template = '{C_const}{cpp_class} *{C_this}obj = static_cast<{C_const}{cpp_class} *>(static_cast<{C_const}void *>({C_this}));'
#                fmt_func.C_object = wformat(template, fmt_func)
            fmt.LUA_this_call = wformat('self->{BBB}->', fmt)  # call method syntax
##            if ctor: add to package else add to class meta-table
                
        else:
            fmt.LUA_this_call = ''  # call function syntax
            self.luaL_Reg_package.append(wformat('{{"{function_name}{function_suffix}", {LUA_name_impl}}},', fmt))

        # call with all arguments
        default_calls.append(
            (len(cpp_call_list),  len(post_parse), ', '.join(cpp_call_list)))

        # If multiple calls, declare return value once
        # Else delare on call line.
        if found_default:
            fmt.rv_asgn = 'rv = '
            LUA_code.append('switch (shroud_nargs) {')
        else:
            fmt.rv_asgn = fmt.rv_decl + ' = '
        need_rv = False

        for nargs, len_post_parse, call_list in default_calls:
            if found_default:
                LUA_code.append('case %d:' % nargs)
                LUA_code.append(1)

            fmt.call_list = call_list
            LUA_code.extend(post_parse[:len_post_parse])

            if is_dtor:
                append_format(LUA_code, 'delete self->{BBB};', fmt)
                append_format(LUA_code, 'self->{BBB} = NULL;', fmt)
            elif CPP_subprogram == 'subroutine':
                line = wformat('{LUA_this_call}{function_name}({call_list});', fmt)
                LUA_code.append(line)
            else:
                need_rv = True
                line = wformat('{rv_asgn}{LUA_this_call}{function_name}({call_list});', fmt)
                LUA_code.append(line)

#            if 'PY_error_pattern' in node:
#                lfmt = util.Options(fmt)
#                lfmt.c_var = fmt.rv
#                lfmt.cpp_var = fmt.rv
#                append_format(LUA_code, self.patterns[node['PY_error_pattern']], lfmt)

            if found_default:
                LUA_code.append('break;')
                LUA_code.append(-1)
        if found_default:
#            LUA_code.append('default:')
#            LUA_code.append(1)
#            LUA_code.append('continue;')  # XXX raise internal error
#            LUA_code.append(-1)
            LUA_code.append('}')
        else:
            need_rv = False

        if need_rv:
            LUA_decl.append(fmt.rv_decl + ';')
        if len(LUA_decl):
            LUA_decl.append('')

        # Compute return value
        if CPP_subprogram == 'function':
            fmt.cpp_var = fmt.rv
            fmt.c_var = wformat(result_typedef.cpp_to_c, fmt)  # if C++
            append_format(LUA_push, result_typedef.LUA_push, fmt)
#            format, vargs = self.intent_out(result_typedef, fmt, PY_code)
#            # Add result to front of result tuple
#            build_format.insert(0, format)
#            build_vargs.insert(0, vargs)

#        LUA_code.extend(post_call)


        body = self.body_lines
        body.extend([
                '',
                wformat('static int {LUA_name_impl}(lua_State *L)', fmt),
                '{',
                1])
        body.extend(LUA_decl)
        body.extend(LUA_code)
        body.extend(LUA_push)    # return values
        body.extend([
                'return {};'.format(len(LUA_push)),
                -1,
                '}'
                ])

    def write_header(self, node):
        options = node['options']
        fmt = node['fmt']
        fname = fmt.LUA_header_filename

        output = []

        # add guard
        guard = fname.replace(".", "_").upper()
        output.extend([
                '#ifndef %s' % guard,
                '#define %s' % guard,
                '#include "lua.h"',
                wformat('int luaopen_{LUA_package_name} (lua_state *{LUA_state});', fmt),
                ])
        output.append('#endif  /* %s */' % guard)
        self.write_output_file(fname, self.config.python_dir, output)

    def append_luaL_Reg(self, output, name, lines):
        """Create luaL_Reg struct"""
        output.extend([
                'static const struct luaL_Reg {} [] = {{'.format(name),
                1,
                ])
        output.extend(self.lines)
        output.extend([
                '{NULL, NULL}   /*sentinel */',
                -1,
                '};',
                ''
                ])

    def write_package(self, node):
        options = node['options']
        fmt = node['fmt']
        fname = fmt.LUA_package_filename

        output = []

        if options.cpp_header:
            for include in options.cpp_header.split():
                output.append('#include "%s"' % include)
        output.append(wformat('#include "{LUA_header_filename}"', fmt))
        self._create_splicer('include', output)

        output.append('')
        self.namespace(node, 'begin', output)
        self._create_splicer('C_definition', output)

        output.extend(self.body_lines)

        self.append_luaL_Reg(output, fmt.LUA_package_reg, self.lua_Reg_package)
        output.extend([
                wformat('int luaopen_{LUA_package_name} (lua_state *{LUA_state}) {{', fmt),
                1,
                wformat('luaL_newLib({LUA_state}, {LUA_package_reg});', fmt),
                'return 1',
                -1,
                '}'
                ])

        self.namespace(node, 'end', output)

        self.write_output_file(fname, self.config.lua_dir, output)
