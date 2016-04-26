"""
Generate Lua module for C++ code.
"""
from __future__ import print_function

import util
from util import wformat, append_format

def add_templates(options):
    options.update(dict(
        LUA_module_filename_template = 'lua{library}module.cpp',
        LUA_header_filename_template = 'lua{library}module.hpp',
        LUA_userdata_type_template = '{LUA_prefix}{cpp_class}_Type',
        LUA_userdata_member_template = 'self',
        LUA_class_reg_template = '{LUA_prefix}{cpp_class}_Reg',
        LUA_metadata_template  = '{cpp_class}.metatable',
        LUA_ctor_name_template = '{cpp_class}',
        LUA_name_template = '{function_name}{function_suffix}',
        LUA_name_impl_template = '{LUA_prefix}{class_name}{underscore_name}{function_suffix}',
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
        fmt_library.LUA_module_name  = fmt_library.library_lower
        fmt_library.LUA_state_var = 'L'
        fmt_library.LUA_module_reg = 'XXX1'
        util.eval_template(top, 'LUA_module_filename')
        util.eval_template(top, 'LUA_header_filename')

        # Variables to accumulate output lines
        self.luaL_Reg_module = []
        self.body_lines = []
        self.class_lines = []
        self.lua_type_structs = []

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
        self.write_module(self.tree)
#        self.write_helper()

    def wrap_class(self, node):
        options = node['options']
        fmt_class = node['fmt']

        fmt_class.LUA_userdata_var = 'SH_this'
        util.eval_template(node, 'LUA_userdata_type')
        util.eval_template(node, 'LUA_userdata_member')
        util.eval_template(node, 'LUA_class_reg')
        util.eval_template(node, 'LUA_metadata')
        util.eval_template(node, 'LUA_ctor_name')

#        self._create_splicer('C_declaration', self.lua_type_structs)
        self.lua_type_structs.append('')
        self.lua_type_structs.append('typedef struct {')
        self.lua_type_structs.append(1)
        append_format(self.lua_type_structs, '{namespace_scope}{cpp_class} * {LUA_userdata_member};', fmt_class)
#        self._create_splicer('C_object', self.lua_type_structs)
        self.lua_type_structs.append(-1)
        self.lua_type_structs.append(wformat('}} {LUA_userdata_type};', fmt_class))

        self.luaL_Reg_class = []

        # wrap methods
        self._push_splicer('method')
        self.wrap_functions(node, node['methods'])
        self.append_luaL_Reg(self.body_lines, fmt_class.LUA_class_reg,
                             self.luaL_Reg_class)

        append_format(self.class_lines, """
/* Create the metatable and put it on the stack. */
luaL_newmetatable({LUA_state_var}, "{LUA_metadata}");
/* Duplicate the metatable on the stack (We now have 2). */
lua_pushvalue({LUA_state_var}, -1);
/* Pop the first metatable off the stack and assign it to __index
 * of the second one. We set the metatable for the table to itself.
 * This is equivalent to the following in lua:
 * metatable = {{}}
 * metatable.__index = metatable
 */
lua_setfield({LUA_state_var}, -2, "__index");
 
/* Set the methods to the metatable that should be accessed via object:func */
luaL_setfuncs({LUA_state_var}, {LUA_class_reg}, 0);""", fmt_class)

        self._pop_splicer('method')

    def wrap_functions(self, cls, functions):
        """Wrap functions for a library or class.
        Create one wrapper for overloaded functions and the
        different variations of default-arguments
        """

        # Find overloaded functions.
        # maintain order, but gather overloads together
        overloaded_methods = {}
        overloads = []
        for function in functions:
            if not function['options'].wrap_lua:
                continue
            name = function['result']['name']
            if name in overloaded_methods:
                overloaded_methods[name].append(function)
            else:
                first = [ function ]
                overloads.append(first)
                overloaded_methods[name] = first

        for overload in overloads:
            self.wrap_function2(cls, overload)

#        for function in functions:
#            self.wrap_function(cls, function)

    def XXXwrap_function(self, cls, node):
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
        util.eval_template(node, 'LUA_name')
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
        if is_dtor:
            fmt.LUA_name = '__gc'

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

        # find class object
        if cls:
            cls_typedef = self.typedef[cls['name']]
            if not is_ctor:
                fmt.c_var = wformat(cls_typedef.LUA_pop, fmt)
                LUA_code.append(
                    wformat('{LUA_userdata_type} * {LUA_userdata_var} = {c_var};', fmt))

        # parse arguments
        # call function based on number of default arguments provided
        default_calls = []   # each possible default call
        found_default = False
        if '_has_default_arg' in node:
            append_format(LUA_decl, 'int SH_nargs = lua_gettop({LUA_state_var});', fmt)

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
#                    lua_pop = wformat(arg_typedef.LUA_pop, fmt)
                    # lua_pop is a C++ expression
                    fmt.c_var = wformat(arg_typedef.LUA_pop, fmt)
                    lua_pop = wformat(arg_typedef.c_to_cpp, fmt)

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
                                'if (SH_nargs > {}) {{'.format(fmt.LUA_index-1),
                                1,
                                '{} = {};'.format(fmt.cpp_var, lua_pop),
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

                    #append_format(LUA_push, arg_typedef.LUA_push, fmt)
                    tmp = wformat(arg_typedef.LUA_push, fmt)
                    LUA_push.append( tmp + ';' )
                   

                # argument for C++ function
                lang = 'cpp_type'
                arg_const = False
                if arg_typedef.base == 'string':
                    # C++ will coerce char * to std::string
                    lang = 'c_type'
                    arg_const = True  # lua_tostring is const
                if attrs.get('reference', False):
                    # convert a reference to a pointer
                    ptr = True
                else:
                    ptr = False

                if lua_pop:
                    decl_suffix = ' = {};'.format(lua_pop)
                else:
                    decl_suffix = ';'
                LUA_decl.append(self.std_c_decl(lang, arg, const=arg_const, ptr=ptr) + decl_suffix)
                
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
#                fmt.ptr=' *' if arg['attrs'].get('ptr', False) else ''
#                append_format(cpp_call_list, arg_typedef.c_to_cpp, fmt)
                cpp_call_list.append(fmt.cpp_var)

        if cls:
#                    template = '{C_const}{cpp_class} *{C_this}obj = static_cast<{C_const}{cpp_class} *>(static_cast<{C_const}void *>({C_this}));'
#                fmt_func.C_object = wformat(template, fmt_func)
            fmt.LUA_this_call = wformat('{LUA_userdata_var}->{LUA_userdata_member}->', fmt)  # call method syntax
            if is_ctor:
                self.luaL_Reg_module.append(wformat('{{"{LUA_ctor_name}{function_suffix}", {LUA_name_impl}}},', fmt))
            else:
                self.luaL_Reg_class.append(wformat('{{"{LUA_name}", {LUA_name_impl}}},', fmt))
                
        else:
            fmt.LUA_this_call = ''  # call function syntax
            self.luaL_Reg_module.append(wformat('{{"{LUA_name}", {LUA_name_impl}}},', fmt))

        # call with all arguments
        default_calls.append(
            (len(cpp_call_list),  len(post_parse), ', '.join(cpp_call_list)))

        # If multiple calls, declare return value once
        # Else delare on call line.
        if found_default:
            fmt.rv_asgn = 'rv = '
            LUA_code.append('switch (SH_nargs) {')
        else:
            fmt.rv_asgn = fmt.rv_decl + ' = '
        need_rv = False

        for nargs, len_post_parse, cpp_call_list in default_calls:
            if found_default:
                LUA_code.append('case %d:' % nargs)
                LUA_code.append(1)

            fmt.cpp_call_list = cpp_call_list
            LUA_code.extend(post_parse[:len_post_parse])

            if is_ctor:
                LUA_code.extend([
                        wformat('{LUA_userdata_type} * {LUA_userdata_var} = ({LUA_userdata_type} *) lua_newuserdata({LUA_state_var}, sizeof(*{LUA_userdata_var}));', fmt),
                        wformat('{LUA_userdata_var}->{LUA_userdata_member} = new {cpp_class}({cpp_call_list});', fmt),
                        '/* Add the metatable to the stack. */',
                        wformat('luaL_getmetatable(L, "{LUA_metadata}");', fmt),
                        '/* Set the metatable on the userdata. */',
                        'lua_setmetatable(L, -2);',
                        ])
            elif is_dtor:
                LUA_code.extend([
                        wformat('delete {LUA_userdata_var}->{LUA_userdata_member};', fmt),
                        wformat('{LUA_userdata_var}->{LUA_userdata_member} = NULL;', fmt),
                        ])
            elif CPP_subprogram == 'subroutine':
                line = wformat('{LUA_this_call}{function_name}({cpp_call_list});', fmt)
                LUA_code.append(line)
            else:
                need_rv = True
                line = wformat('{rv_asgn}{LUA_this_call}{function_name}({cpp_call_list});', fmt)
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
        if CPP_subprogram == 'function' and not is_ctor:
            fmt.cpp_var = fmt.rv
            fmt.c_var = wformat(result_typedef.cpp_to_c, fmt)  # if C++
##            append_format(LUA_push, result_typedef.LUA_push, fmt)
            tmp = wformat(result_typedef.LUA_push, fmt)
            LUA_push.append( tmp + ';' )
#            format, vargs = self.intent_out(result_typedef, fmt, PY_code)
#            # Add result to front of result tuple
#            build_format.insert(0, format)
#            build_vargs.insert(0, vargs)

#        LUA_code.extend(post_call)

        num_return = len(LUA_push)
        if is_ctor:
            # ctor uses lua_newuserdata which leaves on stack (no push).
            num_return += 1

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
                'return {};'.format(num_return),
                -1,
                '}'
                ])

    def wrap_function2(self, cls, overloads):
        """Write a Lua wrapper for a C++ function.

        cls  - class node or None for functions
        overloads - a list of functions to wrap.

        fmt.c_var   - name of variable in PyArg_ParseTupleAndKeywords
        fmt.cpp_var - name of variable in c++ call.
        """

        # First overload defines options
        node = overloads[0]

        function_name = node['result']['name']
        fmt_func = node['fmt']
        fmt = util.Options(fmt_func)
        util.eval_template(node, 'LUA_name')
        util.eval_template(node, 'LUA_name_impl')

        is_ctor  = node['attrs'].get('constructor', False)
        is_dtor  = node['attrs'].get('destructor', False)
        if is_dtor:
            fmt.LUA_name = '__gc'

        if cls:
            fmt.LUA_this_call = wformat('{LUA_userdata_var}->{LUA_userdata_member}->', fmt)  # call method syntax
        else:
            fmt.LUA_this_call = ''  # call function syntax

        # Loop over all overloads and default argument and
        # sort by the number of arguments expected.

        # Create lists of input arguments
        all_calls = []
        maxargs = 0
        for function in overloads:
            nargs = 0
            in_args = []
            out_args = []
            found_default = False
            for arg in function['args']:
                arg_typedef = self.typedef[arg['type']]
                attrs = arg['attrs']

#                # argument for C++ function
#                lang = 'cpp_type'
#                arg_const = False
#                if arg_typedef.base == 'string':
#                    # C++ will coerce char * to std::string
#                    lang = 'c_type'
#                    arg_const = True  # lua_tostring is const
#                if attrs.get('reference', False):
#                    # convert a reference to a pointer
#                    ptr = True
#                else:
#                    ptr = False

                nargs += 1 
#                decl = self.std_c_decl(lang, arg, const=arg_const, ptr=ptr)
#                in_args.append(decl)
                in_args.append(arg)
                if 'default' in attrs:
                    all_calls.append(lua_function(function, in_args[:], out_args))
                    found_default = True
            if not found_default:
                # no defaults, use all arguments
                all_calls.append(lua_function(function, in_args[:], out_args))
            maxargs = max(maxargs, nargs)


        # Gather calls by number of arguments
        by_count = [ [] for i in range(maxargs+1) ]
        for a_call in all_calls:
            by_count[ a_call.nargs ].append(a_call)

        body = self.body_lines
        body.extend([
                '',
                wformat('static int {LUA_name_impl}(lua_State *L)', fmt),
                '{',
                1,
                'int SH_nresult;'
                ])


        if len(all_calls) == 1:
            self.do_function(cls, all_calls[0], fmt)
        else:
            itype_vars = []
            for iarg in range(1,maxargs+1):
                itype_vars.append('SH_itype{}'.format(iarg))
                body.append('int {};'.format(itype_vars[-1]))
            append_format(body, 'int SH_nargs = lua_gettop({LUA_state_var});', fmt)
            body.append('switch (SH_nargs) {')
            for nargs, calls in enumerate(by_count):
                if len(calls) == 0:
                    continue
                body.append('case {}:'.format(nargs))
                body.append(1)
                ifelse = 'if'

                # Local variables to hold argument types
                for iarg in range(nargs):
                    fmt.itype_var = itype_vars[iarg]
                    fmt.iarg = iarg + 1
                    append_format(body, '{itype_var} = lua_type({LUA_state_var}, {iarg});', fmt)

                for call in calls:
                    checks = []
                    for iarg, arg in enumerate(call.inargs):
                        arg_typedef = self.typedef[arg['type']]
                        fmt.itype_var = itype_vars[iarg]
                        fmt.itype = arg_typedef.LUA_type
                        append_format(checks, '{itype_var} == {itype}', fmt)
                    if nargs == 0:
                        self.do_function(cls, call, fmt)
                    elif nargs == 1:
                        body.append('{} ({}) {{'.format(ifelse, checks[0]))
                        body.append(1)
                        self.do_function(cls, call, fmt)
                        body.append(-1)
                        body.append('}')
                    elif nargs == 2:
                        body.append('{} ({} &&'.format(ifelse, checks[0]))
                        body.append(1)
                        body.append('{}) {{'.format(checks[1]))
                        self.do_function(cls, call, fmt)
                        body.append(-1)
                        body.append('}')
                    else:
                        body.append('{} ({} &&'.format(ifelse, checks[0]))
                        body.append(1)
                        for check in checks[1:-1]:
                            body.append('{} &&'.format(check))
                        body.append('{}) {{'.format(checks[-1]))
                        self.do_function(cls, call, fmt)
                        body.append(-1)
                        body.append('}')
                    ifelse = 'else if'
                if nargs > 0:
                    # Trap errors when the argument types do not match
                    body.append('else {')
                    body.append(1)
                    body.append('// raise some error')
                    body.append(-1)
                    body.append('}')
                body.append('break;')
                body.append(-1)
            body.append('default:')
            body.append(1)
            body.append('break;')
            body.append(-1)
            body.append('}')

        body.extend([
                'return SH_nresult;',
                -1,
                '}'
                ])

        # Save pointer to function
        if cls:
            if is_ctor:
                self.luaL_Reg_module.append(wformat('{{"{LUA_ctor_name}{function_suffix}", {LUA_name_impl}}},', fmt))
            else:
                self.luaL_Reg_class.append(wformat('{{"{LUA_name}", {LUA_name_impl}}},', fmt))
                
        else:
            self.luaL_Reg_module.append(wformat('{{"{LUA_name}", {LUA_name_impl}}},', fmt))


    def do_function(self, cls, luafcn, fmt):
        """
        Wrap a single function/overload/default-argument
        variation of a function.

        fmt - local format dictionary
        """
        node = luafcn.function

        if cls:
            cls_function = 'method'
        else:
            cls_function = 'function'
        self.log.write("Lua {0} {1[_decl]}\n".format(cls_function, node))

#        fmt_func = node['fmt']
#        fmt = util.Options(fmt_func)
#        fmt.doc_string = 'documentation'
#        util.eval_template(node, 'LUA_name')
#        util.eval_template(node, 'LUA_name_impl')

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

        # find class object
        if cls:
            cls_typedef = self.typedef[cls['name']]
            if not is_ctor:
                fmt.c_var = wformat(cls_typedef.LUA_pop, fmt)
                LUA_code.append(
                    wformat('{LUA_userdata_type} * {LUA_userdata_var} = {c_var};', fmt))

        # parse arguments
        # call function based on number of default arguments provided
#        default_calls = []   # each possible default call
#        found_default = False
#        if '_has_default_arg' in node:
#            append_format(LUA_decl, 'int SH_nargs = lua_gettop({LUA_state_var});', fmt)

        if True:
            fmt.LUA_index = 1
            # Only process nargs.
            # Each variation of default-arguments produces a new call.
            for iarg in range(luafcn.nargs):
                arg = node['args'][iarg]
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
#                    lua_pop = wformat(arg_typedef.LUA_pop, fmt)
                    # lua_pop is a C++ expression
                    fmt.c_var = wformat(arg_typedef.LUA_pop, fmt)
                    lua_pop = wformat(arg_typedef.c_to_cpp, fmt)
                    fmt.LUA_index += 1 

                if attrs['intent'] in [ 'inout', 'out']:
                    # output variable must be a pointer
                    # XXX - fix up for strings
#                    format, vargs = self.intent_out(arg_typedef, fmt, post_call)
#                    build_format.append(format)
#                    build_vargs.append('*' + vargs)

                    #append_format(LUA_push, arg_typedef.LUA_push, fmt)
                    tmp = wformat(arg_typedef.LUA_push, fmt)
                    LUA_push.append( tmp + ';' )

                # argument for C++ function
                lang = 'cpp_type'
                arg_const = False
                if arg_typedef.base == 'string':
                    # C++ will coerce char * to std::string
                    lang = 'c_type'
                    arg_const = True  # lua_tostring is const
                if attrs.get('reference', False):
                    # convert a reference to a pointer
                    ptr = True
                else:
                    ptr = False

                if lua_pop:
                    decl_suffix = ' = {};'.format(lua_pop)
                else:
                    decl_suffix = ';'
                LUA_decl.append(self.std_c_decl(lang, arg, const=arg_const, ptr=ptr) + decl_suffix)
                
                cpp_call_list.append(fmt.cpp_var)

        # call with all arguments
        nargs = len(cpp_call_list)
        len_post_parse = len(post_parse)
        cpp_call_list = ', '.join(cpp_call_list)

        fmt.rv_asgn = fmt.rv_decl + ' = '
        need_rv = False

        if True:  #for nargs, len_post_parse, cpp_call_list in default_calls:
            fmt.cpp_call_list = cpp_call_list
            LUA_code.extend(post_parse[:len_post_parse])

            if is_ctor:
                LUA_code.extend([
                        wformat('{LUA_userdata_type} * {LUA_userdata_var} = ({LUA_userdata_type} *) lua_newuserdata({LUA_state_var}, sizeof(*{LUA_userdata_var}));', fmt),
                        wformat('{LUA_userdata_var}->{LUA_userdata_member} = new {cpp_class}({cpp_call_list});', fmt),
                        '/* Add the metatable to the stack. */',
                        wformat('luaL_getmetatable(L, "{LUA_metadata}");', fmt),
                        '/* Set the metatable on the userdata. */',
                        'lua_setmetatable(L, -2);',
                        ])
            elif is_dtor:
                LUA_code.extend([
                        wformat('delete {LUA_userdata_var}->{LUA_userdata_member};', fmt),
                        wformat('{LUA_userdata_var}->{LUA_userdata_member} = NULL;', fmt),
                        ])
            elif CPP_subprogram == 'subroutine':
                line = wformat('{LUA_this_call}{function_name}({cpp_call_list});', fmt)
                LUA_code.append(line)
            else:
                need_rv = True
                line = wformat('{rv_asgn}{LUA_this_call}{function_name}({cpp_call_list});', fmt)
                LUA_code.append(line)

#            if 'PY_error_pattern' in node:
#                lfmt = util.Options(fmt)
#                lfmt.c_var = fmt.rv
#                lfmt.cpp_var = fmt.rv
#                append_format(LUA_code, self.patterns[node['PY_error_pattern']], lfmt)

        need_rv = False

        if need_rv:
            LUA_decl.append(fmt.rv_decl + ';')
#        if len(LUA_decl):
#            LUA_decl.append('')

        # Compute return value
        if CPP_subprogram == 'function' and not is_ctor:
            fmt.cpp_var = fmt.rv
            fmt.c_var = wformat(result_typedef.cpp_to_c, fmt)  # if C++
##            append_format(LUA_push, result_typedef.LUA_push, fmt)
            tmp = wformat(result_typedef.LUA_push, fmt)
            LUA_push.append( tmp + ';' )
#            format, vargs = self.intent_out(result_typedef, fmt, PY_code)
#            # Add result to front of result tuple
#            build_format.insert(0, format)
#            build_vargs.insert(0, vargs)

#        LUA_code.extend(post_call)

        num_return = len(LUA_push)
        if is_ctor:
            # ctor uses lua_newuserdata which leaves on stack (no push).
            num_return += 1

        body = self.body_lines
        body.extend(LUA_decl)
        body.extend(LUA_code)
        body.extend(LUA_push)    # return values
        body.append('SH_nresult = {};'.format(num_return))

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
                ])
        util.extern_C(output, 'begin')

        if options.cpp_header:
            for include in options.cpp_header.split():
                output.append('#include "%s"' % include)

        output.append('#include "lua.h"')
        output.extend(self.lua_type_structs)
        output.append('')
        output.append(wformat('int luaopen_{LUA_module_name}(lua_State *{LUA_state_var});', fmt))
        output.append('')
        util.extern_C(output, 'end')
        output.append('#endif  /* %s */' % guard)
        self.write_output_file(fname, self.config.python_dir, output)

    def append_luaL_Reg(self, output, name, lines):
        """Create luaL_Reg struct"""
        output.extend([
                '',
                'static const struct luaL_Reg {} [] = {{'.format(name),
                1,
                ])
        output.extend(lines)
        output.extend([
                '{NULL, NULL}   /*sentinel */',
                -1,
                '};',
                ])

    def write_module(self, node):
        options = node['options']
        fmt = node['fmt']
        fname = fmt.LUA_module_filename

        output = []

        if options.cpp_header:
            for include in options.cpp_header.split():
                output.append('#include "{}"'.format(include))
        output.append(wformat('#include "{LUA_header_filename}"', fmt))

        util.extern_C(output, 'begin')
        output.append('#include "lauxlib.h"')
        util.extern_C(output, 'end')

        self._create_splicer('include', output)

        output.append('')
        self.namespace(node, 'begin', output)
        self._create_splicer('C_definition', output)

        output.extend(self.body_lines)

        self.append_luaL_Reg(output, fmt.LUA_module_reg, self.luaL_Reg_module)
        output.append('')
        util.extern_C(output, 'begin')
        output.extend([
                wformat('int luaopen_{LUA_module_name}(lua_State *{LUA_state_var}) {{', fmt),
                1,
                ])
        output.extend(self.class_lines)
        output.extend([
                '',
                wformat('luaL_newlib({LUA_state_var}, {LUA_module_reg});', fmt),
                'return 1;',
                -1,
                '}'
                ])
        util.extern_C(output, 'end')

        self.namespace(node, 'end', output)

        self.write_output_file(fname, self.config.lua_dir, output)


class lua_function(object):
    """Gather information used to write a wrapper for 
    and overloaded/default-argument function
    """
    def __init__(self, function, inargs, outargs):
        self.function = function
        self.nargs = len(inargs)
        self.inargs = inargs
        self.outargs = outargs
