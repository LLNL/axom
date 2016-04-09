"""
Generate Lua module for C++ code.
"""
from __future__ import print_function

import util
from util import wformat, append_format

def add_templates(options):
    options.update(dict(
        LUA_module_filename_template = 'lua{library}module.cpp',
#        LUA_header_filename_template = 'lua{library}module.hpp',
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
        self.comment = '--'

    def reset_file(self):
        pass

    def wrap_library(self):
        top = self.tree
        options = self.tree['options']
        fmt_library = self.tree['fmt']

        # Format variables
        fmt_library.LUA_prefix          = options.get('LUA_prefix', 'l_')
        fmt_library.LUA_module_name     = fmt_library.library_lower
        util.eval_template(top, 'LUA_module_filename')

        # Variables to accumulate output lines
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

#        self.write_header(self.tree)
        self.write_module(self.tree)
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

##-        result = node['result']
##-        result_type = result['type']
##-        result_is_ptr = result['attrs'].get('ptr', False)
##-        result_is_ref = result['attrs'].get('reference', False)
##-
##-        if node.get('return_this', False):
##-            result_type = 'void'
##-            result_is_ptr = False
##-            CPP_subprogram = 'subroutine'
##-
##-        result_typedef = self.typedef[result_type]
##-        is_ctor  = node['attrs'].get('constructor', False)
##-        is_dtor  = node['attrs'].get('destructor', False)
##-#        is_const = result['attrs'].get('const', False)
##-        if is_ctor:   # or is_dtor:
##-            # XXX - have explicit delete
##-            # need code in __init__ and __del__
##-            return




        body = self.body_lines
        body.extend([
                '',
                wformat('static int {LUA_name_impl}(lua_State *L)', fmt),
                '{',
                1,
                'return 0;',
                -1,
                '}'
                ])

    def write_module(self, node):
        options = node['options']
        fmt = node['fmt']
        fname = fmt.LUA_module_filename

        output = []

        output.extend(self.body_lines)

        self.write_output_file(fname, self.config.lua_dir, output)




        
