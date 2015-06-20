#!/bin/env python3
"""
Generate Python module for C++ code.

Entire library in a single header.
One Extension module per class

"""
from __future__ import print_function

import util
import fwrap_util

import os

wformat = util.wformat

class Wrapp(object):
    """Generate Python bindings.
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
        self.f_type_decl = []
        self.c_interface = []
        self.impl = []         # implementation, after contains
        self.c_interface.append('interface')
        self.c_interface.append(1)

    def _end_output_file(self):
        self.c_interface.append(-1)
        self.c_interface.append('end interface')

    def _begin_class(self):
        pass

    def reset_file(self):
        self.PyMethodBody = []
        self.PyMethodDef = []

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
#        out.append('! splicer push %s' % name)

#X changes for push/pop instead of full paths
    def _pop_splicer(self, name, out):
        # XXX maybe use name for error checking, must pop in reverse order
        self.splicer_stack.pop()
        self.splicer_names.pop()
        if self.splicer_names:
            self.splicer_path = '.'.join(self.splicer_names) + '.'
        else:
            self.splicer_path = ''
#X        out.append('! splicer pop %s' % name)

    def _create_splicer(self, name, out, default=None):
        # The prefix is needed when two different sets of output are being create
        # and they are not in sync.
        # Creating methods and derived types together.
#X        out.append('! splicer begin %s' % name)
        out.append('// splicer begin %s%s' % (self.splicer_path, name))
        if default:
            out.extend(default)
        else:
            out.extend(self.splicer_stack[-1].get(name, []))
#X        out.append('// splicer end %s' % name)
        out.append('// splicer end %s%s' % (self.splicer_path, name))

#####

    def _c_type(self, arg):
        """
        Return the Fortran type, and array attribute
        pass-by-value default

        attributes:
        ptr - True = pass-by-reference
        reference - True = pass-by-reference

        """
        t = []
        typedef = self.typedef.get(arg['type'], None)
        if typedef is None:
            raise RuntimeError("No such type %s" % arg['type'])
        is_ptr = (arg['attrs'].get('ptr', False) or
                  arg['attrs'].get('reference', False))
        intent = arg['attrs'].get('intent', None)
        is_value = arg['attrs'].get('value', False)

        typ = typedef.c_fortran
        if typedef.base == 'string':
            return (typ, True)   # is array
        else:
            #        if arg['attrs'].get('const', False):
            #            t.append('const')
            t.append(typ)
            if is_value:
                t.append(', value')
            if intent:
                t.append(', intent(%s)' % intent.upper())
            return (''.join(t), arg['attrs'].get('array', False))

    def _c_decl(self, arg, name=None):
        """
        Return the Fortran declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ, arr = self._c_type(arg)
        rv = typ + ' :: ' + ( name or arg['name'] )
        if arr:
            rv += '(*)'
        return rv

    def _f_type(self, arg):
        """
        Return the Fortran type, and array attribute
        pass-by-value default

        attributes:
        ptr - True = pass-by-reference
        reference - True = pass-by-reference

        """
        t = []
        typedef = self.typedef.get(arg['type'], None)
        if typedef is None:
            raise RuntimeError("No such type %s" % arg['type'])

        typ = typedef.f_type
        if typedef.base == 'string':
            return (typ, False)  # not array
        else:
            #        if arg['attrs'].get('const', False):
            #            t.append('const')
            t.append(typ)
#            if not (arg['attrs'].get('ptr', False) or
#                    arg['attrs'].get('reference', False)):
#                t.append(', value')
            return (''.join(t), arg['attrs'].get('array', False))

    def _f_decl(self, arg, name=None):
        """
        Return the Fortran declaration.

        If name is not supplied, use name in arg.
        This makes it easy to reproduce the arguments.
        """
        typ, arr = self._f_type(arg)
        rv = typ + ' :: ' + ( name or arg['name'] )
        if arr:
            rv += '(*)'
        return rv

    def wrap_library(self):
        options = self.tree['options']
        fmt_library = self.tree['fmt']

        fmt_library.PY_module_filename = 'python_module.cpp'
        fmt_library.PY_header_filename = 'python_header.cpp'
        self.py_type_object_creation = []
        self.py_type_structs = []

        self._push_splicer('class', [])
        for node in self.tree['classes']:
            name = node['name']
            self.reset_file()
            self.wrap_class(node)
            self.write_extension_type(node)
        self._pop_splicer('class', [])

        self.reset_file()
        if self.tree['functions']:
            self._push_splicer('function', [])
            self._begin_class()
            for node in self.tree['functions']:
                self.wrap_method(None, node)
            self._pop_splicer('function', [])

        self.write_header(self.tree)
        self.write_module(self.tree)

    def wrap_class(self, node):
        self.log.write("class {1[name]}\n".format(self, node))
        name = node['name']
        unname = util.un_camel(name)
        typedef = self.typedef[name]

        options = node['options']
        fmt_class = node['fmt']

        fmt_class.PY_type_filename = wformat('python_{lower_class}.cpp', fmt_class)

        self.py_type_object_creation.append(wformat("""
// {cpp_class}
    PB_Dbnode_Type.tp_new   = PyType_GenericNew;
    PB_Dbnode_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&PB_Dbnode_Type) < 0)
        return RETVAL;
    Py_INCREF(&PB_Dbnode_Type);
    PyModule_AddObject(m, "Dbnode", (PyObject *)&PB_Dbnode_Type);
""", fmt_class))


        self._create_splicer('C_declaration', self.py_type_structs)
        self.py_type_structs.append('typedef struct {')
        self.py_type_structs.append('PyObject_HEAD')
        self._create_splicer('C_object', self.py_type_structs)
        self.py_type_structs.append('} PB_PackageObject;')

        # wrap methods
        self._push_splicer(name, [])
        self._push_splicer('method', [])
        for method in node['methods']:
            self.wrap_method(node, method)
        self._pop_splicer('method', [])

        self._pop_splicer(fmt_class.lower_class, [])

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

        options = node['options']
        fmt_func = node['fmt']
        fmt_func.doc_string = 'documentation'
        fmt_func.PY_func_name    = 'AAA'

        PY_code = []


        self.PyMethodBody.append(wformat("""
static char PB_stop_here__doc__[] =
"{doc_string}"
;

static PyObject *
{PY_func_name}(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{{
""", fmt_func))
        self._create_splicer(fmt_func.F_name_method, self.PyMethodBody, PY_code)
        self.PyMethodBody.append('}')
                                 

        self.PyMethodDef.append( wformat('{{"find_package", (PyCFunction){PY_func_name}, METH_VARARGS|METH_KEYWORDS, PB_find_package__doc__}},', fmt_func))


        return
    #######################
        result = node['result']
        result_type = result['type']
        result_is_ptr = result['attrs'].get('ptr', False)

        if node.get('return_this', False):
            result_type = 'void'
            result_is_ptr = False

        result_typedef = self.typedef[result_type]
        is_ctor  = result['attrs'].get('constructor', False)
        is_dtor  = result['attrs'].get('destructor', False)
        is_const = result['attrs'].get('const', False)

        fmt_func.F_obj = wformat('{F_this}%{F_this}', fmt_func)
        if 'F_this' in options:
            fmt_func.F_this = options.F_this
        if 'F_result' in options:
            fmt_func.F_result = options.F_result
        F_result = fmt_func.F_result
        C_this = fmt_func.C_this

        if 'F_C_name' in options:
            fmt_func.F_C_name = options.F_C_name
        else:
            fmt_func.F_C_name = fmt_func.C_name.lower()

        if cls:
            util.eval_template(options, fmt_func,
                               'F_name_impl', '{lower_class}_{underscore_name}{method_suffix}')
        else:
            util.eval_template(options, fmt_func,
                               'F_name_impl', '{underscore_name}{method_suffix}')
        util.eval_template(options, fmt_func,
                            'F_name_method', '{underscore_name}{method_suffix}')
        util.eval_template(options, fmt_func,
                            'F_name_generic', '{underscore_name}')

        arg_c_names = [ ]
        arg_c_decl = [ ]
        arg_c_call = []      # arguments to C function

        arg_f_names = [ ]
        arg_f_decl = [ ]
        arg_f_use  = [ 'use iso_c_binding' ]  # XXX totally brain dead for now

        # find subprogram type
        # compute first to get order of arguments correct.
        # Add 
        if result_type == 'void' and not result_is_ptr:
            #  void=subroutine   void *=function
            subprogram = 'subroutine'
        else:
            subprogram = 'function'
            fmt_func.F_result_clause = ' result(%s)' % F_result
            if is_const:
                fmt_func.F_pure_clause   = 'pure '
        fmt_func.F_subprogram    = subprogram

        if cls:
            # Add 'this' argument
            if not is_ctor:
                arg_c_names.append(C_this)
                arg_c_decl.append('type(C_PTR), value, intent(IN) :: ' + C_this)
                arg_c_call.append(fmt_func.F_obj)
                arg_f_names.append(fmt_func.F_this)
                if is_dtor:
                    arg_f_decl.append(wformat(
                            'type({F_derived_name}) :: {F_this}',
                            fmt_func))
                else:
                    arg_f_decl.append(wformat(
                            'class({F_derived_name}) :: {F_this}',
                            fmt_func))

        for arg in node.get('args', []):
            # default argument's intent
            # XXX look at const, ptr
            attrs = arg['attrs']
            if 'intent' not in attrs:
                attrs['intent'] = 'in'
            if 'value' not in attrs:
                attrs['value'] = True

            arg_c_names.append(arg['name'])
            arg_c_decl.append(self._c_decl(arg))

            rrr = self.typedef[arg['type']].fortran_to_c
            arg_c_call.append(rrr.format(var=arg['name']))
            arg_f_names.append(arg['name'])
            arg_f_decl.append(self._f_decl(arg))

        # declare function return value after arguments
        # since arguments may be used to compute return value
        # (for example, string lengths)
        if subprogram == 'function':
            if result_typedef.base == 'string':
                # special case returning a string
                rvlen = result['attrs'].get('len', '1')
                fmt_func.rvlen = rvlen
                arg_c_decl.append('type(C_PTR) %s' % F_result)
                arg_f_decl.append(
                    wformat('character(kind=C_CHAR, len={rvlen}) :: {F_result}',
                            fmt_func))
            else:
                # XXX - make sure ptr is set to avoid VALUE
                arg_dict = dict(name=F_result,
                                type=result_type,
                                attrs=dict(ptr=True))
#                arg_c_decl.append(self._c_decl(result, name=F_result))
                arg_c_decl.append(self._c_decl(arg_dict))
                arg_f_decl.append(self._f_decl(result, name=F_result))

        if not is_ctor and not is_dtor:
            # Add method to derived type
            F_name_method = fmt_func.F_name_method
            self.f_type_generic.setdefault(fmt_func.F_name_generic,[]).append(F_name_method)
            self.type_bound_part.append('procedure :: %s => %s' % (
                    F_name_method, fmt_func.F_name_impl))

        fmt_func.F_arg_c_call = ', '.join(arg_c_call)
        fmt_func.F_C_arguments = options.get('F_C_arguments', ', '.join(arg_c_names))
        fmt_func.F_arguments = options.get('F_arguments', ', '.join(arg_f_names))

        # body of function
        splicer_code = self.splicer_stack[-1].get(fmt_func.F_name_method, None)
        if 'F_code' in options:
            F_code = [ options.F_code ]
        elif splicer_code:
            F_code = splicer_code
        else:
            F_code = []
            if is_ctor:
                F_code.append(wformat('{F_result}%{F_this} = {F_C_name}({F_arg_c_call})', fmt_func))
            elif subprogram == 'function':
                fmt = result_typedef.f_return_code
                F_code.append(wformat(fmt, fmt_func))
            else:
                F_code.append(wformat('call {F_C_name}({F_arg_c_call})', fmt_func))
            if is_dtor:
                F_code.append(wformat('{F_this}%{F_this} = C_NULL_PTR', fmt_func))

        c_interface = self.c_interface
        c_interface.append('')
        c_interface.append(wformat(
                '{F_pure_clause}{F_subprogram} {F_C_name}({F_C_arguments}){F_result_clause} &',
                fmt_func))
        c_interface.append(2)  # extra indent for continued line
        c_interface.append(wformat(
                'bind(C, name="{C_name}")',
                fmt_func))
        c_interface.append(-1)
        c_interface.append('use iso_c_binding')
        c_interface.append('implicit none')
        c_interface.extend(arg_c_decl)
        c_interface.append(-1)
        c_interface.append(wformat('end {F_subprogram} {F_C_name}', fmt_func))

        impl = self.impl
        impl.append('')
        impl.append(wformat('{F_subprogram} {F_name_impl}({F_arguments}){F_result_clause}', fmt_func))
        impl.append(1)
        impl.extend(arg_f_use)
        impl.append('implicit none')
        impl.extend(arg_f_decl)
        self._create_splicer(fmt_func.F_name_method, impl, F_code)
        impl.append(-1)
        impl.append(wformat('end {F_subprogram} {F_name_impl}', fmt_func))

    def write_extension_type(self, node):
        fmt = node['fmt']
        fname = fmt.PY_type_filename

        output = []

        self._create_splicer('extra_methods', output)

        output.extend(self.PyMethodBody)

        output.append('static PyMethodDef PB_methods[] = {')
        output.extend(self.PyMethodDef)
        output.append('{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */')
        output.append('};')

        output.append(wformat(PyTypeObject_template, fmt))

        self.write_output_file(fname, self.config.binary_dir, output)

    def write_header(self, node):
        options = node['options']
        fmt = node['fmt']
        fname = fmt.PY_header_filename

        output = []

        output.append("""
/*
 * This is generated code.
 * Any edits must be made between the splicer.begin and splicer.end blocks.
 * All other edits will be lost.
 * Once a block is edited remove the 'UNMODIFIED' on the splicer.begin
 * comment to allow the block to be preserved when it is regenerated.
 */

#ifndef HDR_BASISMODULE
#define HDR_BASISMODULE
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif""")
        
        self._create_splicer('C_declaration', output)

        output.extend(self.py_type_structs)

        output.append("""
extern PyObject *PB_error_obj;

#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_basis
#else
#define MOD_INITBASIS initbasis
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#endif
""")


        self.write_output_file(fname, self.config.binary_dir, output)

    def write_module(self, node):
        options = node['options']
        fmt = node['fmt']
        fname = fmt.PY_module_filename

        fmt.PY_library_doc = 'library documentation'

        output = []
        self._create_splicer('extra_methods', output)
        output.extend(self.PyMethodBody)

        output.append('static PyMethodDef PB_methods[] = {')
        output.extend(self.PyMethodDef)
        output.append('{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */')
        output.append('};')

        


        output.append(wformat(module_begin, fmt))
        self._create_splicer('C_init_locals', output)
        output.append(wformat(module_middle, fmt))
        output.extend(self.py_type_object_creation)
        output.append(wformat(module_middle2, fmt))
        self._create_splicer('C_init_body', output)
        output.append(wformat(module_end, fmt))

        self.write_output_file(fname, self.config.binary_dir, output)

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



#### Python boiler plate


PyTypeObject_template = """
static char PB_Dbnode__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PB_Dbnode_Type = {{
        PyVarObject_HEAD_INIT(NULL, 0)
        "basis.Dbnode",                       /* tp_name */
        sizeof(PB_DbnodeObject),         /* tp_basicsize */
        0,                              /* tp_itemsize */
        /* Methods to implement standard operations */
        (destructor)0,                 /* tp_dealloc */
        (printfunc)0,                   /* tp_print */
        (getattrfunc)0,                 /* tp_getattr */
        (setattrfunc)0,                 /* tp_setattr */
#ifdef IS_PY3K
        0,                               /* tp_reserved */
#else
        (cmpfunc)0,                     /* tp_compare */
#endif
        (reprfunc)PB_Dbnode_repr,                    /* tp_repr */
        /* Method suites for standard classes */
        0,                              /* tp_as_number */
        0,                              /* tp_as_sequence */
        0,                              /* tp_as_mapping */
        /* More standard operations (here for binary compatibility) */
        (hashfunc)0,                    /* tp_hash */
        (ternaryfunc)0,                 /* tp_call */
        (reprfunc)0,                    /* tp_str */
        (getattrofunc)0,                /* tp_getattro */
        (setattrofunc)0,                /* tp_setattro */
        /* Functions to access object as input/output buffer */
        0,                              /* tp_as_buffer */
        /* Flags to define presence of optional/expanded features */
        Py_TPFLAGS_DEFAULT,             /* tp_flags */
        PB_Dbnode__doc__,         /* tp_doc */
        /* Assigned meaning in release 2.0 */
        /* call function for all accessible objects */
        (traverseproc)0,                /* tp_traverse */
        /* delete references to contained objects */
        (inquiry)0,                     /* tp_clear */
        /* Assigned meaning in release 2.1 */
        /* rich comparisons */
        (richcmpfunc)0,                 /* tp_richcompare */
        /* weak reference enabler */
        0,                              /* tp_weaklistoffset */
        /* Added in release 2.2 */
        /* Iterators */
        (getiterfunc)0,                 /* tp_iter */
        (iternextfunc)0,                /* tp_iternext */
        /* Attribute descriptor and subclassing stuff */
        0,                             /* tp_methods */
        0,                              /* tp_members */
        PB_Dbnode_getset,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)PB_Dbnode_init,                   /* tp_init */
        (allocfunc)0,                  /* tp_alloc */
        (newfunc)0,                    /* tp_new */
        (freefunc)0,                   /* tp_free */
        (inquiry)0,                     /* tp_is_gc */
        0,                              /* tp_bases */
        0,                              /* tp_mro */
        0,                              /* tp_cache */
        0,                              /* tp_subclasses */
        0,                              /* tp_weaklist */
        (destructor)0,                 /* tp_del */
        0,                              /* tp_version_tag */
#ifdef IS_PY3K
        (destructor)0,                  /* tp_finalize */
#endif
}};
"""


module_begin = """
/*
 * init{lower_library} - Initialization function for the module
 * *must* be called init{lower_library}
 */
static char PB__doc__[] =
"{PY_library_doc}"
;

struct module_state {{
    PyObject *error;
}};

#ifdef IS_PY3K
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

#ifdef IS_PY3K
static int basis_traverse(PyObject *m, visitproc visit, void *arg) {{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}}

static int basis_clear(PyObject *m) {{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}}

static struct PyModuleDef moduledef = {{
    PyModuleDef_HEAD_INIT,
    "{lower_library}", /* m_name */
    PB__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    PB_methods, /* m_methods */
    NULL, /* m_reload */
    {lower_library}_traverse, /* m_traverse */
    {lower_library}_clear, /* m_clear */
    NULL  /* m_free */
}};

#define RETVAL m
#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_{lower_library}(void)

#else
#define RETVAL
#define INITERROR return

PyMODINIT_FUNC
init{lower_library}(void)
#endif
{{
    PyObject *m = NULL;
"""

module_middle = """

    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("{lower_library}", PB_methods,
                       PB__doc__,
                       (PyObject*)NULL,PYTHON_API_VERSION);
#endif
    if (m == NULL)
        return RETVAL;
    struct module_state *st = GETSTATE(m);
"""

module_middle2 = """
    PB_error_obj = PyErr_NewException("{lower_library}.Error", NULL, NULL);
    if (PB_error_obj == NULL)
        return RETVAL;
    st->error = PB_error_obj;
    PyModule_AddObject(m, "Error", st->error);
"""

module_end = """
    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module {lower_library}");
    return RETVAL;
}}
"""
