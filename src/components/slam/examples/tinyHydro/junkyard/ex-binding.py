import sys
import pybindgen

mod = pybindgen.Module('ex')
mod.add_include('"example.h"')
mod.add_function('fact', pybindgen.retval('int'),
                 [pybindgen.param('int', 'n')])
struct = mod.add_struct('Foo')
struct.add_instance_attribute('a', 'int')
struct.add_instance_attribute('b', 'double')

from pybindgen import *
klass = mod.add_class('Darray')
klass.add_constructor([param('int','size',default_value='10')])
klass.add_method('get', retval('double'), [param('int', 'i')])
klass.add_method('set', None, [param('int', 'i'),
                               param('double', 'x'),])
klass.add_method('size', retval('int'), [])

mod.generate(sys.stdout)
