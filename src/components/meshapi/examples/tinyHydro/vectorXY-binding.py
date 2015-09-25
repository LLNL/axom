import sys
from os.path import expanduser
sys.path.append(expanduser('~')+'/python/PyBindGen-0.17.0')
import pybindgen

mod = pybindgen.Module('vectorXY')
mod.add_include('"VectorXY.h"')
klass = mod.add_class('VectorXY')
klass.add_instance_attribute('x', 'double')
klass.add_instance_attribute('y', 'double')

klass.add_constructor([])

klass.add_constructor([pybindgen.param('double', 'x'),
                       pybindgen.param('double', 'y')])

klass.add_method('add', pybindgen.retval('VectorXY'),
                 [pybindgen.param('VectorXY', 'a')])

klass.add_method('accum', None,
                 [pybindgen.param('VectorXY', 'b')])

klass.add_method('sub', pybindgen.retval('VectorXY'),
                 [pybindgen.param('VectorXY', 'b')])

klass.add_method('elim', None,
                 [pybindgen.param('VectorXY', 'b')])

klass.add_method('mult', pybindgen.retval('VectorXY'),
                 [pybindgen.param('double', 's')])

klass.add_method('scale', None,
                 [pybindgen.param('double', 's')])

klass.add_method('cross', pybindgen.retval('double'),
                 [pybindgen.param('VectorXY', 'v')])

klass.add_method('mag2', pybindgen.retval('double'), [])
klass.add_method('mag' , pybindgen.retval('double'), [])

mod.generate(sys.stdout)
