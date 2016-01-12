import sys
from  pybindgen import *

mod = Module('PolygonMeshXY')
mod.add_struct('VectorXY', import_from_module='vectorXY')
mod.add_include('"PolygonMeshXY.hpp"')
klass = mod.add_class('PolygonMeshXY')

klass.add_constructor([param('int', 'kmax'),
                       param('int', 'lmax'),
                       param('double', 'xmin', default_value = '0.0'),
                       param('double', 'xmax', default_value = '1.0'),
                       param('double', 'ymin', default_value = '0.0'),
                       param('double', 'ymax', default_value = '1.0')])

klass.add_method('getPos', retval('VectorXY'), [param('int', 'i')])

klass.add_method('setPos', None, [param('int', 'i'), param('VectorXY', 'v'),])

klass.add_method('getZonePos', retval('VectorXY'), [param('int', 'i')])
klass.add_method('zoneVol', retval('double'), [param('int', 'i')])

klass.add_method('zoneNumNodes', retval('int'), [param('int', 'i')])
klass.add_method('zoneNode', retval('int'), [param('int', 'iz'), param('int', 'in')])
klass.add_method('zoneNodePos', retval('VectorXY'), [param('int', 'iz'), param('int', 'in')])

klass.add_method('meshAverageKLZMemOrderA', retval('VectorXY'), [])

klass.add_instance_attribute('nzones', 'int')
klass.add_instance_attribute('nnodes', 'int')
klass.add_instance_attribute('timeElapsed', 'double')

mod.generate(sys.stdout)
