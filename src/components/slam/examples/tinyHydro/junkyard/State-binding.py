import sys
from pybindgen import *

mod = Module('State')
mod.add_include('"State.hpp"')

mod.add_struct('VectorXY', import_from_module='vectorXY')
mod.add_struct('PolygonMeshXY', import_from_module='PolygonMeshXY')
mod.add_struct('Part', import_from_module='Part')

klass = mod.add_class('State')

klass.add_constructor([param('PolygonMeshXY', 'mesh')])

klass.add_method('u', retval('VectorXY'),
                 [param('int', 'i')])

klass.add_method('setU', None,
                 [param('int', 'i'),
                  param('VectorXY', 'val')])

klass.add_method('partBegin', retval('Part *',reference_existing_object=True),[])
klass.add_method('getPart', retval('Part *',reference_existing_object=False),[param('int', 'i')])
klass.add_method('addPart', None, [param('Part *', 'part',transfer_ownership=False)])
klass.add_instance_attribute('nParts', 'int')
# klass.add_method('averageRho', retval('double'), [param('int', 'i')])
# klass.add_method('totalE', retval('double'), [param('int', 'i')])

# klass.add_method('operator+=', retval('State'),
#                  [param('State', 'b')])

# klass.add_method('operator*=', retval('State'),
#                  [param('double', 's')])

mod.generate(sys.stdout)
