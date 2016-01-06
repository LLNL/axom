import sys
from pybindgen import *

mod = Module('Hydro')
mod.add_include('"HydroC.hpp"')

mod.add_class('PolygonMeshXY', import_from_module='PolygonMeshXY')
mod.add_class('State', import_from_module='State')
mod.add_class('VectorXY', import_from_module='vectorXY')

klass = mod.add_class('Hydro')

# using a pointer in c'tor avoids the State object ever being destroyed
klass.add_constructor([param('State *', 's',transfer_ownership=False)])

klass.add_method('step', None, [param('double', 'dt')])
klass.add_method('steps', None, [param('int', 'numSteps')])
klass.add_method('newDT', retval('double'), [])
klass.add_method('advance', None, [param('double', 'stopTime')])
klass.add_method('setState', None, [param('State','s')])
klass.add_method('getState', retval('State *', is_const=False, reference_existing_object=False), [])
klass.add_method('calcQ', None, [param('State', 'state')])
klass.add_method('getQ', 'double', [param('int', 'i')])
klass.add_method('calcForce', None, [param('const State &', 'state')])
klass.add_method('getForce', 'VectorXY', [param('int', 'i')])
klass.add_method('initialize', None, [])
klass.add_method('setBC', None, [param('char *','boundary'),
                                 param('double','xVel'),
                                 param('double','yVel')])
klass.add_method('numBCnodes', retval('int'), [param('int', 'i')])
klass.add_method('bcNode', retval('int'), [param('int', 'bc'), param('int', 'n')])
klass.add_method('totalEnergy', retval('double'), [])
klass.add_instance_attribute('printCycle', 'int')
klass.add_instance_attribute('cycle', 'int')
klass.add_instance_attribute('time', 'double')
klass.add_instance_attribute('Cq', 'double')
klass.add_instance_attribute('Cl', 'double')
klass.add_instance_attribute('cfl', 'double')
# can't quite get state as an attribute yet
# klass.add_instance_attribute('state', retval('State *', is_const=False, reference_existing_object=True),
#                              # 'State *',
#                              is_const=False, getter='getState')

mod.generate(sys.stdout)
