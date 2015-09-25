import sys
from os.path import expanduser
sys.path.append(expanduser('~')+'/python/PyBindGen-0.17.0')
from pybindgen import *

mod = Module('Part')
mod.add_include('"Part.h"')

klass = mod.add_class('Part')

mod.add_container('std::vector<int>', 'int', 'vector') # declare a container only once

klass.add_constructor([param('std::vector<int>', 'zoneList'),
                       param('double', 'gamma', default_value='5.0/3.0')])

klass.add_method('rho', retval('double'),  [param('int', 'i')])
klass.add_method('e', retval('double'),  [param('int', 'i')])

klass.add_method('setRho', None, [param('int', 'i'), param('double', 'val')])
klass.add_method('setE', None, [param('int', 'i'), param('double', 'val')])

klass.add_method('meshZone', retval('int'), [param('int', 'i')])

klass.add_instance_attribute('gamma', 'double')
klass.add_instance_attribute('nzones', 'int')

# klass.add_method('operator+=', retval('Part'),
#                  [param('Part', 'b')])

# klass.add_method('operator*=', retval('Part'),
#                  [param('double', 's')])

mod.generate(sys.stdout)
