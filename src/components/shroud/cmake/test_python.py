#
# Test if modules required for shroud are available 
#
from __future__ import print_function
import sys

missing = []

try:
    import yaml
except ImportError:
    missing.append('yaml')
try:
    import parsley
except ImportError:
    missing.append('parsley')

if missing:
    print('Missing modules ' + ', '.join(missing))
    sys.exit(1)
    
