#!/usr/bin/env python3

# gen-multidom-structured-mesh.py
# Write a simple multidomain structured blueprint mesh for testing.

# This script requires a conduit installation configured with python3 and hdf5.
# Make sure PYTHONPATH includes /path/to/conduit/install/python-modules

try:
  import conduit
  import conduit.blueprint
  import conduit.relay
except ModuleNotFoundError as e:
  print(f'{e}\nMake sure your PYTHONPATH includes /path/to/conduit/install/python-modules\nConduit must be configured with python and hdf5.')
  exit(-1)

import numpy as np

def i_c(s):
  '''Integer coordinates.'''
  return list(map(int, s.split(',')))
def f_c(s):
  '''Floating point coordinates.'''
  return list(map(float, s.split(',')))

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
ps = ArgumentParser(description='Write a blueprint multidomain unstructured mesh.',
                    formatter_class=ArgumentDefaultsHelpFormatter)
ps.add_argument('--useList', action='store_true', help='Put domains in a list instead of a map')
ps.add_argument('-ml', type=f_c, default=(0.,0.), help='Mesh lower coordinates')
ps.add_argument('-mu', type=f_c, default=(1.,1.), help='Mesh upper coordinates')
ps.add_argument('-ms', type=i_c, default=(3,3), help='Logical size of mesh (cells)')
ps.add_argument('-dc', type=i_c, default=(1,1), help='Domain counts in each index direction')
ps.add_argument('-o', '--output', type=str, default='mdmesh', help='Output file base name')
opts,unkn = ps.parse_known_args()
print(opts, unkn)
if(unkn):
  print("Unrecognized arguments:", *unkn)
  quit(1)

dim = len(opts.dc)

if dim not in (2,3)	\
   or len(opts.ms) != dim	\
   or len(opts.ml) != dim	\
   or len(opts.mu) != dim:
  raise RuntimeError('dc, ms, ml and mu options must have the same dimensions (2 or 3)')

def scale_structured_domain(n, startCoord, endCoord):
  '''This function scales and shift a blueprint structured domain after
  it has been created.  There's no way to specify the physical extent
  of a domain using conduit.blueprint.mesh.examples.basic, as far as I
  can tell.

  '''
  domPhysicalSize = np.array(endCoord) - np.array(startCoord)
  xyz = 'xyz'
  assert(n['topologies/mesh/type'] == 'structured')
  ndim = n['coordsets/coords/values'].number_of_children()
  assert(len(startCoord) >= ndim)
  assert(len(domPhysicalSize) >= ndim)
  for d in range(ndim):
    coords = n['coordsets/coords/values'][d]
    minC, maxC = min(coords), max(coords)
    curRange = maxC - minC
    shift = startCoord[d] - minC
    coords = (coords - minC) * domPhysicalSize[d]/curRange + startCoord[d]
    n['coordsets/coords/values'][xyz[d]] = coords

domType = 'structured'

domCounts = opts.dc if dim == 3 else (*opts.dc, 1)
meshSize = opts.ms if dim == 3 else (*opts.ms, 1)
meshLower = opts.ml if dim == 3 else (*opts.ml, 0)
meshUpper = opts.mu if dim == 3 else (*opts.mu, 0)

# Convert to np.array to use element-wise arithmetic.
domCounts = np.array(domCounts, dtype=np.int)
meshSize = np.array(meshSize, dtype=np.int)
meshLower = np.array(meshLower)
meshUpper = np.array(meshUpper)

domPhysicalSize = (meshUpper - meshLower)/domCounts
cellPhysicalSize = (meshUpper - meshLower)/meshSize

domSize = meshSize//domCounts
domSizeRem = meshSize % domCounts
print(f'meshSize={meshSize} cells, domCounts={domCounts} domSize={domSize} domSizeRem={domSizeRem}')

def domain_size(di, dj, dk):
  rval = np.array(domSize)
  rval += (di,dj,dk) < domSizeRem
  rval[0:dim] += 1
  return rval

def domain_index_begin(di, dj, dk):
  '''Compute first cell index of domain (di, dj, dk).'''
  idx = np.array( (di, dj, dk) )
  std = np.array(domSize) * (di, dj, dk)
  extra = np.where( idx < domSizeRem, idx, domSizeRem )
  begin = std + extra
  #print(di, dj, dk, f'begin={begin}, std={std}, extra={extra}')
  return begin

mdMesh = conduit.Node()
for dk in range(domCounts[2]):
  for dj in range(domCounts[1]):
    for di in range(domCounts[0]):
      if opts.useList:
        dom = mdMesh.append()
      else:
        domName = f'domain_{di:1d}_{dj:1d}' if len(opts.dc) == 2 else f'domain_{di:1d}_{dj:1d}_{dk:1d}'
        dom = mdMesh[domName]

      cellStart = domain_index_begin(di, dj, dk)
      cellEnd = domain_index_begin(di+1, dj+1, dk+1)
      pointCounts = cellEnd - cellStart + (1, 1, 1 if dim == 3 else 0)
      #print(di, dj, dk, f'{cellStart} -> {cellEnd}, {pointCounts}')

      conduit.blueprint.mesh.examples.basic(domType, *pointCounts, dom)

      domLower = meshLower + cellStart * cellPhysicalSize
      domUpper = meshLower + cellEnd * cellPhysicalSize
      scale_structured_domain(dom, domLower, domUpper)

# print(mdMesh)

info = conduit.Node()
if not conduit.blueprint.mesh.verify(mdMesh, info):
    print(info)

conduit.relay.io.blueprint.save_mesh(mdMesh,opts.output, "hdf5")
print(f'Wrote mesh {opts.output}')
