from spack import *
import shutil

import os
from os.path import join as pjoin

class SparsehashHeaders(Package):
    """An extremely memory-efficient hash_map implementation. 2 bits/entry
       overhead! The SparseHash library contains several hash-map
       implementations, including implementations that optimize for space or
       speed.

       These hashtable implementations are similar in API to SGI's hash_map
       class and the tr1 unordered_map class, but with different performance
       characteristics. It's easy to replace hash_map or unordered_map by
       sparse_hash_map or dense_hash_map in C++ code. 
    """
    homepage = "https://code.google.com/p/sparsehash"
    url      = "https://github.com/sparsehash/sparsehash/archive/sparsehash-2.0.2.tar.gz"


    version('2.0.2', '3d946fb1eabf2a1990bb50b04de0fa12')    

    def install(self, spec, prefix):
        # simply install the headers
        configure("--prefix=%s" % prefix)
        make("src/sparsehash/internal/sparseconfig.h")
        mkdirp(prefix.include)
        shutil.copytree('src/sparsehash', pjoin(prefix.include,"sparsehash"))
