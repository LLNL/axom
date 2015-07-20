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
    url      = "https://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz"
    list_url = "https://code.google.com/p/sparsehash/downloads/list/"
    list_depth = 2

    version('2.0.2', '1db92ed7f257d9b5f14a309d75e8a1d4')    

    def install(self, spec, prefix):
        # simply install the headers
        mkdirp(prefix.include)
        shutil.copytree('src/sparsehash', pjoin(prefix.include,"sparsehash"))
