###############################################################################
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
###############################################################################

import sidre
ds = sidre.DataStore()

print ds
print dir(ds)

grp = ds.getRoot()
print "Group", grp
#print "Dir of Group", dir(grp)

print "Name of Group", grp.getName()

print "Name of non-existent view", grp.getViewName(100)



#xx = sidre.DataGroup()
#print xx

