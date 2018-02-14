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

import unittest
import sidre

class SidreSmoke(unittest.TestCase):

  def test_create_datastore(self):
      ds = sidre.DataStore()
      ds.delete()
      self.assertTrue( True )

  def XXtest_valid_invalid(self):
      ds = sidre.DataStore()

      idx = 3

      self.assertTrue(idx != InvalidIndex)

      name = "foo"
      self.assertTrue(nameIsValid(name))

      root = ds.getRoot()

      self.assertTrue(root.getGroupName(idx) == InvalidName)
      self.assertTrue(root.getGroupIndex(name) == InvalidIndex)

      ds.delete()


if __name__ == '__main__':
    unittest.main()
