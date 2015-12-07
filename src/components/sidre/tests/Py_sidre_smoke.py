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
