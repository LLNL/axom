"""
Test utility module
"""
from __future__ import print_function

from shroud import util

import unittest

class UtilCase(unittest.TestCase):

    def test_un_camel(self):
        self.assertEqual(util.un_camel('incrementCount'),
                         'increment_count')
        self.assertEqual(util.un_camel('local_function1'),
                         'local_function1')
        self.assertEqual(util.un_camel('getHTTPResponseCode'),
                         'get_http_response_code')


class OptionCase(unittest.TestCase):

    def setUp(self):
        self.lev0 = util.Options(None, a=1, b=2, c=3)
        self.lev1 = util.Options(self.lev0, x=100, y=1, z=102)

    def test_access01(self):
        # 'a' accessable from both
        self.assertEqual(self.lev0.a, 1)
        self.assertEqual(self.lev1.a, 1)

        # 'z' only accessable from lev1
        with self.assertRaises(AttributeError):
            self.lev0.z
        self.assertEqual(self.lev1.z, 102)

    def test_access02(self):
        """set and access"""
        self.lev0.c2 = 32
        self.assertEqual(self.lev0.c2, 32)

    def test_get01(self):
        self.assertEqual(self.lev1.get('a', 'notfound'), 1)
        self.assertEqual(self.lev1.get('nosuch', 'notfound'), 'notfound')

    def test_in(self):
         self.assertIn('a', self.lev0)
         self.assertIn('a', self.lev1)

         self.assertNotIn('z', self.lev0)
         self.assertIn('z', self.lev1)

         self.assertNotIn('nosuch', self.lev1)

    def test_setdefault(self):
        lev1 = self.lev1
        self.assertNotIn('yyy', lev1)
        lev1.setdefault('yyy', 'yyyvalue')
        self.assertIn('yyy', lev1)
        self.assertEqual(lev1.yyy, 'yyyvalue')


if __name__ == '__main__':
    unittest.main()
