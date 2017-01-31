#
# Shroud tests
#
import test_parse_decl

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(test_parse_decl.CheckDeclCase))
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
