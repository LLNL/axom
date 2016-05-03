# ${CMAKE_COMMAND} -E create_symlink
#
# test the tutorial module
#

import unittest
import tutorial

class NotTrue:
    """Test bool arguments errors"""
    def __bool__(self):
        raise NotImplementedError
 
class Tutorial(unittest.TestCase):
    """Test tutorial problem"""
     
    def XXsetUp(self):
        """ Setting up for the test """
        print "FooTest:setUp_:begin"
        ## do something...
        print "FooTest:setUp_:end"
     
    def XXtearDown(self):
        """Cleaning up after the test"""
        print "FooTest:tearDown_:begin"
        ## do something...
        print "FooTest:tearDown_:end"
     
    # test routine A
    def testFunction1(self):
        tutorial.Function1()

    def testFunction2(self):
        rv_double = tutorial.Function2(1.0, 4)
        self.assertEqual(rv_double, 5.0)#, "A is not equal to B")

    def testFunction3(self):
        rv_logical = tutorial.Function3(False)
        self.assertTrue(rv_logical)

        self.assertRaises(TypeError, tutorial.Function3, 0)
#rv_logical = tutorial.Function3(NotTrue())

    def testFunction4a(self):
        rv_char = tutorial.Function4a("dog", "cat")
        self.assertEqual(rv_char, "dogcat")
#
#    call function4b("dog", "cat", rv_char)
#    call assert_true( rv_char == "dogcat")
#
    def testFunction5(self):
        rv_double = tutorial.Function5()
        self.assertAlmostEqual(rv_double, 13.1415)
        rv_double = tutorial.Function5(1.0)
        self.assertAlmostEqual(rv_double, 11.0)

        # XXX fix bool argument
        rv_double = tutorial.Function5(1.0, False)
        self.assertAlmostEqual(rv_double, 1.0)
#, 1.d0)
#
    def testFunction6(self):
        tutorial.Function6("name")
        self.assertEqual(tutorial.LastFunctionCalled(), "Function6(string)")

        tutorial.Function6(1)
        self.assertEqual(tutorial.LastFunctionCalled(), "Function6(int)")

        self.assertRaises(TypeError, tutorial.Function6, 1.0)


#
#    call function7(1)
#    call assert_true(last_function_called() == "Function7<int>")
#    call function7(10.d0)
#    call assert_true(last_function_called() == "Function7<double>")
#
#    ! return values set by calls to function7
#    rv_integer = function8_int()
#    call assert_true(rv_integer == 1)
#    rv_double = function8_double()
#    call assert_true(rv_double == 10.d0)
#
#    call function9(1.0)
#    call assert_true(.true.)
#    call function9(1.d0)
#    call assert_true(.true.)
#
#    call function10()
#    call assert_true(.true.)
#    call function10("foo", 1.0e0)
#    call assert_true(.true.)
#    call function10("bar", 2.0d0)
#    call assert_true(.true.)
#
#    call sum(5, [1,2,3,4,5], rv_int)
#    call assert_true(rv_int .eq. 15)
#
#    rv_int = overload1(10)
#    call assert_true(rv_int .eq. 10)
#    rv_int = overload1(1d0, 10)
#    call assert_true(rv_int .eq. 10)
#
#    rv_int = overload1(10, 11, 12)
#    call assert_true(rv_int .eq. 142)
#    rv_int = overload1(1d0, 10, 11, 12)
#    call assert_true(rv_int .eq. 142)
#
#    rv_int = typefunc(2)
#    call assert_true(rv_int .eq. 2)
#
#    rv_int = enumfunc(1)
#    call assert_true(rv_int .eq. 2)
#
#  end subroutine test_functions
#
#  subroutine test_class1
#    type(class1) obj
#
#    obj = class1_new()
#    call assert_true(c_associated(obj%voidptr), "class1_new")
#
#    call obj%method1()
#    call assert_true(.true.)
#
#    call obj%delete()
#    call assert_true(.not. c_associated(obj%voidptr), "class1_delete")
#  end subroutine test_class1
#

# creating a new test suite
newSuite = unittest.TestSuite()
 
# adding a test case
newSuite.addTest(unittest.makeSuite(Tutorial))

if __name__ == "__main__":
    unittest.main()
