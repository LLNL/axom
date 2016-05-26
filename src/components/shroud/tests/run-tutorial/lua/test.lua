-- test tutorial module

local tutorial = require "tutorial"
local rv_int, rv_double, rv_logical, rv_char

tutorial.Function1()
print(tutorial.LastFunctionCalled())

rv_double = tutorial.Function2(1.0, 4)
print(tutorial.LastFunctionCalled(), rv_double)


rv_logical = tutorial.Function3(false)
print(tutorial.LastFunctionCalled(), rv_logical)

rv_char = tutorial.Function4a("dog", "cat")
print(tutorial.LastFunctionCalled(), rv_char)

rv_char = tutorial.Function4b("dog", "cat")
print(tutorial.LastFunctionCalled(), rv_char)

rv_double = tutorial.Function5()
-- 13.1415
print(tutorial.LastFunctionCalled(), rv_double)
rv_double = tutorial.Function5(11.0)
print(tutorial.LastFunctionCalled(), rv_double)
-- 11.0
rv_double = tutorial.Function5(11.0, false)
print(tutorial.LastFunctionCalled(), rv_double)
-- 1.0

tutorial.Function6("name")
print(tutorial.LastFunctionCalled())
tutorial.Function6(1)
print(tutorial.LastFunctionCalled())

--[[
tutorial.Function7(1)
print(tutorial.LastFunctionCalled())
tutorial.Function7(10.0)
print(tutorial.LastFunctionCalled())
--]]


tutorial.Function10()
print(tutorial.LastFunctionCalled())
tutorial.Function10("foo", 1.0)
print(tutorial.LastFunctionCalled())


rv_int = tutorial.overload1(10)
print(tutorial.LastFunctionCalled(), rv_int)
-- This should call overload (double type, int num)
-- but instead calls (int num, int offset)
-- since there is only one number type
rv_int = tutorial.overload1(1.0, 10)
print(tutorial.LastFunctionCalled(), rv_int)

rv_int = tutorial.overload1(10, 11, 12)
print(tutorial.LastFunctionCalled(), rv_int)
rv_int = tutorial.overload1(1.0, 10, 11, 12)
print(tutorial.LastFunctionCalled(), rv_int)

-- rv_int = tutorial.overload1("no such overload")


-- call a class
local obj = tutorial.Class1()
obj:Method1()
print(tutorial.LastFunctionCalled())

--XXX    tutorial.useclass(obj)
