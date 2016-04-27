-- test tutorial module

local tutorial = require "tutorial"
local rv_int, rv_double

tutorial.Function1()
print(tutorial.LastFunctionCalled())

rv_double = tutorial.Function2(1.0, 4)
print(tutorial.LastFunctionCalled(), rv_double)

rv_int = tutorial.overload1(10)
print(tutorial.LastFunctionCalled(), rv_int)
rv_int = tutorial.overload1(1.0, 10)
print(tutorial.LastFunctionCalled(), rv_int)

rv_int = tutorial.overload1(10, 11, 12)
print(tutorial.LastFunctionCalled(), rv_int)
rv_int = tutorial.overload1(1.0, 10, 11, 12)
print(tutorial.LastFunctionCalled(), rv_int)


-- call a class
local obj = tutorial.Class1()
obj:Method1()
print(tutorial.LastFunctionCalled())

