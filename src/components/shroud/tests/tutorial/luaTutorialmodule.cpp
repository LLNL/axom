// luaTutorialmodule.cpp
// This is generated code, do not edit
#include "tutorial.hpp"
#include "luaTutorialmodule.hpp"
#ifdef __cplusplus
extern "C" {
#endif
#include "lauxlib.h"
#ifdef __cplusplus
}
#endif
// splicer begin include
// splicer end include

namespace tutorial {
// splicer begin C_definition
// splicer end C_definition

static int l_class1_new(lua_State *L)
{
    // splicer begin class.Class1.method.new
    l_Class1_Type * SH_this = (l_Class1_Type *) lua_newuserdata(L, sizeof(*SH_this));
    SH_this->self = new Class1();
    /* Add the metatable to the stack. */
    luaL_getmetatable(L, "Class1.metatable");
    /* Set the metatable on the userdata. */
    lua_setmetatable(L, -2);
    return 1;
    // splicer end class.Class1.method.new
}

static int l_class1_delete(lua_State *L)
{
    // splicer begin class.Class1.method.__gc
    l_Class1_Type * SH_this = (l_Class1_Type *)luaL_checkudata(L, 1, "Class1.metatable");
    delete SH_this->self;
    SH_this->self = NULL;
    return 0;
    // splicer end class.Class1.method.__gc
}

static int l_class1_method1(lua_State *L)
{
    // splicer begin class.Class1.method.Method1
    l_Class1_Type * SH_this = (l_Class1_Type *)luaL_checkudata(L, 1, "Class1.metatable");
    SH_this->self->Method1();
    return 0;
    // splicer end class.Class1.method.Method1
}

static const struct luaL_Reg l_Class1_Reg [] = {
    {"__gc", l_class1_delete},
    {"Method1", l_class1_method1},
    {NULL, NULL}   /*sentinel */
};

static int l_function1(lua_State *L)
{
    // splicer begin function.Function1
    Function1();
    return 0;
    // splicer end function.Function1
}

static int l_function2(lua_State *L)
{
    // splicer begin function.Function2
    double arg1 = lua_tonumber(L, 1);
    int arg2 = lua_tointeger(L, 2);
    double rv = Function2(arg1, arg2);
    lua_pushnumber(L, rv);
    return 1;
    // splicer end function.Function2
}

static int l_function3(lua_State *L)
{
    // splicer begin function.Function3
    bool arg = lua_toboolean(L, 1);
    bool rv = Function3(arg);
    lua_pushboolean(L, rv);
    return 1;
    // splicer end function.Function3
}

static int l_function4a(lua_State *L)
{
    // splicer begin function.Function4a
    const char * arg1 = lua_tostring(L, 1);
    const char * arg2 = lua_tostring(L, 2);
    const std::string rv = Function4a(arg1, arg2);
    lua_pushstring(L, rv.c_str());
    return 1;
    // splicer end function.Function4a
}

static int l_function4b(lua_State *L)
{
    // splicer begin function.Function4b
    const char * arg1 = lua_tostring(L, 1);
    const char * arg2 = lua_tostring(L, 2);
    const std::string & rv = Function4b(arg1, arg2);
    lua_pushstring(L, rv.c_str());
    return 1;
    // splicer end function.Function4b
}

static int l_function5(lua_State *L)
{
    // splicer begin function.Function5
    int SH_nresult;
    int SH_nargs = lua_gettop(L);
    int SH_itype1 = lua_type(L, 1);
    int SH_itype2 = lua_type(L, 2);
    switch (SH_nargs) {
    case 0:
        {
            double rv = Function5();
            lua_pushnumber(L, rv);
            SH_nresult = 1;
        }
        break;
    case 1:
        if (SH_itype1 == LUA_TNUMBER) {
            double arg1 = lua_tonumber(L, 1);
            double rv = Function5(arg1);
            lua_pushnumber(L, rv);
            SH_nresult = 1;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    case 2:
        if (SH_itype1 == LUA_TNUMBER &&
            SH_itype2 == LUA_TBOOLEAN) {
            double arg1 = lua_tonumber(L, 1);
            bool arg2 = lua_toboolean(L, 2);
            double rv = Function5(arg1, arg2);
            lua_pushnumber(L, rv);
            SH_nresult = 1;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    default:
        luaL_error(L, "error with arguments");
        break;
    }
    return SH_nresult;
    // splicer end function.Function5
}

static int l_function6(lua_State *L)
{
    // splicer begin function.Function6
    int SH_nresult;
    int SH_nargs = lua_gettop(L);
    int SH_itype1 = lua_type(L, 1);
    switch (SH_nargs) {
    case 1:
        if (SH_itype1 == LUA_TSTRING) {
            const char * name = lua_tostring(L, 1);
            Function6(name);
            SH_nresult = 0;
        }
        else if (SH_itype1 == LUA_TNUMBER) {
            int indx = lua_tointeger(L, 1);
            Function6(indx);
            SH_nresult = 0;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    default:
        luaL_error(L, "error with arguments");
        break;
    }
    return SH_nresult;
    // splicer end function.Function6
}

static int l_function9(lua_State *L)
{
    // splicer begin function.Function9
    double arg = lua_tonumber(L, 1);
    Function9(arg);
    return 0;
    // splicer end function.Function9
}

static int l_function10(lua_State *L)
{
    // splicer begin function.Function10
    int SH_nresult;
    int SH_nargs = lua_gettop(L);
    int SH_itype1 = lua_type(L, 1);
    int SH_itype2 = lua_type(L, 2);
    switch (SH_nargs) {
    case 0:
        {
            Function10();
            SH_nresult = 0;
        }
        break;
    case 2:
        if (SH_itype1 == LUA_TSTRING &&
            SH_itype2 == LUA_TNUMBER) {
            const char * name = lua_tostring(L, 1);
            double arg2 = lua_tonumber(L, 2);
            Function10(name, arg2);
            SH_nresult = 0;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    default:
        luaL_error(L, "error with arguments");
        break;
    }
    return SH_nresult;
    // splicer end function.Function10
}

static int l_overload1(lua_State *L)
{
    // splicer begin function.overload1
    int SH_nresult;
    int SH_nargs = lua_gettop(L);
    int SH_itype1 = lua_type(L, 1);
    int SH_itype2 = lua_type(L, 2);
    int SH_itype3 = lua_type(L, 3);
    int SH_itype4 = lua_type(L, 4);
    switch (SH_nargs) {
    case 1:
        if (SH_itype1 == LUA_TNUMBER) {
            int num = lua_tointeger(L, 1);
            int rv = overload1(num);
            lua_pushinteger(L, rv);
            SH_nresult = 1;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    case 2:
        if (SH_itype1 == LUA_TNUMBER &&
            SH_itype2 == LUA_TNUMBER) {
            int num = lua_tointeger(L, 1);
            int offset = lua_tointeger(L, 2);
            int rv = overload1(num, offset);
            lua_pushinteger(L, rv);
            SH_nresult = 1;
        }
        else if (SH_itype1 == LUA_TNUMBER &&
            SH_itype2 == LUA_TNUMBER) {
            double type = lua_tonumber(L, 1);
            int num = lua_tointeger(L, 2);
            int rv = overload1(type, num);
            lua_pushinteger(L, rv);
            SH_nresult = 1;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    case 3:
        if (SH_itype1 == LUA_TNUMBER &&
            SH_itype2 == LUA_TNUMBER &&
            SH_itype3 == LUA_TNUMBER) {
            int num = lua_tointeger(L, 1);
            int offset = lua_tointeger(L, 2);
            int stride = lua_tointeger(L, 3);
            int rv = overload1(num, offset, stride);
            lua_pushinteger(L, rv);
            SH_nresult = 1;
        }
        else if (SH_itype1 == LUA_TNUMBER &&
            SH_itype2 == LUA_TNUMBER &&
            SH_itype3 == LUA_TNUMBER) {
            double type = lua_tonumber(L, 1);
            int num = lua_tointeger(L, 2);
            int offset = lua_tointeger(L, 3);
            int rv = overload1(type, num, offset);
            lua_pushinteger(L, rv);
            SH_nresult = 1;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    case 4:
        if (SH_itype1 == LUA_TNUMBER &&
            SH_itype2 == LUA_TNUMBER &&
            SH_itype3 == LUA_TNUMBER &&
            SH_itype4 == LUA_TNUMBER) {
            double type = lua_tonumber(L, 1);
            int num = lua_tointeger(L, 2);
            int offset = lua_tointeger(L, 3);
            int stride = lua_tointeger(L, 4);
            int rv = overload1(type, num, offset, stride);
            lua_pushinteger(L, rv);
            SH_nresult = 1;
        }
        else {
            luaL_error(L, "error with arguments");
        }
        break;
    default:
        luaL_error(L, "error with arguments");
        break;
    }
    return SH_nresult;
    // splicer end function.overload1
}

static int l_typefunc(lua_State *L)
{
    // splicer begin function.typefunc
    TypeID arg = lua_tointeger(L, 1);
    TypeID rv = typefunc(arg);
    lua_pushinteger(L, rv);
    return 1;
    // splicer end function.typefunc
}

static int l_enumfunc(lua_State *L)
{
    // splicer begin function.enumfunc
    EnumTypeID arg = static_cast<EnumTypeID>(lua_tointeger(L, 1));
    EnumTypeID rv = enumfunc(arg);
    lua_pushinteger(L, static_cast<int>(rv));
    return 1;
    // splicer end function.enumfunc
}

static int l_last_function_called(lua_State *L)
{
    // splicer begin function.LastFunctionCalled
    const std::string & rv = LastFunctionCalled();
    lua_pushstring(L, rv.c_str());
    return 1;
    // splicer end function.LastFunctionCalled
}

static const struct luaL_Reg XXX1 [] = {
    {"Class1", l_class1_new},
    {"Function1", l_function1},
    {"Function2", l_function2},
    {"Function3", l_function3},
    {"Function4a", l_function4a},
    {"Function4b", l_function4b},
    {"Function5", l_function5},
    {"Function6", l_function6},
    {"Function9", l_function9},
    {"Function10", l_function10},
    {"overload1", l_overload1},
    {"typefunc", l_typefunc},
    {"enumfunc", l_enumfunc},
    {"LastFunctionCalled", l_last_function_called},
    {NULL, NULL}   /*sentinel */
};

#ifdef __cplusplus
extern "C" {
#endif
int luaopen_tutorial(lua_State *L) {
    
    /* Create the metatable and put it on the stack. */
    luaL_newmetatable(L, "Class1.metatable");
    /* Duplicate the metatable on the stack (We now have 2). */
    lua_pushvalue(L, -1);
    /* Pop the first metatable off the stack and assign it to __index
     * of the second one. We set the metatable for the table to itself.
     * This is equivalent to the following in lua:
     * metatable = {}
     * metatable.__index = metatable
     */
    lua_setfield(L, -2, "__index");
     
    /* Set the methods to the metatable that should be accessed via object:func */
    luaL_setfuncs(L, l_Class1_Reg, 0);
    
    luaL_newlib(L, XXX1);
    return 1;
}
#ifdef __cplusplus
}
#endif
}  // namespace tutorial
