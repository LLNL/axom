// luaTutorialpackage.cpp
// This is generated code, do not edit
#include "tutorial.hpp"
#include "luaTutorialpackage.hpp"
#include "lauxlib.h"
// splicer begin include
// splicer end include

namespace tutorial {
// splicer begin C_definition
// splicer end C_definition

static int l_function1(lua_State *L)
{
    Function1();
    return 0;
}

static int l_function2(lua_State *L)
{
    double arg1 = lua_tonumber(L, 1);
    int arg2 = lua_tointeger(L, 2);
    
    double rv = Function2(arg1, arg2);
    lua_pushnumber(L, rv);
    return 1;
}

static int l_function3(lua_State *L)
{
    bool arg = lua_toboolean(L, 1);
    
    bool rv = Function3(arg);
    lua_pushboolean(L, rv);
    return 1;
}

static int l_function4a(lua_State *L)
{
    const char * arg1 = lua_tostring(L, 1);
    const char * arg2 = lua_tostring(L, 2);
    
    const std::string rv = Function4a(arg1, arg2);
    lua_pushstring(L, rv.c_str());
    return 1;
}

static int l_function4b(lua_State *L)
{
    const char * arg1 = lua_tostring(L, 1);
    const char * arg2 = lua_tostring(L, 2);
    
    const std::string & rv = Function4b(arg1, arg2);
    lua_pushstring(L, rv.c_str());
    return 1;
}

static int l_function5_arg1_arg2(lua_State *L)
{
    int shroud_nargs = lua_gettop(L);
    double arg1;
    bool arg2;
    double rv;
    
    if (shroud_nargs > 0) {
        arg1 = lua_tonumber(L, 1);
    }
    if (shroud_nargs > 1) {
        arg2 = lua_toboolean(L, 2);
    }
    switch (shroud_nargs) {
    case 0:
        rv = Function5();
        break;
    case 1:
        rv = Function5(arg1);
        break;
    case 2:
        rv = Function5(arg1, arg2);
        break;
    }
    lua_pushnumber(L, rv);
    return 1;
}

static int l_function6_from_name(lua_State *L)
{
    const char * name = lua_tostring(L, 1);
    
    Function6(name);
    return 0;
}

static int l_function6_from_index(lua_State *L)
{
    int indx = lua_tointeger(L, 1);
    
    Function6(indx);
    return 0;
}

static int l_function9(lua_State *L)
{
    double arg = lua_tonumber(L, 1);
    
    Function9(arg);
    return 0;
}

static int l_function10_0(lua_State *L)
{
    Function10();
    return 0;
}

static int l_function10_1(lua_State *L)
{
    const char * name = lua_tostring(L, 1);
    double arg2 = lua_tonumber(L, 2);
    
    Function10(name, arg2);
    return 0;
}

static int l_overload1_num_offset_stride(lua_State *L)
{
    int shroud_nargs = lua_gettop(L);
    int num = lua_tointeger(L, 1);
    int offset;
    int stride;
    int rv;
    
    if (shroud_nargs > 1) {
        offset = lua_tointeger(L, 2);
    }
    if (shroud_nargs > 2) {
        stride = lua_tointeger(L, 3);
    }
    switch (shroud_nargs) {
    case 1:
        rv = overload1(num);
        break;
    case 2:
        rv = overload1(num, offset);
        break;
    case 3:
        rv = overload1(num, offset, stride);
        break;
    }
    lua_pushinteger(L, rv);
    return 1;
}

static int l_overload1_5(lua_State *L)
{
    int shroud_nargs = lua_gettop(L);
    double type = lua_tonumber(L, 1);
    int num = lua_tointeger(L, 2);
    int offset;
    int stride;
    int rv;
    
    if (shroud_nargs > 2) {
        offset = lua_tointeger(L, 3);
    }
    if (shroud_nargs > 3) {
        stride = lua_tointeger(L, 4);
    }
    switch (shroud_nargs) {
    case 2:
        rv = overload1(type, num);
        break;
    case 3:
        rv = overload1(type, num, offset);
        break;
    case 4:
        rv = overload1(type, num, offset, stride);
        break;
    }
    lua_pushinteger(L, rv);
    return 1;
}

static int l_typefunc(lua_State *L)
{
    TypeID arg = lua_tointeger(L, 1);
    
    TypeID rv = typefunc(arg);
    lua_pushinteger(L, rv);
    return 1;
}

static int l_enumfunc(lua_State *L)
{
    EnumTypeID arg = static_cast<EnumTypeID>(lua_tointeger(L, 1));
    
    EnumTypeID rv = enumfunc(arg);
    lua_pushinteger(L, static_cast<int>(rv));
    return 1;
}

static int l_last_function_called(lua_State *L)
{
    const std::string & rv = LastFunctionCalled();
    lua_pushstring(L, rv.c_str());
    return 1;
}
static const struct luaL_Reg XXX1 [] = {
    {"Function1", l_function1},
    {"Function2", l_function2},
    {"Function3", l_function3},
    {"Function4a", l_function4a},
    {"Function4b", l_function4b},
    {"Function5_arg1_arg2", l_function5_arg1_arg2},
    {"Function6_from_name", l_function6_from_name},
    {"Function6_from_index", l_function6_from_index},
    {"Function9", l_function9},
    {"Function10_0", l_function10_0},
    {"Function10_1", l_function10_1},
    {"overload1_num_offset_stride", l_overload1_num_offset_stride},
    {"overload1_5", l_overload1_5},
    {"typefunc", l_typefunc},
    {"enumfunc", l_enumfunc},
    {"LastFunctionCalled", l_last_function_called},
    {NULL, NULL}   /*sentinel */
};

int luaopen_tutorial (lua_State *L) {
    luaL_newlib(L, XXX1);
    return 1;
}
}  // namespace tutorial
