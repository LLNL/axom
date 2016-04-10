// luaTutorialpackage.cpp
// This is generated code, do not edit

static int l_function1(lua_State *L)
{
    Function1();
    return 0;
}

static int l_function2(lua_State *L)
{
    double arg1 = lua_tonumber(L, 0);
    int arg2 = lua_tointeger(L, 0);
    
    double rv = Function2(arg1, arg2);
    lua_pushnumber(L, rv)
    return 1;
}

static int l_sum(lua_State *L)
{
    int len = lua_tointeger(L, 0);
    int * values = lua_tointeger(L, 0);
    int * result;
    
    Sum(len, values, result);
    lua_pushinteger(L, result)
    return 1;
}

static int l_function3(lua_State *L)
{
    bool arg = lua_toboolean(L, 0);
    
    bool rv = Function3(arg);
    lua_pushboolean(L, rv)
    return 1;
}

static int l_function4a(lua_State *L)
{
    char * arg1 = lua_tostring(L, 0);
    char * arg2 = lua_tostring(L, 0);
    
    const std::string rv = Function4a(arg1, arg2);
    lua_pushstring(L, rv.c_str())
    return 1;
}

static int l_function4b(lua_State *L)
{
    char * arg1 = lua_tostring(L, 0);
    char * arg2 = lua_tostring(L, 0);
    
    const std::string & rv = Function4b(arg1, arg2);
    lua_pushstring(L, rv.c_str())
    return 1;
}

static int l_function5_arg1_arg2(lua_State *L)
{
    int shroud_nargs = lua_gettop(L);
    double arg1 = lua_tonumber(L, 0);
    bool arg2 = lua_toboolean(L, 0);
    double rv;
    
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
    lua_pushnumber(L, rv)
    return 1;
}

static int l_function6_from_name(lua_State *L)
{
    char * name = lua_tostring(L, 0);
    
    Function6(name);
    return 0;
}

static int l_function6_from_index(lua_State *L)
{
    int indx = lua_tointeger(L, 0);
    
    Function6(indx);
    return 0;
}

static int l_function9(lua_State *L)
{
    double arg = lua_tonumber(L, 0);
    
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
    char * name = lua_tostring(L, 0);
    double arg2 = lua_tonumber(L, 0);
    
    Function10(name, arg2);
    return 0;
}

static int l_overload1_num_offset_stride(lua_State *L)
{
    int shroud_nargs = lua_gettop(L);
    int num = lua_tointeger(L, 0);
    int offset = lua_tointeger(L, 0);
    int stride = lua_tointeger(L, 0);
    int rv;
    
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
    lua_pushinteger(L, rv)
    return 1;
}

static int l_overload1_5(lua_State *L)
{
    int shroud_nargs = lua_gettop(L);
    double type = lua_tonumber(L, 0);
    int num = lua_tointeger(L, 0);
    int offset = lua_tointeger(L, 0);
    int stride = lua_tointeger(L, 0);
    int rv;
    
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
    lua_pushinteger(L, rv)
    return 1;
}

static int l_typefunc(lua_State *L)
{
    TypeID arg = lua_tointeger(L, 0);
    
    TypeID rv = typefunc(arg);
    lua_pushinteger(L, rv)
    return 1;
}

static int l_enumfunc(lua_State *L)
{
    EnumTypeID arg = lua_tointeger(L, 0);
    
    EnumTypeID rv = enumfunc(static_cast<EnumTypeID>(arg));
    lua_pushinteger(L, static_cast<int>(rv))
    return 1;
}

static int l_last_function_called(lua_State *L)
{
    const std::string & rv = LastFunctionCalled();
    lua_pushstring(L, rv.c_str())
    return 1;
}
