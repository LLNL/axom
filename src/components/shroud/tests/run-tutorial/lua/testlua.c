#include <stdio.h>
#include <string.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

//#include <luaTutorialmodule.hpp>
int luaopen_tutorial(lua_State *L);

int main (void)
{
    lua_State *L = luaL_newstate();   /* opens Lua */
    luaL_openlibs(L);
    luaL_requiref(L, "tutorial", luaopen_tutorial, 1);    

#if 0
    char buff[256];
    int error;
    while (fgets(buff, sizeof(buff), stdin) != NULL) {
        error = luaL_loadbuffer(L, buff, strlen(buff), "line") ||
	    lua_pcall(L, 0, 0, 0);
        if (error) {
	    fprintf(stderr, "%s", lua_tostring(L, -1));
	    lua_pop(L, 1);  /* pop error message from the stack */
        }
    }
#else
    if (luaL_dofile(L, "test.lua")) {
	luaL_error(L, "error running script: %s", lua_tostring(L, -1));
    }
#endif
    
    lua_close(L);
    return 0;
}
