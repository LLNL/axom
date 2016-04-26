// luaTutorialmodule.hpp
// This is generated code, do not edit
#ifndef LUATUTORIALMODULE_HPP
#define LUATUTORIALMODULE_HPP
#ifdef __cplusplus
extern "C" {
#endif
#include "tutorial.hpp"
#include "lua.h"

typedef struct {
    tutorial::Class1 * self;
} FFF;

int luaopen_tutorial(lua_State *L);

#ifdef __cplusplus
}
#endif
#endif  /* LUATUTORIALMODULE_HPP */
