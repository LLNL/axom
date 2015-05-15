// wrapDataBuffer.h
// For C users and C++ implementation

#ifndef WRAPDATABUFFER_H
#define WRAPDATABUFFER_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef EXAMPLE_WRAPPER_IMPL
typedef void DS_databuffer;
#else
struct s_DS_databuffer;
typedef struct s_DS_databuffer DS_databuffer;
#endif

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATABUFFER_H
