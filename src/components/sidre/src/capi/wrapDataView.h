// wrapDataView.h
// For C users and C++ implementation

#ifndef WRAPDATAVIEW_H
#define WRAPDATAVIEW_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef EXAMPLE_WRAPPER_IMPL
typedef void DS_dataview;
#else
struct s_DS_dataview;
typedef struct s_DS_dataview DS_dataview;
#endif

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATAVIEW_H
