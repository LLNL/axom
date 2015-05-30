//
// DataTypes.h
//
// Data enums used by C and Fortran
#ifndef SIDRE_DATATYPES_H
#define SIDRE_DATATYPES_H

//
// Copied from src/components/conduit/src/conduit/DataType.h
//  Must be kept in sync with datatypes.h

typedef enum
    {
//        EMPTY_T = 0, // empty (default type)
//        OBJECT_T,    // object
//        LIST_T,      // list
        ATK_INT8_T = 3,  // int8
        ATK_INT16_T,     // int16
        ATK_INT32_T,     // int32
        ATK_INT64_T,     // int64
        ATK_UINT8_T,     // int8
        ATK_UINT16_T,    // uint16
        ATK_UINT32_T,    // uint32
        ATK_UINT64_T,    // uint64
        ATK_FLOAT32_T,   // float32
        ATK_FLOAT64_T,   // float64
        ATK_CHAR8_STR_T, // char8 string (incore c-string)
    } ATK_TypeEnum;

#endif  // SIDRE_DATATYPES_H
