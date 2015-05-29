//
// DataTypes.h
//
// Data enums used by C and Fortran

//
// Copied from src/components/conduit/src/conduit/DataType.h
//  Must be kept in sync with datatypes.h

typedef enum
    {
//        EMPTY_T = 0, // empty (default type)
//        OBJECT_T,    // object
//        LIST_T,      // list
        INT8_T = 3,  // int8
        INT16_T,     // int16
        INT32_T,     // int32
        INT64_T,     // int64
        UINT8_T,     // int8
        UINT16_T,    // uint16
        UINT32_T,    // uint32
        UINT64_T,    // uint64
        FLOAT32_T,   // float32
        FLOAT64_T,   // float64
        CHAR8_STR_T, // char8 string (incore c-string)
    } ATK_TypeID;
