!
! DataTypes.h
!
! Data enums used by C and Fortran

!
! Must be kept in sync with datatypes.h
!
enum, bind(C)
   enumerator :: INT8_T = 3
   enumerator :: INT16_T
   enumerator :: INT32_T
   enumerator :: INT64_T
   enumerator :: UINT8_T
   enumerator :: UINT16_T
   enumerator :: UINT32_T
   enumerator :: UINT64_T
   enumerator :: FLOAT32_T
   enumerator :: FLOAT64_T 
   enumerator :: CHAR8_STR_T
end enum
