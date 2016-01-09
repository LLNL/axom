! wrapfslic.f
! This is generated code, do not edit
!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
!>
!! \file wrapfslic.f
!! \brief Shroud generated wrapper for SLIC library
!<
module slic_mod
    use fstr_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none
    
    
    interface
        
        subroutine initialize() &
                bind(C, name="SLIC_initialize")
            use iso_c_binding
            implicit none
        end subroutine initialize
        
        function is_initialized() &
                result(rv) &
                bind(C, name="SLIC_is_initialized")
            use iso_c_binding
            implicit none
            logical(C_BOOL) :: rv
        end function is_initialized
        
        subroutine finalize() &
                bind(C, name="SLIC_finalize")
            use iso_c_binding
            implicit none
        end subroutine finalize
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module slic_mod
