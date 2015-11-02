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
module slic_mod
    use fstr_mod
    implicit none
    
    
    interface
        
        subroutine atk_initialize() &
                bind(C, name="ATK_initialize")
            use iso_c_binding
            implicit none
        end subroutine atk_initialize
        
        function atk_is_initialized() result(rv) &
                bind(C, name="ATK_is_initialized")
            use iso_c_binding
            implicit none
            logical(C_BOOL) :: rv
        end function atk_is_initialized
        
        subroutine atk_finalize() &
                bind(C, name="ATK_finalize")
            use iso_c_binding
            implicit none
        end subroutine atk_finalize
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    subroutine initialize()
        use iso_c_binding
        implicit none
        ! splicer begin initialize
        call atk_initialize()
        ! splicer end initialize
    end subroutine initialize
    
    function is_initialized() result(rv)
        use iso_c_binding
        implicit none
        logical :: rv
        ! splicer begin is_initialized
        rv = atk_is_initialized()
        ! splicer end is_initialized
    end function is_initialized
    
    subroutine finalize()
        use iso_c_binding
        implicit none
        ! splicer begin finalize
        call atk_finalize()
        ! splicer end finalize
    end subroutine finalize
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module slic_mod
