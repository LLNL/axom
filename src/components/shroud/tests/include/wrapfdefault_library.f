! wrapfdefault_library.f
! This is generated code, do not edit
!>
!! \file wrapfdefault_library.f
!! \brief Shroud generated wrapper for default_library library
!<
module default_library_mod
    implicit none


    interface

        subroutine function1() &
                bind(C, name="DEF_function1")
            implicit none
        end subroutine function1

    end interface

contains


end module default_library_mod
