!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
module databuffer_mod
    use fstr_mod
    
    type databuffer
        type(C_PTR) obj
    contains
        procedure :: get_uid => databuffer_get_uid
    end type databuffer
    
    interface
        
        function atk_databuffer_get_uid(self) result(rv) bind(C, name="ATK_databuffer_get_uid")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT) :: rv
        end function atk_databuffer_get_uid
    end interface

contains
    
    function databuffer_get_uid(obj) result(rv)
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: rv
        ! splicer begin
        rv = atk_databuffer_get_uid(obj%obj)
        ! splicer end
    end function databuffer_get_uid

end module databuffer_mod
