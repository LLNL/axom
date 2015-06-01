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
        procedure :: get_index => databuffer_get_index
        procedure :: declare => databuffer_declare
        procedure :: allocate => databuffer_allocate
        procedure :: get_data => databuffer_get_data
    end type databuffer
    
    interface
        
        function atk_databuffer_get_index(self) result(rv) bind(C, name="ATK_databuffer_get_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT) :: rv
        end function atk_databuffer_get_index
        
        function atk_databuffer_declare(self, type, len) result(rv) bind(C, name="ATK_databuffer_declare")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT), value :: type
            integer(C_LONG), value :: len
            type(C_PTR) :: rv
        end function atk_databuffer_declare
        
        function atk_databuffer_allocate(self) result(rv) bind(C, name="ATK_databuffer_allocate")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_databuffer_allocate
        
        function atk_databuffer_get_data(self) result(rv) bind(C, name="ATK_databuffer_get_data")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_databuffer_get_data
    end interface

contains
    
    function databuffer_get_index(obj) result(rv)
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: rv
        ! splicer begin
        rv = atk_databuffer_get_index(obj%obj)
        ! splicer end
    end function databuffer_get_index
    
    function databuffer_declare(obj, type, len) result(rv)
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(databuffer) :: rv
        ! splicer begin
        rv%obj = atk_databuffer_declare(obj%obj, type, len)
        ! splicer end
    end function databuffer_declare
    
    function databuffer_allocate(obj) result(rv)
        implicit none
        class(databuffer) :: obj
        type(databuffer) :: rv
        ! splicer begin
        rv%obj = atk_databuffer_allocate(obj%obj)
        ! splicer end
    end function databuffer_allocate
    
    function databuffer_get_data(obj) result(rv)
        implicit none
        class(databuffer) :: obj
        type(C_PTR) :: rv
        ! splicer begin
        rv = atk_databuffer_get_data(obj%obj)
        ! splicer end
    end function databuffer_get_data

end module databuffer_mod
