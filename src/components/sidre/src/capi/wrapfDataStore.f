!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
module datastore_mod
    use fstr_mod
    use databuffer_mod, only : databuffer
    use datagroup_mod, only : datagroup
    
    type datastore
        type(C_PTR) obj
    contains
        procedure :: create_buffer => datastore_create_buffer
        procedure :: destroy_buffer => datastore_destroy_buffer
        procedure :: get_root => datastore_get_root
    end type datastore
    
    interface
        
        function atk_datastore_new() result(rv) bind(C, name="ATK_datastore_new")
            use iso_c_binding
            implicit none
            type(C_PTR) :: rv
        end function atk_datastore_new
        
        subroutine atk_datastore_delete(self) bind(C, name="ATK_datastore_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
        end subroutine atk_datastore_delete
        
        function atk_datastore_create_buffer(self) result(rv) bind(C, name="ATK_datastore_create_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_datastore_create_buffer
        
        subroutine atk_datastore_destroy_buffer(id) bind(C, name="ATK_datastore_destroy_buffer")
            use iso_c_binding
            implicit none
            integer(C_INT), value :: id
        end subroutine atk_datastore_destroy_buffer
        
        function atk_datastore_get_root(self) result(rv) bind(C, name="ATK_datastore_get_root")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_datastore_get_root
    end interface

contains
    
    function datastore_new() result(rv)
        implicit none
        type(datastore) :: rv
        ! splicer begin
        rv%obj = atk_datastore_new()
        ! splicer end
    end function datastore_new
    
    subroutine datastore_delete(obj)
        implicit none
        type(datastore) :: obj
        ! splicer begin
        call atk_datastore_delete(obj%obj)
        obj%obj = C_NULL_PTR
        ! splicer end
    end subroutine datastore_delete
    
    function datastore_create_buffer(obj) result(rv)
        implicit none
        class(datastore) :: obj
        type(databuffer) :: rv
        ! splicer begin
        rv%obj = atk_datastore_create_buffer(obj%obj)
        ! splicer end
    end function datastore_create_buffer
    
    subroutine datastore_destroy_buffer(id)
        implicit none
        integer(C_INT) :: id
        ! splicer begin
        call atk_datastore_destroy_buffer(id)
        ! splicer end
    end subroutine datastore_destroy_buffer
    
    function datastore_get_root(obj) result(rv)
        implicit none
        class(datastore) :: obj
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = atk_datastore_get_root(obj%obj)
        ! splicer end
    end function datastore_get_root

end module datastore_mod
