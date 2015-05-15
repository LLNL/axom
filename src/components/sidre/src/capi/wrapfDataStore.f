module datastore_mod
    use fstr_mod
    use databuffer_mod, only : databuffer
    use datagroup_mod, only : datagroup
    
    type datastore
        type(C_PTR) obj
    contains
        procedure :: create_buffer => datastore_create_buffer
        procedure :: get_root => datastore_get_root
    end type datastore
    
    interface
        
        function ds_datastore_new() result(rv) bind(C, name="DS_datastore_new")
            use iso_c_binding
            implicit none
            type(C_PTR) :: rv
        end function ds_datastore_new
        
        subroutine ds_datastore_delete(self) bind(C, name="DS_datastore_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
        end subroutine ds_datastore_delete
        
        function ds_datastore_create_buffer(self) result(rv) bind(C, name="DS_datastore_create_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function ds_datastore_create_buffer
        
        function ds_datastore_get_root(self) result(rv) bind(C, name="DS_datastore_get_root")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function ds_datastore_get_root
    end interface

contains
    
    function datastore_new() result(rv)
        implicit none
        type(datastore) :: rv
        ! splicer begin
        rv%obj = ds_datastore_new()
        ! splicer end
    end function datastore_new
    
    subroutine datastore_delete(obj)
        implicit none
        type(datastore) :: obj
        ! splicer begin
        call ds_datastore_delete(obj%obj)
        obj%obj = C_NULL_PTR
        ! splicer end
    end subroutine datastore_delete
    
    function datastore_create_buffer(obj) result(rv)
        implicit none
        class(datastore) :: obj
        type(databuffer) :: rv
        ! splicer begin
        rv%obj = ds_datastore_create_buffer(obj%obj)
        ! splicer end
    end function datastore_create_buffer
    
    function datastore_get_root(obj) result(rv)
        implicit none
        class(datastore) :: obj
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = ds_datastore_get_root(obj%obj)
        ! splicer end
    end function datastore_get_root

end module datastore_mod
