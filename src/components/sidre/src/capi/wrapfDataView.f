!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
module dataview_mod
    use fstr_mod
    use datagroup_mod, only : datagroup
    
    type dataview
        type(C_PTR) obj
    contains
        procedure :: get_owning_group => dataview_get_owning_group
    end type dataview
    
    interface
        
        function atk_dataview_get_owning_group(self) result(rv) bind(C, name="ATK_dataview_get_owning_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_owning_group
    end interface

contains
    
    function dataview_get_owning_group(obj) result(rv)
        implicit none
        class(dataview) :: obj
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = atk_dataview_get_owning_group(obj%obj)
        ! splicer end
    end function dataview_get_owning_group

end module dataview_mod
