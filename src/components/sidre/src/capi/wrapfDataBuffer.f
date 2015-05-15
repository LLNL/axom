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
    end type databuffer
    
    interface
    end interface

contains

end module databuffer_mod
