"""
Hide some literal Fortran code in this module
"""
from __future__ import print_function

import os

fstr_module = """
module fstr_mod

  interface fstr
     module procedure fstr_ptr, fstr_arr
  end interface

  interface
     pure function strlen_ptr(s) result(result) bind(c,name="strlen")
       use, intrinsic :: iso_c_binding
       integer(c_int) :: result
       type(c_ptr), value, intent(in) :: s
     end function strlen_ptr
  end interface

contains

  ! Count the characters in a null-terminated array.
  pure function strlen_arr(s)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char
    character(kind=c_char, len=1), intent(in) :: s(*)
    integer :: i, strlen_arr
    i=1
    do
       if (s(i) == c_null_char) exit
       i = i+1
    enddo
    strlen_arr = i-1
  end function strlen_arr

  ! Convert a null-terminated C "char *" pointer to a Fortran string.
  function fstr_ptr(s) result(fs)
    use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_f_pointer
    type(c_ptr), intent(in) :: s
    character(kind=c_char, len=strlen_ptr(s)) :: fs
    character(kind=c_char), pointer :: cptr(:)
    integer :: i
    call c_f_pointer(s, cptr, [len(fs)])
    do i=1, len(fs)
       fs(i:i) = cptr(i)
    enddo
  end function fstr_ptr

  ! Convert a null-terminated array of characters to a Fortran string.
  function fstr_arr(s) result(fs)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char
    character(kind=c_char, len=1), intent(in) :: s(*)
    character(kind=c_char, len=strlen_arr(s)) :: fs
    integer :: i
    do i = 1, len(fs)
       fs(i:i) = s(i)
    enddo
  end function fstr_arr

  ! Convert a fortran string in 's' to a null-terminated array of characters.
  pure function cstr(s)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char
    character(len=*), intent(in) :: s
    character(kind=c_char, len=1) :: cstr(len_trim(s)+1)
    integer :: i
    if (len_trim(s) > 0) cstr = [ (s(i:i), i=1,len_trim(s)) ]
    cstr(len_trim(s)+1) = c_null_char
  end function cstr

end module fstr_mod
"""

def write_runtime(tree, config):
    fname = 'fstr_mod.f'
    fp = open(os.path.join(config.binary_dir, fname), 'w')

    for line in tree.get('copyright', []):
        if line:
            fp.write('! ' + line + '\n')
        else:
            fp.write('!\n')

    fp.write(fstr_module)
    fp.close()
    print("Wrote", fname)
