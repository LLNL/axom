
!-------------------------------------------------------------------------------
! Simple test program to determine if the c preprocessor flag defined in the
! build system variable ATK_PREPROCESS_FORTRAN is working
!-------------------------------------------------------------------------------
program preprocessor_tester

#if 0
  Random stuff that should be removed if file is preprocessed
#endif

  print *, "Fortran's C preprocessor is working!"

end program
