
!-------------------------------------------------------------------------------
! Simple test program to determine if the c preprocessor flag defined in the
! build system variable ATK_PREPROCESS_FORTRAN is working
!
! Note: Our policy is for Fortran files that require the C preprocessort
!       to have the file extension *.F  (note uppercase F)
!       while files that do not require the preprocessor should have
!       the extension *.f  (note: lowercase f).
!       The current file has a lower case extension, but is compiled
!       with the ATK_PREPROCESS_FORTRAN flag to invoke the preprocessor.
!-------------------------------------------------------------------------------
program preprocessor_tester

#if 0
  Random stuff that should be removed if file is preprocessed
#endif

  print *, "Fortran's C preprocessor is working!"

end program
