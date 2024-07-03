module sina_functions

  interface
    
    subroutine create_document_and_record(id) 
      character(*) id  
    end subroutine create_document_and_record
    
  end interface
  
  interface
    
    subroutine sina_add_file(file_nm, mime_type) 
      character(*) file_nm  
      character(*) mime_type
    end subroutine sina_add_file
    
  end interface
  
  interface
    
    subroutine write_sina_document(file_nm) 
      character(*) file_nm  
    end subroutine write_sina_document
    
  end interface

  interface sina_add
  
    subroutine sina_add_long(key, value, units, tags)
      character(*) key
      integer (KIND=8) value
      character(*) units
      character(*) tags
    end subroutine sina_add_long
    
    subroutine sina_add_logical(key, value, units, tags)
      character(*) key
      logical value
      character(*) units
      character(*) tags
    end subroutine sina_add_logical
    
    subroutine sina_add_int(key, value, units, tags)
      character(*) key
      integer value
      character(*) units
      character(*) tags
    end subroutine sina_add_int
    
    subroutine sina_add_double(key, value, units, tags)
      character(*) key
      double precision value
      character(*) units
      character(*) tags
    end subroutine sina_add_double

    subroutine sina_add_float(key, value, units, tags)
      character(*) key
      real value
      character(*) units
      character(*) tags
    end subroutine sina_add_float

    subroutine sina_add_string(key, value, units, tags)
      character(*) key
      character(*) value
      character(*) units
      character(*) tags
    end subroutine sina_add_string
    
  end interface

  interface sina_add_curveset
  
    subroutine sina_add_curveset(name)
      character(*) name
    end subroutine sina_add_curveset
    
  end interface
  
  interface sina_add_curve
  
    subroutine sina_add_curve_double(name, curve, values, n, independent)
      character(*) name
      character(*) curve
      double precision values(n)
      integer n
      logical independent
    end subroutine sina_add_curve_double
    
    subroutine sina_add_curve_float(name, curve, values, n, independent)
      character(*) name
      character(*) curve
      real values(n)
      integer n
      logical independent
    end subroutine sina_add_curve_float
    
    subroutine sina_add_curve_int(name, curve, values, n, independent)
      character(*) name
      character(*) curve
      integer (KIND=4), dimension(n) :: values
      integer n
      logical independent
    end subroutine sina_add_curve_int
    
    subroutine sina_add_curve_long(name, curve, values, n, independent)
      character(*) name
      character(*) curve
      integer (KIND=8), dimension(n) :: values
      integer n
      logical independent
    end subroutine sina_add_curve_long
    
  end interface
  
end module 