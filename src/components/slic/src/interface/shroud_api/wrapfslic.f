! wrapfslic.f
! This is generated code, do not edit
!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
!>
!! \file wrapfslic.f
!! \brief Shroud generated wrapper for SLIC library
!<
module slic_mod
    use fstr_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none
    
    
    interface
        
        subroutine initialize() &
                bind(C, name="SLIC_initialize")
            use iso_c_binding
            implicit none
        end subroutine initialize
        
        function slic_is_initialized() &
                result(rv) &
                bind(C, name="SLIC_is_initialized")
            use iso_c_binding
            implicit none
            logical(C_BOOL) :: rv
        end function slic_is_initialized
        
        subroutine finalize() &
                bind(C, name="SLIC_finalize")
            use iso_c_binding
            implicit none
        end subroutine finalize
        
        subroutine slic_set_abort_on_assert(willAbort) &
                bind(C, name="SLIC_set_abort_on_assert")
            use iso_c_binding
            implicit none
            logical(C_BOOL), value, intent(IN) :: willAbort
        end subroutine slic_set_abort_on_assert
        
        function slic_get_abort_on_assert() &
                result(rv) &
                bind(C, name="SLIC_get_abort_on_assert")
            use iso_c_binding
            implicit none
            logical(C_BOOL) :: rv
        end function slic_get_abort_on_assert
        
        subroutine slic_set_abort_on_error(willAbort) &
                bind(C, name="SLIC_set_abort_on_error")
            use iso_c_binding
            implicit none
            logical(C_BOOL), value, intent(IN) :: willAbort
        end subroutine slic_set_abort_on_error
        
        function slic_get_abort_on_error() &
                result(rv) &
                bind(C, name="SLIC_get_abort_on_error")
            use iso_c_binding
            implicit none
            logical(C_BOOL) :: rv
        end function slic_get_abort_on_error
        
        subroutine slic_activate_logger(name) &
                bind(C, name="SLIC_activate_logger")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine slic_activate_logger
        
        subroutine slic_activate_logger_bufferify(name, Lname) &
                bind(C, name="SLIC_activate_logger_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine slic_activate_logger_bufferify
        
        subroutine set_logging_msg_level(level) &
                bind(C, name="SLIC_set_logging_msg_level")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: level
        end subroutine set_logging_msg_level
        
        subroutine slic_log_message(level, message, fileName, line, filter) &
                bind(C, name="SLIC_log_message")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: level
            character(kind=C_CHAR), intent(IN) :: message(*)
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: line
            logical(C_BOOL), value, intent(IN) :: filter
        end subroutine slic_log_message
        
        subroutine slic_log_message_bufferify(level, message, Lmessage, fileName, LfileName, line, filter) &
                bind(C, name="SLIC_log_message_bufferify")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: level
            character(kind=C_CHAR), intent(IN) :: message(*)
            integer(C_INT), value, intent(IN) :: Lmessage
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: LfileName
            integer(C_INT), value, intent(IN) :: line
            logical(C_BOOL), value, intent(IN) :: filter
        end subroutine slic_log_message_bufferify
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    function is_initialized() result(rv)
        use iso_c_binding
        implicit none
        logical :: rv
        ! splicer begin is_initialized
        rv = slic_is_initialized()
        ! splicer end is_initialized
    end function is_initialized
    
    subroutine set_abort_on_assert(willAbort)
        use iso_c_binding
        implicit none
        logical, value, intent(IN) :: willAbort
        logical(C_BOOL) tmp_willAbort
        tmp_willAbort = willAbort  ! coerce to C_BOOL
        ! splicer begin set_abort_on_assert
        call slic_set_abort_on_assert(tmp_willAbort)
        ! splicer end set_abort_on_assert
    end subroutine set_abort_on_assert
    
    function get_abort_on_assert() result(rv)
        use iso_c_binding
        implicit none
        logical :: rv
        ! splicer begin get_abort_on_assert
        rv = slic_get_abort_on_assert()
        ! splicer end get_abort_on_assert
    end function get_abort_on_assert
    
    subroutine set_abort_on_error(willAbort)
        use iso_c_binding
        implicit none
        logical, value, intent(IN) :: willAbort
        logical(C_BOOL) tmp_willAbort
        tmp_willAbort = willAbort  ! coerce to C_BOOL
        ! splicer begin set_abort_on_error
        call slic_set_abort_on_error(tmp_willAbort)
        ! splicer end set_abort_on_error
    end subroutine set_abort_on_error
    
    function get_abort_on_error() result(rv)
        use iso_c_binding
        implicit none
        logical :: rv
        ! splicer begin get_abort_on_error
        rv = slic_get_abort_on_error()
        ! splicer end get_abort_on_error
    end function get_abort_on_error
    
    subroutine activate_logger(name)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: name
        ! splicer begin activate_logger
        call slic_activate_logger_bufferify(  &
            name,  &
            len_trim(name))
        ! splicer end activate_logger
    end subroutine activate_logger
    
    subroutine log_message(level, message, fileName, line, filter)
        use iso_c_binding
        implicit none
        integer(C_INT), value, intent(IN) :: level
        character(*), intent(IN) :: message
        character(*), intent(IN) :: fileName
        integer(C_INT), value, intent(IN) :: line
        logical, value, intent(IN) :: filter
        logical(C_BOOL) tmp_filter
        tmp_filter = filter  ! coerce to C_BOOL
        ! splicer begin log_message
        call slic_log_message_bufferify(  &
            level,  &
            message,  &
            len_trim(message),  &
            fileName,  &
            len_trim(fileName),  &
            line,  &
            tmp_filter)
        ! splicer end log_message
    end subroutine log_message
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module slic_mod
