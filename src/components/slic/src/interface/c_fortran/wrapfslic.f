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
! splicer begin file_top
! splicer end file_top
module slic_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none


    interface

        subroutine slic_initialize() &
                bind(C, name="SLIC_initialize")
            use iso_c_binding
            implicit none
        end subroutine slic_initialize

        function c_is_initialized() &
                result(rv) &
                bind(C, name="SLIC_is_initialized")
            use iso_c_binding
            implicit none
            logical(C_BOOL) :: rv
        end function c_is_initialized

        subroutine slic_finalize() &
                bind(C, name="SLIC_finalize")
            use iso_c_binding
            implicit none
        end subroutine slic_finalize

        subroutine c_create_logger(name, imask) &
                bind(C, name="SLIC_create_logger")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            character(kind=C_CHAR), value, intent(IN) :: imask
        end subroutine c_create_logger

        subroutine c_create_logger_bufferify(name, Lname, imask) &
                bind(C, name="SLIC_create_logger_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            character(kind=C_CHAR), value, intent(IN) :: imask
        end subroutine c_create_logger_bufferify

        function c_activate_logger(name) &
                result(rv) &
                bind(C, name="SLIC_activate_logger")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: rv
        end function c_activate_logger

        function c_activate_logger_bufferify(name, Lname) &
                result(rv) &
                bind(C, name="SLIC_activate_logger_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: rv
        end function c_activate_logger_bufferify

        subroutine c_get_active_logger_name_bufferify(name, Lname) &
                bind(C, name="SLIC_get_active_logger_name_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine c_get_active_logger_name_bufferify

        subroutine slic_set_logging_msg_level(level) &
                bind(C, name="SLIC_set_logging_msg_level")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: level
        end subroutine slic_set_logging_msg_level

        subroutine c_log_message(level, message, fileName, line, filter) &
                bind(C, name="SLIC_log_message")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: level
            character(kind=C_CHAR), intent(IN) :: message(*)
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: line
            logical(C_BOOL), value, intent(IN) :: filter
        end subroutine c_log_message

        subroutine c_log_message_bufferify(level, message, Lmessage, fileName, LfileName, line, filter) &
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
        end subroutine c_log_message_bufferify

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains

    function slic_is_initialized() result(rv)
        use iso_c_binding, only : C_BOOL
        implicit none
        logical :: rv
        ! splicer begin is_initialized
        rv = c_is_initialized()
        ! splicer end is_initialized
    end function slic_is_initialized

    subroutine slic_create_logger(name, imask)
        use iso_c_binding, only : C_INT
        implicit none
        character(*), intent(IN) :: name
        character, value, intent(IN) :: imask
        ! splicer begin create_logger
        call c_create_logger_bufferify(  &
            name,  &
            len_trim(name, kind=C_INT),  &
            imask)
        ! splicer end create_logger
    end subroutine slic_create_logger

    function slic_activate_logger(name) result(rv)
        use iso_c_binding, only : C_BOOL, C_INT
        implicit none
        character(*), intent(IN) :: name
        logical :: rv
        ! splicer begin activate_logger
        rv = c_activate_logger_bufferify(  &
            name,  &
            len_trim(name, kind=C_INT))
        ! splicer end activate_logger
    end function slic_activate_logger

    subroutine slic_get_active_logger_name(name)
        use iso_c_binding, only : C_INT
        implicit none
        character(*), intent(OUT) :: name
        ! splicer begin get_active_logger_name
        call c_get_active_logger_name_bufferify(  &
            name,  &
            len(name, kind=C_INT))
        ! splicer end get_active_logger_name
    end subroutine slic_get_active_logger_name

    subroutine slic_log_message(level, message, fileName, line, filter)
        use iso_c_binding, only : C_BOOL, C_INT
        implicit none
        integer(C_INT), value, intent(IN) :: level
        character(*), intent(IN) :: message
        character(*), intent(IN) :: fileName
        integer(C_INT), value, intent(IN) :: line
        logical, value, intent(IN) :: filter
        logical(C_BOOL) tmp_filter
        tmp_filter = filter  ! coerce to C_BOOL
        ! splicer begin log_message
        call c_log_message_bufferify(  &
            level,  &
            message,  &
            len_trim(message, kind=C_INT),  &
            fileName,  &
            len_trim(fileName, kind=C_INT),  &
            line,  &
            tmp_filter)
        ! splicer end log_message
    end subroutine slic_log_message

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module slic_mod
