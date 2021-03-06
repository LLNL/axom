! wrapfslic.f
! This file is generated by Shroud 0.12.2. Do not edit.
!
! Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)
!>
!! \file wrapfslic.f
!! \brief Shroud generated wrapper for slic namespace
!<
! splicer begin file_top
! splicer end file_top
module axom_slic
    use iso_c_binding, only : C_INT
    ! splicer begin module_use
    ! splicer end module_use
    implicit none

    ! splicer begin module_top
    ! splicer end module_top

    !  enum axom::slic::message::Level
    integer(C_INT), parameter :: message_error = 0
    integer(C_INT), parameter :: message_warning = 1
    integer(C_INT), parameter :: message_info = 2
    integer(C_INT), parameter :: message_debug = 3
    integer(C_INT), parameter :: message_num_levels = 4

    interface

        subroutine slic_initialize() &
                bind(C, name="SLIC_initialize")
            implicit none
        end subroutine slic_initialize

        function c_is_initialized() &
                result(SHT_rv) &
                bind(C, name="SLIC_is_initialized")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL) :: SHT_rv
        end function c_is_initialized

        subroutine c_create_logger(name, imask) &
                bind(C, name="SLIC_create_logger")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            character(kind=C_CHAR), value, intent(IN) :: imask
        end subroutine c_create_logger

        subroutine c_create_logger_bufferify(name, Lname, imask) &
                bind(C, name="SLIC_create_logger_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            character(kind=C_CHAR), value, intent(IN) :: imask
        end subroutine c_create_logger_bufferify

        function c_activate_logger(name) &
                result(SHT_rv) &
                bind(C, name="SLIC_activate_logger")
            use iso_c_binding, only : C_BOOL, C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: SHT_rv
        end function c_activate_logger

        function c_activate_logger_bufferify(name, Lname) &
                result(SHT_rv) &
                bind(C, name="SLIC_activate_logger_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: SHT_rv
        end function c_activate_logger_bufferify

        subroutine c_get_active_logger_name_bufferify(name, Nname) &
                bind(C, name="SLIC_get_active_logger_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: name(*)
            integer(C_INT), value, intent(IN) :: Nname
        end subroutine c_get_active_logger_name_bufferify

        function slic_get_logging_msg_level() &
                result(SHT_rv) &
                bind(C, name="SLIC_get_logging_msg_level")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT) :: SHT_rv
        end function slic_get_logging_msg_level

        subroutine slic_set_logging_msg_level(level) &
                bind(C, name="SLIC_set_logging_msg_level")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: level
        end subroutine slic_set_logging_msg_level

        subroutine c_set_abort_on_error(status) &
                bind(C, name="SLIC_set_abort_on_error")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL), value, intent(IN) :: status
        end subroutine c_set_abort_on_error

        subroutine slic_enable_abort_on_error() &
                bind(C, name="SLIC_enable_abort_on_error")
            implicit none
        end subroutine slic_enable_abort_on_error

        subroutine slic_disable_abort_on_error() &
                bind(C, name="SLIC_disable_abort_on_error")
            implicit none
        end subroutine slic_disable_abort_on_error

        function c_is_abort_on_errors_enabled() &
                result(SHT_rv) &
                bind(C, name="SLIC_is_abort_on_errors_enabled")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL) :: SHT_rv
        end function c_is_abort_on_errors_enabled

        subroutine c_set_abort_on_warning(status) &
                bind(C, name="SLIC_set_abort_on_warning")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL), value, intent(IN) :: status
        end subroutine c_set_abort_on_warning

        subroutine slic_enable_abort_on_warning() &
                bind(C, name="SLIC_enable_abort_on_warning")
            implicit none
        end subroutine slic_enable_abort_on_warning

        subroutine slic_disable_abort_on_warning() &
                bind(C, name="SLIC_disable_abort_on_warning")
            implicit none
        end subroutine slic_disable_abort_on_warning

        function c_is_abort_on_warnings_enabled() &
                result(SHT_rv) &
                bind(C, name="SLIC_is_abort_on_warnings_enabled")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL) :: SHT_rv
        end function c_is_abort_on_warnings_enabled

        subroutine c_log_message(level, message, fileName, line, filter) &
                bind(C, name="SLIC_log_message")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: level
            character(kind=C_CHAR), intent(IN) :: message(*)
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: line
            logical(C_BOOL), value, intent(IN) :: filter
        end subroutine c_log_message

        subroutine c_log_message_bufferify(level, message, Lmessage, &
                fileName, LfileName, line, filter) &
                bind(C, name="SLIC_log_message_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: level
            character(kind=C_CHAR), intent(IN) :: message(*)
            integer(C_INT), value, intent(IN) :: Lmessage
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: LfileName
            integer(C_INT), value, intent(IN) :: line
            logical(C_BOOL), value, intent(IN) :: filter
        end subroutine c_log_message_bufferify

        subroutine slic_finalize() &
                bind(C, name="SLIC_finalize")
            implicit none
        end subroutine slic_finalize

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains

    function slic_is_initialized() &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        logical :: SHT_rv
        ! splicer begin function.is_initialized
        SHT_rv = c_is_initialized()
        ! splicer end function.is_initialized
    end function slic_is_initialized

    subroutine slic_create_logger(name, imask)
        use iso_c_binding, only : C_INT
        character(len=*), intent(IN) :: name
        character, value, intent(IN) :: imask
        ! splicer begin function.create_logger
        call c_create_logger_bufferify(name, len_trim(name, kind=C_INT), &
            imask)
        ! splicer end function.create_logger
    end subroutine slic_create_logger

    function slic_activate_logger(name) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        character(len=*), intent(IN) :: name
        logical :: SHT_rv
        ! splicer begin function.activate_logger
        SHT_rv = c_activate_logger_bufferify(name, &
            len_trim(name, kind=C_INT))
        ! splicer end function.activate_logger
    end function slic_activate_logger

    subroutine slic_get_active_logger_name(name)
        use iso_c_binding, only : C_INT
        character(len=*), intent(OUT) :: name
        ! splicer begin function.get_active_logger_name
        call c_get_active_logger_name_bufferify(name, &
            len(name, kind=C_INT))
        ! splicer end function.get_active_logger_name
    end subroutine slic_get_active_logger_name

    subroutine slic_set_abort_on_error(status)
        use iso_c_binding, only : C_BOOL
        logical, value, intent(IN) :: status
        ! splicer begin function.set_abort_on_error
        logical(C_BOOL) SH_status
        SH_status = status  ! coerce to C_BOOL
        call c_set_abort_on_error(SH_status)
        ! splicer end function.set_abort_on_error
    end subroutine slic_set_abort_on_error

    function slic_is_abort_on_errors_enabled() &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        logical :: SHT_rv
        ! splicer begin function.is_abort_on_errors_enabled
        SHT_rv = c_is_abort_on_errors_enabled()
        ! splicer end function.is_abort_on_errors_enabled
    end function slic_is_abort_on_errors_enabled

    subroutine slic_set_abort_on_warning(status)
        use iso_c_binding, only : C_BOOL
        logical, value, intent(IN) :: status
        ! splicer begin function.set_abort_on_warning
        logical(C_BOOL) SH_status
        SH_status = status  ! coerce to C_BOOL
        call c_set_abort_on_warning(SH_status)
        ! splicer end function.set_abort_on_warning
    end subroutine slic_set_abort_on_warning

    function slic_is_abort_on_warnings_enabled() &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        logical :: SHT_rv
        ! splicer begin function.is_abort_on_warnings_enabled
        SHT_rv = c_is_abort_on_warnings_enabled()
        ! splicer end function.is_abort_on_warnings_enabled
    end function slic_is_abort_on_warnings_enabled

    subroutine slic_log_message(level, message, fileName, line, filter)
        use iso_c_binding, only : C_BOOL, C_INT
        integer(C_INT), value, intent(IN) :: level
        character(len=*), intent(IN) :: message
        character(len=*), intent(IN) :: fileName
        integer(C_INT), value, intent(IN) :: line
        logical, value, intent(IN) :: filter
        ! splicer begin function.log_message
        logical(C_BOOL) SH_filter
        SH_filter = filter  ! coerce to C_BOOL
        call c_log_message_bufferify(level, message, &
            len_trim(message, kind=C_INT), fileName, &
            len_trim(fileName, kind=C_INT), line, SH_filter)
        ! splicer end function.log_message
    end subroutine slic_log_message

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module axom_slic
