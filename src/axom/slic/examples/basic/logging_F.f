! Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

program slic_example
  ! This is a comment line; it is ignored by the compiler
  use iso_c_binding
  use axom_slic
  implicit none

  ! initialize slic environment
  !call slic_initialize()

  ! set the message level
  !call slic_set_logging_msg_level( message_debug )

  !call slic_disable_abort_on_error()
  ! register log stream (this doesn't exist???)
  !call slic_add_stream_to_all_levels( 'console' )
  !call slic_add_stream_to_all_levels( 'file:output.log' )

  ! log messages
  !call slic_log_message( slic_warning, "This is a warning", __FILE__, __LINE__ )
  call slic_log_message(message_debug, 'Here is a debug message!', 'logging_F.f', 33, .false.)
  call slic_log_message(message_info, 'Here is an info mesage!', 'logging_F.f', 33, .false.)
  call slic_log_message(message_warning, 'Here is a warning!', 'logging_F.f', 33, .false.)
  call slic_log_message(message_error, 'Here is an error message!', 'logging_F.f', 33, .false.)

  ! finalize slic environment
  call slic_finalize()

end program slic_example