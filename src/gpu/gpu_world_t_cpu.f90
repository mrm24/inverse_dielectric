! Copyright (C) 2020-2024 GreenX library
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied. See the License for the specific language governing
! permissions and limitations under the License.

!> @file
!> This file contains things to handle host-device communication. CPU backend; so does nothing.

!> Module containing things to handle host-device communication
module idiel_gpu_world_t
    
    use iso_c_binding
    use idiel_constants, only: i32, i32, aip

    implicit none

    private
    public  gpu_world_t

    !> Type to handle the GPU Magma queue and world
    type gpu_world_t
        !> Device id taking care of control
        integer :: device = -1
        !> GPU queue
        type(C_ptr), private  :: queue    = C_null_ptr
    contains
        procedure, public :: init, finish, is_queue_set, get_queue, syncronize, get_device
    end type gpu_world_t

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                   CPU WORLD                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> This subroutine inits the dummy object
    !> @param[in] this - the GPU magma world to initialize
    !> @param[in] device - (optional) the device that should hold the queue (default: 0, which is the CPU)
    subroutine init(this, device)
        class(gpu_world_t), intent(inout) :: this
        integer, optional,  intent(in)    :: device
    end subroutine init

    !> This subroutine finishes the dummy object
    !> @param[in] this - the queue to finish
    subroutine finish(this)
        class(gpu_world_t), intent(inout) :: this
    end subroutine finish

    !> This function returns true 
    !> @param[in] this - the world to check if has inited queue
    pure function is_queue_set(this) result(answer)
        class(gpu_world_t), intent(in) :: this
        logical :: answer
        answer = .true.
    end function is_queue_set

    !> This is a getter for the queue, it returns the null pointer
    !> @param[in] this - the world object from which the queue is retrieved
    function get_queue(this) result(queue)

        class(gpu_world_t), target, intent(in) :: this
        type(C_ptr), pointer :: queue
        queue => this%queue

    end function get_queue 

    !> This syncronizes the queue
    !> @param[in] this - the world object from which the queue is syncronized
    subroutine syncronize(this)
        class(gpu_world_t), target, intent(in) :: this
    end subroutine syncronize

    !> This provides device id
    !> @param[in] this - return the associated device id
    pure function get_device(this) result(device)
        class(gpu_world_t), intent(in) :: this
        integer :: device
        device = this%device
    end function get_device

end module idiel_gpu_world_t
