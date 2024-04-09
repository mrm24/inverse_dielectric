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
!> This file contains things to handle host-device communication. GPU backend.

!> Module containing things to handle host-device communication
module idiel_gpu_world_t
    
    use iso_c_binding
#if defined(NVIDIAGPU) || defined(AMDGPU)
    use magma2
    use omp_lib
#endif
#if defined(INTELGPU)
    use omp_lib
#endif
    use idiel_constants, only: i32, i32, aip
    use m_gpu_register_fortran, only: device_host_register

    implicit none

    private
    public  gpu_world_t

    !> Type to handle the GPU Magma queue and world
    type gpu_world_t
        !> Device id taking care of control
        integer, private :: device = -1
        !> GPU queue
        type(C_ptr), private  :: queue = C_null_ptr
        !> The number of teams in the GPU
        integer, private :: num_teams
        !> Device host register
        type(device_host_register) :: register
    contains
        procedure, public :: init, finish, is_queue_set, get_queue, syncronize, get_device, get_num_teams
    end type gpu_world_t

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                   GPU WORLD                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> This subroutine inits the MAGMA world and the queue
    !> @param[in] this - the GPU magma world to initialize
    !> @param[in] device - (optional) the device that should hold the queue (if not compute)
    subroutine init(this, device)
            
        class(gpu_world_t), intent(inout) :: this
        integer, optional,  intent(in)    :: device

        integer :: num_teams

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Init MAGMA
        call magma_init()
#endif  

        ! Set device
        if (present(device)) then
            this%device = device
        else
#if defined(NVIDIAGPU) || defined(AMDGPU)
            call magma_get_device(this%device)
#endif 
#if defined(INTELGPU)
            this%device = omp_get_default_device()
            call omp_set_default_device(this%device)
#endif
        end if

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Init the GPU queue
        call magma_queue_create(this%device, this%queue)
#endif

        ! Init the number of teams 
        !$omp target teams map(from: num_teams)
        num_teams = omp_get_num_teams()
        !$omp end target teams
        this%num_teams = num_teams
    
        ! Init register
        call this%register%init()

    end subroutine init

    !> This subroutine finishes the MAGMA world and the queue
    !> @param[in] this - the queue to finish
    subroutine finish(this)
        
        class(gpu_world_t), intent(inout) :: this

#if defined(NVIDIAGPU) || defined(AMDGPU)
        ! Destroy queue
        call magma_queue_destroy(this%queue)
        
        ! Destroy world
        call magma_finalize()
#endif
        ! Clean register
        call this%register%finish()
    
    end subroutine finish

    !> This function returns true if the queue is inited
    !> @param[in] this - the world to check if has inited queue
    pure function is_queue_set(this) result(answer)

        class(gpu_world_t), intent(in) :: this
        logical :: answer

        answer = C_associated(this%queue)

    end function is_queue_set

    !> This is a getter for the queue
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
#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_queue_sync(this%queue)
#endif
    end subroutine syncronize

    !> This provides device id
    !> @param[in] this - return the associated device id
    pure function get_device(this) result(device)
        class(gpu_world_t), intent(in) :: this
        integer :: device
        device = this%device
    end function get_device
    
    !> This provides the number of teams
    !> @param[in] this - return the number of teams of the device
    pure function get_num_teams(this) result(num_teams)
        class(gpu_world_t), intent(in) :: this
        integer :: num_teams
        num_teams = this%num_teams
    end function get_num_teams

end module idiel_gpu_world_t
