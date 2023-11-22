! Copyright 2023 EXCITING developers
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
!> This file contains a module for exposing the inverse_dielectric_t to python
!> 

!> This module defines the subroutines for exposing the inverse_dielectric_t to python
module m_inverse_dielectric_py

    use iso_c_binding
    use m_inverse_dielectric, only: inverse_dielectric_t
    
    implicit none
    
    public
    
contains
    
    !> Allocates inverse_dielectric_t
    subroutine allocate_inverse_dielectric_t(object_ptr)
        !f2py integer(8), intent(out) :: object_ptr
        type(c_ptr), intent(out) :: object_ptr
        type(inverse_dielectric_t), pointer :: object
        allocate(object)
        object_ptr = C_loc(object)
    end subroutine allocate_inverse_dielectric_t
    
    !> Deallocates the inverse_dielectric_t
    subroutine allocate_inverse_dielectric_t(object_ptr)
        !f2py integer(8), intent(in) :: object_ptr
        type(c_ptr), intent(in) :: object_ptr
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        deallocate(object)
    end subroutine allocate_inverse_dielectric_t
    
    subroutine init_common(object_ptr, lattice, redpos, elements, nq)
        !f2py integer(8), intent(in) :: object_ptr
        type(c_ptr), intent(in) :: object_ptr
        real(r64), intent(in) :: lattice(3,3)
        real(r64), allocatable, intent(in) :: redpos(:,:) 
        integer(r64), allocatable, intent(in) :: elements(:)
        integer(i64), intent(in) :: nq(3)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%init_common(lattice, redpos, elements, nq)
        
    end subroutine init_common
    
    subroutine set_dielectric_blocks(object_ptr, h, wl, wu, ib)
        !f2py integer(8), intent(in) :: object_ptr
        type(c_ptr), intent(in) :: object_ptr
        complex(r64), target, intent(in)      :: h(:,:)
        complex(r64), target, intent(in)      :: wl(:,:)
        complex(r64), target, intent(in)      :: wu(:,:)
        complex(r64), target, optional, intent(in) :: ib(:,:)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        if (present(ib)) then
            call object%set_dielectric_blocks(h, wl, wu, ib)
        else
            call object%set_dielectric_blocks(h, wl, wu)
        end if
        
    end set_dielectric_blocks
    
    subroutine compute_anisotropic_avg(object_ptr)
        !f2py integer(8), intent(in) :: object_ptr
        type(c_ptr), intent(in) :: object_ptr
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%compute_anisotropic_avg()
        
    end subroutine compute_anisotropic_avg
    
    subroutine invert_body(object_ptr, body)
        !f2py integer(8), intent(in) :: object_ptr
        type(c_ptr), intent(in) :: object_ptr
        complex(r64), allocatable, intent(in) :: body(:,:)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%compute_anisotropic_avg()
    end subroutine invert_body
    
end module m_inverse_dielectric_py
