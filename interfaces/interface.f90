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
!> This file contains a module for exposing the inverse_dielectric_t to C
!> 

!> This module defines the subroutines for exposing the inverse_dielectric_t to C
module m_inverse_dielectric_f90

    use iso_c_binding
    use m_constants, only: i64, r64
    use m_inverse_dielectric, only: inverse_dielectric_t
    
    implicit none
    
    public
    
contains
    
    !> Allocates inverse_dielectric_t
    subroutine allocate_inverse_dielectric_t(object_ptr) bind(C)
        type(c_ptr), intent(out) :: object_ptr
        type(inverse_dielectric_t), pointer :: object
        allocate(object)
        object_ptr = C_loc(object)
    end subroutine allocate_inverse_dielectric_t
    
    !> Deallocates the inverse_dielectric_t
    subroutine deallocate_inverse_dielectric_t(object_ptr) bind(C)
        type(c_ptr), intent(in) :: object_ptr
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        deallocate(object)
    end subroutine deallocate_inverse_dielectric_t
    
    subroutine init_common(object_ptr, lattice,  redpos, elements, nq) bind(C)
        type(c_ptr), intent(in)    :: object_ptr
        real(r64), intent(in)      :: lattice(3,3)
        real(r64), intent(in)      :: redpos(:,:) 
        integer(r64), intent(in)   :: elements(:)
        integer(i64), intent(in)   :: nq(3)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%init_common(lattice, redpos, elements, nq)
        
    end subroutine init_common
    
    subroutine set_dielectric_blocks_full(object_ptr, h, wl, wu, ib)
        type(c_ptr), intent(in) :: object_ptr
        complex(r64), target, intent(in)      :: h(:,:)
        complex(r64), target, intent(in)      :: wl(:,:)
        complex(r64), target, intent(in)      :: wu(:,:)
        complex(r64), target, intent(in)      :: ib(:,:)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%set_dielectric_blocks(h, wl, wu, ib)
        
    end subroutine set_dielectric_blocks_full
    
    subroutine set_dielectric_blocks_partial(object_ptr, h, wl, wu)
        type(c_ptr), intent(in) :: object_ptr
        complex(r64), target, intent(in)      :: h(:,:)
        complex(r64), target, intent(in)      :: wl(:,:)
        complex(r64), target, intent(in)      :: wu(:,:)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%set_dielectric_blocks(h, wl, wu)
        
    end subroutine set_dielectric_blocks_partial


    subroutine compute_anisotropic_avg(object_ptr)
        type(c_ptr), intent(in) :: object_ptr
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%compute_anisotropic_avg()
        
    end subroutine compute_anisotropic_avg
    
    subroutine invert_body(object_ptr, body)
        type(c_ptr), intent(in) :: object_ptr
        complex(r64), allocatable, intent(in) :: body(:,:)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)

        call object%invert_body(body)

    end subroutine invert_body

    function get_n_basis(object_ptr) result(nbasis) bind(C)
        type(c_ptr), intent(in) :: object_ptr
        integer(i64) :: nbasis
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)

        nbasis = object%get_n_basis()

    end function get_n_basis

    subroutine get_inverted_blocks(object_ptr, inv_head, inv_wingL, inv_wingU, inv_body) bind(C)
        type(c_ptr),  intent(in) :: object_ptr
        complex(r64), intent(inout) :: inv_head
        complex(r64), intent(inout) :: inv_wingL(:)
        complex(r64), intent(inout) :: inv_wingU(:)
        complex(r64), intent(inout) :: inv_body(:,:)
        
        type(inverse_dielectric_t), pointer :: object
        call C_F_pointer(object_ptr, object)
        
        call object%get_inverted_blocks(inv_head, inv_wingL, inv_wingU, inv_body)
        
    end subroutine get_inverted_blocks
    
end module m_inverse_dielectric_f90
