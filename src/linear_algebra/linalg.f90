! Copyright 2023 Exciting developers
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
!> This file contain subroutines of linear algebra using MAGMA (GPU) or LAPACK (CPU) like library

!> This module contains the subroutines for the linear algebra
module m_linalg

    use m_constants, only: i32, i64, r64
#ifdef USE_GPU
    use iso_c_binding
    use magma2
    use m_gpu_magma_t, only: linalg_world_t, linalg_obj_t
#else
    use m_cpu_magma_t, only: linalg_world_t, linalg_obj_t
#endif
    
    implicit none

contains

    !> Inverts the A matrix and save the result to inverse_A
    !> @param[in] A - the complex matrix to invert
    !> @param[in] inverse_A - the inverse of A
    subroutine inverse_complex_LU(A, inverse_A)
        
        complex(r64), contiguous, intent(in)            :: A(:,:)
        complex(r64), target, allocatable, intent(out)  :: inverse_A(:,:)

        ! Locals
        type(linalg_world_t) :: world 
        type(linalg_obj_t)   :: ref_A
        type(linalg_obj_t)   :: ref_work

        ! LAPACK
        integer(i32) :: info, lwork
        integer(i32), allocatable :: ipiv(:)
        complex(r64), allocatable :: work(:)

#if !defined(USE_GPU)
        ! External
        external :: zgetri, zgetrf
#endif
        ! Allocate
        allocate(inverse_A, source=A)
        lwork = size(A,1)
        allocate(work(lwork), ipiv(lwork))

        ! If GPU init the queu otherwise this does nothing
        call world%init()

        ! In case of GPU transfer if not this call does nothing
        call ref_A%allocate_gpu(inverse_A)
        call ref_A%transfer_cpu_gpu(inverse_A, world%get_queue())
        call ref_work%allocate_gpu(work)

        ! Perform here the inversion
#ifdef USE_GPU
        call magma_zgetrf_gpu(ref_A%rows(), ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_zgetrf_gpu"
        call magma_zgetri_gpu(ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, ref_work%gpu_ptr(), lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_zgetri_gpu"
#else
        call zgetrf(ref_A%rows(), ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetrf"
        call zgetri(ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, ref_work%gpu_ptr(), lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetri"
#endif

        ! If GPU retrieve otherwise this does nothing
        call ref_A%transfer_gpu_cpu(inverse_A, world%get_queue())
        
        ! Clean
        call ref_A%destroy(world%get_queue())
        call ref_work%destroy(world%get_queue())
        call world%finish()

    end subroutine inverse_complex_LU

end module m_linalg