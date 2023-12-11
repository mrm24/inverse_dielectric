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
module idiel_linalg

    use idiel_constants, only: i32, i64, r64, zzero, zone
    use iso_c_binding
#ifdef USE_GPU
    use magma2
    use idiel_gpu_magma_t, only: linalg_world_t, linalg_obj_t
#else
    use idiel_cpu_magma_t, only: linalg_world_t, linalg_obj_t
#endif
    
    implicit none

contains

    !> Inverts the A matrix and save the result to inverse_A
    !> @param[in] A - the complex matrix to invert
    !> @param[in] inverse_A - the inverse of A
    !> @param[in] world - the linalg_world_t handler 
    subroutine inverse_complex_LU(A, inverse_A, world)
        
        complex(r64), contiguous, intent(in)            :: A(:,:)
        complex(r64), target, allocatable, intent(out)  :: inverse_A(:,:)
        type(linalg_world_t), intent(inout)             :: world

        ! Locals
        type(linalg_obj_t)   :: ref_A
        type(linalg_obj_t)   :: ref_work

        ! LAPACK
        integer(i32) :: info, lwork, nb, n
        integer(i32), allocatable :: ipiv(:)
        complex(r64), allocatable :: work(:)

#if !defined(USE_GPU)
        ! External
        external :: zgetri, zgetrf
#endif  

        ! Some constants
        n = size(A,1)
        ! Allocate
        allocate(inverse_A, source = A)

#ifdef USE_GPU
        nb = magma_get_zgetri_nb(n)
#else
        nb = 64
#endif
        lwork = nb * size(A,1)
        allocate(ipiv(size(A,1)))
        allocate(work(lwork))

        ! In case of GPU transfer if not this call does nothing
        call ref_A%allocate_gpu(inverse_A)
        call ref_A%transfer_cpu_gpu(inverse_A, world)
        call ref_work%allocate_gpu(work)

        ! Perform here the inversion (memory overflow is possible)
#ifdef USE_GPU
        !call magma_zgetrf(n, n, inverse_A, n, ipiv, info)
        call magma_zgetrf_gpu(ref_A%rows(), ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_zgetrf"
        call magma_zgetri_gpu(ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, ref_work%gpu_ptr(), lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_zgetri_gpu"
#else   
        call zgetrf(ref_A%rows(), ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetrf"
        call zgetri(ref_A%rows(), ref_A%gpu_ptr(), ref_A%rows(), ipiv, ref_work%gpu_ptr(), lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetri"
#endif

        ! If GPU retrieve otherwise this does nothing
        call ref_A%transfer_gpu_cpu(inverse_A, world)
        
        ! Clean
        call world%syncronize()
        call ref_A%destroy()
        call ref_work%destroy()

    end subroutine inverse_complex_LU

    !> It computes the auxiliary S auxiliary vector and the macroscopic dielectric matrix
    subroutine compute_S_and_L(Binv, wingL, S, ref_S, L, ref_L, world)
        complex(r64), intent(in)               :: Binv(:,:)
        complex(r64), intent(in)               :: wingL(:,:)
        complex(r64), allocatable, intent(out) :: S(:,:)
        type(linalg_obj_t), intent(out)        :: ref_S
        complex(r64), intent(inout)            :: L(:,:)
        type(linalg_obj_t), intent(out)        :: ref_L
        type(linalg_world_t), intent(inout)   :: world

        type(linalg_obj_t) :: ref_Binv
        type(linalg_obj_t) :: ref_wingL
        complex(r64) :: head(3,3)

        integer(i32) :: nb 
        
        nb = size(Binv,1)

        allocate(S, mold=wingL)
        head = L

        ! Allocate and transfer to the GPU/CPU
        call ref_Binv%allocate_gpu(Binv)
        call ref_Binv%transfer_cpu_gpu(Binv, world)
        call ref_wingL%allocate_gpu(wingL)
        call ref_wingL%transfer_cpu_gpu(wingL, world)
        call ref_L%allocate_gpu(L)
        call ref_L%transfer_cpu_gpu(L, world)
        call ref_S%allocate_gpu(S)

#ifdef USE_GPU
        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, zone, ref_Binv%gpu_ptr(), &
                         nb, ref_wingL%gpu_ptr(), nb, zzero, ref_S%gpu_ptr(), nb, world%get_queue())
                
        call magma_zgemm(MagmaConjTrans, MagmaNoTrans, 3, 3, nb, -zone, ref_wingL%gpu_ptr(), &
                        nb, ref_S%gpu_ptr(), nb, zone, ref_L%gpu_ptr(), 3, world%get_queue())
#else   
        call zgemm('n', 'n',nb, 3, nb, zone, ref_Binv%gpu_ptr(), &
                   nb, ref_wingL%gpu_ptr(), nb, zzero, ref_S%gpu_ptr(), nb)
        call zgemm('c', 'n',  3, 3, nb, -zone, ref_wingL%gpu_ptr(), &
                   nb, ref_S%gpu_ptr(), nb, zone, ref_L%gpu_ptr(), 3)
#endif     

        call ref_S%transfer_gpu_cpu(S, world)
        call ref_L%transfer_gpu_cpu(L, world)
        call world%syncronize()
        call ref_Binv%destroy()
        call ref_wingL%destroy()

    end subroutine compute_S_and_L

    !> It computes q^T \cdot L \cdot q
    !> @param[in] L     - macroscopic dielectric matrix
    !> @param[in] q     - the object q
    !> @param[in] ref_q - the linear algebra object containing q
    !> @param[in] world - the linalg_world_t handler 
    !> @returns qLq
    function compute_qLq(ref_L, q, ref_q, world) result(qLq)

        type(linalg_obj_t),   intent(inout)   :: ref_L
        complex(r64),         intent(in)      :: q(:,:)
        type(linalg_obj_t),   intent(inout)   :: ref_q
        type(linalg_world_t), intent(inout)   :: world

        complex(r64), allocatable :: qLq(:)

        complex(r64), allocatable :: Lq(:,:)
        integer(i32) :: nr 
        type(linalg_obj_t)        :: ref_Lq
        integer(i64) :: i

        nr = size(q,2)
        allocate(Lq, mold=q)
        call ref_Lq%allocate_gpu(Lq)

#ifdef USE_GPU
        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, 3, nr, 3, zone, ref_L%gpu_ptr(), & 
                         3, ref_q%gpu_ptr(), 3, zzero, ref_Lq%gpu_ptr(), 3, world%get_queue())
#else   
        call zgemm('n', 'n', 3, nr, 3, zone, L%gpu_ptr(), 3, ref_q%gpu_ptr(), 3, zzero, ref_Lq%gpu_ptr(), 3)
#endif  
        call ref_Lq%transfer_gpu_cpu(Lq, world)
        call world%syncronize()
        call ref_Lq%destroy()

        allocate(qLq(size(q,2)))

        !$omp parallel shared(qLq, Lq, q) private(i)
        !$omp do
        do i = 1, size(q,2)
            qLq(i) = dot_product(q(:,i),Lq(:,i))
        end do
        !$omp end do
        !$omp end parallel

    end function compute_qLq


    !> It computes Sq
    !> @param[in] ref_S - the linear algebra object containing S
    !> @param[in] ref_q - the linear algebra object containing q
    !> @param[in] world - the linalg_world_t handler 
    !> @returns qS
    function compute_qS(ref_S, ref_q, world) result(qS)
        type(linalg_obj_t),   intent(inout)   :: ref_S
        type(linalg_obj_t),   intent(inout)   :: ref_q
        type(linalg_world_t), intent(inout)   :: world

        complex(r64), allocatable :: qS(:,:)
        type(linalg_obj_t) :: ref_qS
        integer(i32) :: nr, nb 

        nr = ref_q%cols()
        nb = ref_S%rows()
        allocate(qS(nr,nb))
        call ref_qS%allocate_gpu(qS)

#ifdef USE_GPU
        call magma_zgemm(MagmaTrans, MagmaTrans, nr, nb, 3, zone, ref_q%gpu_ptr(), & 
                        3, ref_S%gpu_ptr(), nb, zzero, ref_qS%gpu_ptr(), nr, world%get_queue())
#else   
        call zgemm('t', 't', nr, nr, nb, 3, zone, ref_q%gpu_ptr(), & 
                  3, ref_S%gpu_ptr(), nb, zzero, ref_qS%gpu_ptr(), nr)
#endif  
        call ref_qS%transfer_gpu_cpu(qS, world)
        call world%syncronize()
        call ref_qS%destroy()

    end function compute_qS

end module idiel_linalg
