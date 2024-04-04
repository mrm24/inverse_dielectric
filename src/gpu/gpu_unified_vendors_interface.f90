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
!> This file contains unified calls to GPU accelerated routines for NVIDIA, AMD, and INTEL cards
!> NOTE: AMD and NVIDIA support MPI-device aware calls through MAGMA library

!> Module containing unified calls to device accelerated routines
#if defined(INTELGPU)
    include "mkl_omp_offload.f90"
#endif
module idiel_unique_gpu_interface

! The iso_c_binding is included inside the pragmas because its positioning differs in both cases
#if defined(NVIDIAGPU) || defined(AMDGPU)
    use omp_lib
    use iso_c_binding
    use magma2
#endif
#if defined(INTELGPU)
    use omp_lib
    use iso_c_binding
    use onemkl_lapack_omp_offload_lp64
    use onemkl_blas_omp_offload_lp64
#endif
    use idiel_constants, only: i32, r32, r64
    use idiel_gpu_world_t, only: gpu_world_t

    implicit none

contains

    !!!!!!!!!!!!!!!  SINGLE PRECISION !!!!!!!!!!!!!!

    !> Complex single precision LU decomposition.
    !> All but ipiv needs to be allocated in the GPU.
    subroutine cgetrf_gpu(m, n, dA, lda, ipiv, info, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        type(C_ptr),  intent(inout)                     :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(gpu_world_t), intent(inout)                :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_cgetrf_gpu(m, n, dA, lda, ipiv, info)
#endif
#if defined(INTELGPU)
        complex(r32), pointer :: A(:,:)

        !We associate A with the device pointer dA, pointing to GPU memory.
        !This association allows the Fortran code to access and manipulate the data
        !on the GPU without explicitly transferring it back and forth between the CPU and GPU.
        call c_f_pointer(dA, A, [lda,n])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch
        call cgetrf(m, n, A, lda, ipiv, info)
        !$omp end target data

        nullify(A)
#endif
    end subroutine cgetrf_gpu

    !> Complex single inversion using LU decomposition
    !> ipiv must be allocated in the GPU
    subroutine cgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(C_ptr),  intent(inout)                     :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(C_ptr),  intent(inout)                     :: dwork
        integer(i32), intent(in)                        :: lwork
        type(gpu_world_t), intent(inout)                :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_cgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info)
#endif
#if defined(INTELGPU)

        complex(r32), pointer :: A(:,:), work(:)

        call c_f_pointer(dA, A, [lda*n])
        call c_f_pointer(dwork, work, [lwork])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch
        call cgetri(n, A, lda, ipiv, work, lwork, info)
        !$omp end target data

        nullify(A, work)
#endif

    end subroutine cgetri_gpu

    !> It provides the appropiate block size of cgetrf function
    function get_cgetri_nb_gpu(n) result(nblock)
        integer(i32), intent(in) :: n

        integer :: nblock

#if defined(NVIDIAGPU) || defined(AMDGPU)
        nblock = magma_get_cgetri_nb(n)
#endif
#if defined (INTELGPU)
        nblock = 64_i32
#endif

    end function get_cgetri_nb_gpu

    !> Complex double precision matrix-matrix product
    subroutine cgemm_gpu(transa, transb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world)
        character, intent(in)         :: transa
        character, intent(in)         :: transb
        integer(i32), intent(in)      :: m
        integer(i32), intent(in)      :: n
        integer(i32), intent(in)      :: k
        complex(r32), intent(in)      :: alpha
        type(c_ptr),  intent(in)      :: da
        integer(i32), intent(in)      :: lda
        type(c_ptr),  intent(in)      :: db
        integer(i32), intent(in)      :: ldb
        complex(r32), intent(in)      :: beta
        type(c_ptr),  intent(inout)   :: dc
        integer(i32), intent(in)      :: ldc
        type(gpu_world_t), intent(in) :: world

        integer(c_int) :: opA, opB

#if defined(NVIDIAGPU) || defined(AMDGPU)

        if (transa == 'N' .or. transa == 'n') opA = MagmaNoTrans
        if (transa == 'T' .or. transa == 't') opA = MagmaTrans
        if (transa == 'C' .or. transa == 'c') opA = MagmaConjTrans

        if (transb == 'N' .or. transb == 'n') opB = MagmaNoTrans
        if (transb == 'T' .or. transb == 't') opB = MagmaTrans
        if (transb == 'C' .or. transb == 'c') opB = MagmaConjTrans

        call magma_cgemm(opA, opB, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r32), pointer :: A(:,:), B(:,:), C(:,:)
        integer(c_int) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(da, A, [lda,ka])
        call c_f_pointer(db, B, [ldb,kb])
        call c_f_pointer(dc, C, [ldc,n])

        !$omp target data
        !$omp dispatch
        call cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        !$omp end target data

        nullify(A, B, C)

#endif
    end subroutine cgemm_gpu

#if !defined(USE_SINGLE_PRECISION)
    !> Complex double precision LU decomposition.
    !> All but ipiv needs to be allocated in the GPU.
    subroutine zgetrf_gpu(m, n, dA, lda, ipiv, info, world)
        integer(i32), intent(in)                        :: m
        integer(i32), intent(in)                        :: n
        type(C_ptr),  intent(inout)                     :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(gpu_world_t), intent(inout)                :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zgetrf_gpu(m, n, dA, lda, ipiv, info)
#endif
#if defined(INTELGPU)
        complex(r64), pointer :: A(:,:)
        ! We associate A with the device pointer dA, pointing to GPU memory.
        ! This association allows the Fortran code to access and manipulate the data
        ! on the GPU without explicitly transferring it back and forth between the CPU and GPU.
        call c_f_pointer(dA, A, [lda,n])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch
        call zgetrf(m, n, A, lda, ipiv, info)
        !$omp end target data

        nullify(A)
#endif
    end subroutine zgetrf_gpu

    !> Complex double inversion using LU decomposition
    !> ipiv must be allocated in the GPU
    subroutine zgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info, world)
        integer(i32), intent(in)                        :: n
        type(C_ptr),  intent(inout)                     :: dA
        integer(i32), intent(in)                        :: lda
        integer(i32), contiguous, target, intent(inout) :: ipiv(:)
        integer(i32), intent(out)                       :: info
        type(C_ptr),  intent(inout)                     :: dwork
        integer(i32), intent(in)                        :: lwork
        type(gpu_world_t), intent(inout)                :: world

#if defined(NVIDIAGPU) || defined(AMDGPU)
        call magma_zgetri_gpu(n, dA, lda, ipiv, dwork, lwork, info)
#endif
#if defined(INTELGPU)
        complex(r64), pointer :: A(:,:), work(:)

        call c_f_pointer(dA, A, [lda*n])
        call c_f_pointer(dwork, work, [lwork])

        !$omp target data map(tofrom: info, ipiv)
        !$omp dispatch
        call zgetri(n, A, lda, ipiv, work, lwork, info)
        !$omp end target data

        nullify(A, work)
#endif

    end subroutine zgetri_gpu

    !> It provides the appropiate block size of zgetrf function
    function get_zgetri_nb_gpu(n) result(nblock)
        integer(i32), intent(in) :: n

        integer :: nblock

#if defined(NVIDIAGPU) || defined(AMDGPU)
        nblock = magma_get_zgetri_nb(n)
#endif
#if defined (INTELGPU)
        nblock = 64_i32
#endif

    end function get_zgetri_nb_gpu

    !> Complex double precision matrix-matrix product
    subroutine zgemm_gpu(transa, transb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world)
        character, intent(in)         :: transa
        character, intent(in)         :: transb
        integer(i32), intent(in)      :: m
        integer(i32), intent(in)      :: n
        integer(i32), intent(in)      :: k
        complex(r64), intent(in)      :: alpha
        type(c_ptr),  intent(inout)   :: da
        integer(i32), intent(in)      :: lda
        type(c_ptr),  intent(inout)   :: db
        integer(i32), intent(in)      :: ldb
        complex(r64), intent(in)      :: beta
        type(c_ptr),  intent(inout)   :: dc
        integer(i32), intent(in)      :: ldc
        type(gpu_world_t), intent(in) :: world

        integer(c_int) :: opA, opB
        complex(r64), pointer :: A(:), B(:), C(:)
        integer(c_int) :: ka, kb
#if defined(NVIDIAGPU) || defined(AMDGPU)

        if (transa == 'N' .or. transa == 'n') opA = MagmaNoTrans
        if (transa == 'T' .or. transa == 't') opA = MagmaTrans
        if (transa == 'C' .or. transa == 'c') opA = MagmaConjTrans

        if (transb == 'N' .or. transb == 'n') opB = MagmaNoTrans
        if (transb == 'T' .or. transb == 't') opB = MagmaTrans
        if (transb == 'C' .or. transb == 'c') opB = MagmaConjTrans

        call magma_zgemm(opA, opB, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, world%get_queue())
#endif
#if defined(INTELGPU)

        complex(r64), pointer :: A(:,:), B(:,:), C(:,:)
        integer(c_int) :: ka, kb

        if (transa == 'N' .or. transa == 'n') ka = k
        if (transa == 'T' .or. transa == 't') ka = m
        if (transa == 'C' .or. transa == 'c') ka = m

        if (transb == 'N' .or. transb == 'n') kb = n
        if (transb == 'T' .or. transb == 't') kb = k
        if (transb == 'C' .or. transb == 'c') kb = k

        call c_f_pointer(da, A, [lda,ka])
        call c_f_pointer(db, B, [ldb,kb])
        call c_f_pointer(dc, C, [ldc,n])

        !$omp target data
        !$omp dispatch
        call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        !$omp end target data

        nullify(A, B, C)

#endif
    end subroutine zgemm_gpu

#endif

end module idiel_unique_gpu_interface
