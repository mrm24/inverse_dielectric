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
!> This file contain subroutines containing linear algebra using MAGMA (GPU) or LAPACK (CPU) like library
!> the routines are specifically written down case by case so the GPU sending and retrieving 
!> is minimized

!> This module contains the subroutines that build up vectors and blocks using linear algebra
module idiel_linalg

    use idiel_constants, only: i32, i64, r64, zzero, zone, fourpi
    use iso_c_binding
#ifdef USE_GPU
    use omp_lib
    use magma2
    use idiel_gpu_magma_t, only: linalg_world_t
#else
    use idiel_cpu_magma_t, only: linalg_world_t
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

        ! LAPACK
        integer(i32) :: info, lwork, nb, n
        integer(i32), allocatable :: ipiv(:)
        complex(r64), target, allocatable :: work(:)

#ifdef USE_GPU
        type(C_ptr) :: dA_ptr, dwork_ptr
#else
        external :: zgetri, zgetrf
#endif
        ! Some constants
        n = size(A,1)
        ! Allocate
        allocate(inverse_A, source = A)

#ifdef USE_GPU

        nb = magma_get_zgetri_nb(n)
        lwork = nb * size(A,1)
        allocate(ipiv(size(A,1)))
        allocate(work(lwork))

        !$omp target enter data map(to: inverse_A) map(alloc: work)

        ! Obtain device pointers
        dA_ptr    = omp_get_mapped_ptr(C_loc(inverse_A), world%get_device())
        dwork_ptr = omp_get_mapped_ptr(C_loc(work), world%get_device())

        call magma_zgetrf_gpu(n, n, dA_ptr, n, ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_zgetrf"
        call magma_zgetri_gpu(n, dA_ptr, n, ipiv, dwork_ptr, lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_zgetri_gpu"
        
        !$omp target update from(inverse_A)
        !$omp target exit data map(delete: inverse_A, work)
#else

        nb = 64
        lwork = nb * size(A,1)
        allocate(ipiv(size(A,1)))
        allocate(work(lwork))

        call zgetrf(n, n, inverse_A, n, ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetrf"
        call zgetri(n, inverse_A, n, ipiv, work, lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetri"
#endif

        deallocate(work, ipiv)

    end subroutine inverse_complex_LU

    !> It computes the auxiliary S and T vectors and the macroscopic dielectric matrix
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] S       - the auxiliary vector related to lower wing
    !> @param[in] wingU    - the upper wing
    !> @param[out] T       - the auxiliary vector related to upper wing
    !> @param[inout] L     - the macroscopic dielectric ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager   
    subroutine compute_auxiliary_and_macroscopic_3d_general(Binv, wingL, S, wingU, T, L, world)

        complex(r64), target, intent(inout)                        :: Binv(:,:)
        complex(r64), target, intent(inout)                        :: wingL(:,:)
        complex(r64), target, allocatable, intent(out)             :: S(:,:)
        complex(r64), target, intent(inout)                        :: wingU(:,:)
        complex(r64), target, allocatable, intent(out)             :: T(:,:)
        complex(r64), target, intent(inout)                        :: L(:,:)
        type(linalg_world_t), intent(inout)                        :: world

        integer(i32) :: nb

#ifdef USE_GPU
        type(C_ptr) :: L_dptr, S_dptr, T_dptr, Binv_dptr, wingL_dptr, wingU_dptr
#else
        ! External
        external :: zgemm
#endif

        nb = size(Binv,1)
        allocate(S(nb,3))
        allocate(T(3,nb))
#ifdef USE_GPU
        !$omp target enter data map(to: Binv, wingL, wingU, L) map(alloc: S, T)
        L_dptr     = omp_get_mapped_ptr(C_loc(L)    , world%get_device())
        S_dptr     = omp_get_mapped_ptr(C_loc(S)    , world%get_device())
        T_dptr     = omp_get_mapped_ptr(C_loc(T)    , world%get_device())
        Binv_dptr  = omp_get_mapped_ptr(C_loc(Binv) , world%get_device())
        wingL_dptr = omp_get_mapped_ptr(C_loc(wingL), world%get_device())
        wingU_dptr = omp_get_mapped_ptr(C_loc(wingU), world%get_device())

        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, zone, Binv_dptr, &
                    nb, wingL_dptr, nb, zzero, S_dptr, nb, world%get_queue())
        call magma_zgemm(MagmaTrans, MagmaNoTrans, 3, nb, nb, zone, wingU_dptr, &
                    nb, Binv_dptr, nb, zzero, T_dptr, 3, world%get_queue())
        call magma_zgemm(MagmaTrans, MagmaNoTrans, 3, 3, nb, -zone, wingU_dptr, &
                    nb, S_dptr, nb, zone, L_dptr, 3, world%get_queue())
        call world%syncronize()
        !$omp target update from(L, S)
        !$omp target exit data map(delete: Binv, wingU, wingL, L, S)
#else
        call zgemm('n', 'n', nb,  3, nb,  zone, Binv,  nb, wingL, nb, zzero, S, nb)
        call zgemm('t', 'n',  3, nb, nb,  zone, wingU, nb, Binv, nb, zzero, T, 3)
        call zgemm('t', 'n',  3,  3, nb, -zone, wingU, nb, S,    nb, zone,  L, 3)
#endif

    end subroutine compute_auxiliary_and_macroscopic_3d_general


    !> It computes the auxiliary S auxiliary vector and the macroscopic dielectric matrix
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] S       - the auxiliary vector related to lower wing
    !> @param[inout] L     - the macroscopic dielectric ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager
    subroutine compute_auxiliary_and_macroscopic_3d_hermitian(Binv, wingL, S, L, world)

        complex(r64), target, intent(inout)                        :: Binv(:,:)
        complex(r64), target, intent(inout)                        :: wingL(:,:)
        complex(r64), target, allocatable, intent(out)             :: S(:,:)
        complex(r64), target, intent(inout)                        :: L(:,:)
        type(linalg_world_t), intent(inout)                        :: world

        integer(i32) :: nb

#ifdef USE_GPU
        type(C_ptr) :: L_dptr, S_dptr, T_dptr, Binv_dptr, wingL_dptr
#else
        ! External
        external :: zgemm
#endif

        nb = size(Binv,1)
        allocate(S(nb,3))

#ifdef USE_GPU
        !$omp target enter data map(to: Binv, wingL, L)  map(alloc: S)
        L_dptr     = omp_get_mapped_ptr(C_loc(L)    , world%get_device())
        S_dptr     = omp_get_mapped_ptr(C_loc(S)    , world%get_device())
        Binv_dptr  = omp_get_mapped_ptr(C_loc(Binv) , world%get_device())
        wingL_dptr = omp_get_mapped_ptr(C_loc(wingL), world%get_device())

        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, zone, Binv_dptr, &
                        nb, wingL_dptr, nb, zzero, S_dptr, nb, world%get_queue())
        call magma_zgemm(MagmaConjTrans, MagmaNoTrans, 3, 3, nb, -zone, wingL_dptr, &
                        nb, S_dptr, nb, zone, L_dptr, 3, world%get_queue())
        call world%syncronize()
        !$omp target update from(L, S)
        !$omp target exit data map(delete: Binv, wingL, L, S)
#else
        call zgemm('n', 'n', nb, 3, nb,  zone, Binv,  nb, wingL, nb, zzero, S, nb)
        call zgemm('c', 'n',  3, 3, nb, -zone, wingL, nb, S,     nb, zone,  L, 3)
#endif


    end subroutine compute_auxiliary_and_macroscopic_3d_hermitian

    !> It computes 1 / q^T \cdot L \cdot q 
    !> @param[in] L      - macroscopic dielectric matrix (if GPU it must be loaded)
    !> @param[in] q      - the object q (if GPU it must be loaded)
    !> @param[in] world  - the linalg_world_t handler 
    !> @returns invqLq (if GPU it will also provide the GPU object)
    subroutine compute_inverse_head(L, q, world, invqLq)

        complex(r64),    target, allocatable, intent(inout) :: L(:,:)
        complex(r64),    target, allocatable, intent(in)    :: q(:,:)
        type(linalg_world_t), intent(inout)                 :: world
        complex(r64), allocatable, intent(inout)            :: invqLq(:)

        complex(r64), target, allocatable :: Lq(:,:)
        integer(i32) :: nr 
        integer(i64) :: i

#ifdef USE_GPU
        type(C_ptr) :: q_dptr, Lq_dptr, L_dptr 
#else
        ! External
        external :: zgemm
#endif 

        nr = size(q,2)
        allocate(Lq, mold=q)

#ifdef USE_GPU
        !$omp target enter data map(to: q, L) map(alloc: Lq)
        L_dptr   = omp_get_mapped_ptr(C_loc(L) , world%get_device())
        Lq_dptr  = omp_get_mapped_ptr(C_loc(Lq), world%get_device())
        q_dptr   = omp_get_mapped_ptr(C_loc(q) , world%get_device())

        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, 3, nr, 3, zone, L_dptr, & 
                         3, q_dptr, 3, zzero, Lq_dptr, 3, world%get_queue())
        call world%syncronize()
#else   
        call zgemm('n', 'n', 3, nr, 3, zone, L, 3, q, 3, zzero, Lq, 3)
#endif  

        if (allocated(invqLq)) deallocate(invqLq)
        allocate(invqLq(size(q,2)))
#ifdef USE_GPU 
        !$omp target teams distribute parallel do private(i) shared(invqLq, Lq, q) map(tofrom: invqLq)
        do i = 1, size(q,2)
                invqLq(i) = 1.0_r64 / dot_product(q(:,i),Lq(:,i))
        end do
        !$omp end target teams distribute parallel do
        !$omp target exit data map(delete: Lq, L, q)
#else       
        !$omp parallel shared(invqLq, Lq, q) private(i)
        !$omp do
        do i = 1, size(q,2)
                invqLq(i) = 1.0_r64 / dot_product(q(:,i),Lq(:,i))
        end do
        !$omp end do
        !$omp end parallel
#endif

        deallocate(L)

    end subroutine compute_inverse_head


    !> It computes Sq term of  \frac{ \mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}}
    !> @param[in] q - the linear algebra object containing q
    !> @param[in] S - S vector
    !> @param[in] world - the linalg_world_t handler 
    !> @param[out] qS   - (in GPU if compiled for)
    subroutine compute_inverse_wingL(q, S, world, qS)
        complex(r64), target, allocatable,  intent(inout)   :: q(:,:)
        complex(r64), target, allocatable,  intent(inout)   :: S(:,:)
        type(linalg_world_t),       intent(inout)           :: world
        complex(r64), target, allocatable,  intent(inout)   :: qS(:,:)

        integer(i32) :: nr, nb 

#if USE_GPU
        type(C_ptr) :: qS_dptr, S_dptr, q_dptr
#else
        ! External
        external :: zgemm
#endif 

        nr = size(q,2)
        nb = size(S,1)
        allocate(qS(nr,nb))

#ifdef USE_GPU
        !$omp target enter data map(always, alloc: qS) map(to: S, q)
        q_dptr   = omp_get_mapped_ptr(C_loc(q)  , world%get_device())
        qS_dptr  = omp_get_mapped_ptr(C_loc(qS) , world%get_device())
        S_dptr   = omp_get_mapped_ptr(C_loc(S)  , world%get_device())

        call magma_zgemm(MagmaTrans, MagmaTrans, nr, nb, 3, zone, q_dptr, & 
                        3, S_dptr, nb, zzero, qS_dptr, nr, world%get_queue())
        call world%syncronize()
        !$omp target update from(qS)
        !$omp target exit data map(delete: S, q, qS)
#else   
        call zgemm('T', 'T', nr, nb, 3, zone, q, 3, S, nb, zzero, qS, nr)
#endif  
        deallocate(S)

    end subroutine compute_inverse_wingL

    !> It computes Tq term of \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}}
    !> @param[in] ref_q - the linear algebra object containing q
    !> @param[in] ref_S - the linear algebra object containing T
    !> @param[in] world - the linalg_world_t handler 
    !> @param[out] qT   - (in GPU if compiled for)
    subroutine compute_inverse_wingU(q, T, world, qT)
        complex(r64), target, allocatable,  intent(inout)  :: q(:,:)
        complex(r64), target, allocatable,  intent(inout)  :: T(:,:)
        type(linalg_world_t),               intent(inout)  :: world
        complex(r64), target, allocatable,  intent(out)    :: qT(:,:)

        integer(i32) :: nr, nb 

#if USE_GPU
        type(C_ptr) :: qT_dptr, T_dptr, q_dptr
#else
        ! External
        external :: zgemm
#endif 

        nr = size(q,2)
        nb = size(T,2)
        allocate(qT(nr,nb))
#ifdef USE_GPU
        !$omp target enter data map(alloc: qT) map(to: T, q)
        q_dptr   = omp_get_mapped_ptr(C_loc(q)  , world%get_device())
        qT_dptr  = omp_get_mapped_ptr(C_loc(qT) , world%get_device())
        T_dptr   = omp_get_mapped_ptr(C_loc(T)  , world%get_device())

        call magma_zgemm(MagmaTrans, MagmaNoTrans, nr, nb, 3, zone, q_dptr, & 
                        3, T_dptr, 3, zzero, qT_dptr, nr, world%get_queue())
        call world%syncronize()
        !$omp target update from(qT)
        !$omp target exit data map(delete: T, q, qT)
#else   
        call zgemm('T', 'N', nr, nb, 3, zone, q, 3, T, 3, zzero, qT, nr)
#endif  
        deallocate(T)

    end subroutine compute_inverse_wingU

    !> It computes the auxiliary ag, bg vectors and the A tensor
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] ag      - the auxiliary vector related to lower wing
    !> @param[in] wingU    - the upper wing
    !> @param[out] bg      - the auxiliary vector related to upper wing
    !> @param[inout] A     - the A ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager   
    subroutine compute_auxiliary_and_A_2d_general(Binv, wingL, ag, wingU, bg, A, world)
        complex(r64), target, intent(inout)              :: Binv(:,:)
        complex(r64), target, intent(inout)              :: wingL(:,:)
        complex(r64), target, allocatable, intent(out)   :: ag(:,:)
        complex(r64), target, intent(inout)              :: wingU(:,:)
        complex(r64), target, allocatable, intent(out)   :: bg(:,:)
        complex(r64), target, intent(inout)              :: A(:,:)
        type(linalg_world_t), intent(inout)              :: world

        integer :: nb, i

#ifdef USE_GPU
        type(C_ptr) :: A_dptr, ag_dptr, bg_dptr, Binv_dptr, wingL_dptr, wingU_dptr
#else
        ! External
        external :: zgemm
#endif

        nb = size(Binv,1)
        allocate(ag(nb,3))

        ! Set arrays to the ones required by the formalism
        do i = 1, 3
            A(i,i) = zone - A(i,i)
        end do

        wingL = -cmplx(1.0_r64/sqrt(fourpi), 0.0_r64, r64) * wingL
        A     = -cmplx(1.0_r64/fourpi, 0.0_r64, r64) * A
        wingU =  -cmplx(1.0_r64/sqrt(fourpi), 0.0_r64, r64) * wingU
        allocate(bg(3,nb))

#ifdef USE_GPU
        !$omp target enter data map(to: A, wingL, Binv) map(alloc: ag)
        !$omp target enter data map(to: wingU) map(alloc: bg)
        A_dptr     = omp_get_mapped_ptr(C_loc(A)    , world%get_device())
        ag_dptr    = omp_get_mapped_ptr(C_loc(ag)   , world%get_device())
        bg_dptr    = omp_get_mapped_ptr(C_loc(bg)   , world%get_device())
        Binv_dptr  = omp_get_mapped_ptr(C_loc(Binv) , world%get_device())
        wingL_dptr = omp_get_mapped_ptr(C_loc(wingL), world%get_device())
        wingU_dptr = omp_get_mapped_ptr(C_loc(wingU), world%get_device())

        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, -zone, Binv_dptr, &
                    nb, wingL_dptr, nb, zzero, ag_dptr, nb, world%get_queue())
        call magma_zgemm(MagmaTrans, MagmaNoTrans, 3, nb, nb, -zone, wingU_dptr, &
                    nb, Binv_dptr, nb, zzero, bg_dptr, 3, world%get_queue())
        call magma_zgemm(MagmaTrans, MagmaNoTrans, 3, 3, nb, -zone, wingU_dptr, &
                    nb, ag_dptr, nb, zone, A_dptr, 3, world%get_queue())
        call world%syncronize()
        !$omp target update from(A)
        !$omp target exit data map(delete: Binv, wingU, wingL)
#else
        call zgemm('n', 'n', nb,  3, nb,  -zone, Binv,  nb, wingL, nb, zzero, ag, nb)
        call zgemm('t', 'n',  3, nb, nb,  -zone, wingU, nb, Binv, nb, zzero, bg, 3)
        call zgemm('t', 'n',  3,  3, nb,  -zone, wingU, nb, ag,   nb, zone,  A, 3)
#endif

    end subroutine compute_auxiliary_and_A_2d_general

    !> It computes the auxiliary ag vector and the A tensor
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] ag      - the auxiliary vector related to lower wing
    !> @param[inout] A     - the A ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager
    subroutine compute_auxiliary_and_A_2d_hermitian(Binv, wingL, ag, A, world)
        complex(r64), target, intent(inout)                      :: Binv(:,:)
        complex(r64), target, intent(inout)                      :: wingL(:,:)
        complex(r64), target, allocatable, intent(out)           :: ag(:,:)
        complex(r64), target, intent(inout)                      :: A(:,:)
        type(linalg_world_t), intent(inout)                      :: world

        integer :: nb, i

#ifdef USE_GPU
        type(C_ptr) :: A_dptr, ag_dptr, bg_dptr, Binv_dptr, wingL_dptr
#else
        ! External
        external :: zgemm
#endif

        nb = size(Binv,1)
        allocate(ag(nb,3))

        ! Set arrays to the ones required by the formalism
        do i = 1, 3
            A(i,i) = zone - A(i,i)
        end do

        wingL = -cmplx(1.0_r64/sqrt(fourpi), 0.0_r64, r64) * wingL
        A     = -cmplx(1.0_r64/fourpi, 0.0_r64, r64) * A

#ifdef USE_GPU
        !$omp target enter data map(to: A, wingL, Binv) map(alloc: ag)
        A_dptr     = omp_get_mapped_ptr(C_loc(A)    , world%get_device())
        ag_dptr    = omp_get_mapped_ptr(C_loc(ag)   , world%get_device())
        Binv_dptr  = omp_get_mapped_ptr(C_loc(Binv) , world%get_device())
        wingL_dptr = omp_get_mapped_ptr(C_loc(wingL), world%get_device())

        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, -zone, Binv_dptr, &
                        nb, wingL_dptr, nb, zzero, ag_dptr, nb, world%get_queue())
        call magma_zgemm(MagmaConjTrans, MagmaNoTrans, 3, 3, nb, -zone, wingL_dptr, &
                        nb, ag_dptr, nb, zone, A_dptr, 3, world%get_queue())
        call world%syncronize()
        !$omp target update from(A)
        !$omp target exit data map(delete: Binv, wingL)
#else
        call zgemm('n', 'n', nb, 3, nb, -zone, Binv,  nb, wingL, nb, zzero, ag, nb)
        call zgemm('c', 'n',  3, 3, nb, -zone, wingL, nb, ag,     nb, zone,  A, 3)
#endif

    end subroutine compute_auxiliary_and_A_2d_hermitian


    !> It computes q^T \cdot A \cdot q 
    !> @param[in] A      - A tensor (if GPU it must be loaded)
    !> @param[in] q      - the object q (if GPU it must be loaded)
    !> @param[in] world  - the linalg_world_t handler 
    !> @param[out] qAq   - (if GPU it will also provide the GPU object)
    subroutine compute_qAq(A, q, world, qAq)

        complex(r64), target, allocatable, intent(inout)   :: A(:,:)
        complex(r64), target, allocatable, intent(in)      :: q(:,:)
        type(linalg_world_t), intent(inout)                :: world
        complex(r64), allocatable, intent(out)             :: qAq(:)

        complex(r64), target, allocatable :: Aq(:,:)
        integer(i32) :: nr 
        integer(i64) :: i

#ifdef USE_GPU
        type(C_ptr) :: q_dptr, Aq_dptr, A_dptr 
#else
        ! External
        external :: zgemm
#endif 

        nr = size(q,2)
        allocate(Aq, mold=q)

#ifdef USE_GPU
        !$omp target enter data map(alloc: Aq)
        A_dptr   = omp_get_mapped_ptr(C_loc(A) , world%get_device())
        Aq_dptr  = omp_get_mapped_ptr(C_loc(Aq), world%get_device())
        q_dptr   = omp_get_mapped_ptr(C_loc(q) , world%get_device())
        call magma_zgemm(MagmaNoTrans, MagmaNoTrans, 3, nr, 3, zone, A_dptr, & 
                         3, q_dptr, 3, zzero, Aq_dptr, 3, world%get_queue())
        call world%syncronize()
#else   
        call zgemm('n', 'n', 3, nr, 3, zone, A, 3, q, 3, zzero, Aq, 3)
#endif  

        allocate(qAq(size(q,2)))
#ifdef USE_GPU 
        !$omp target enter data map(alloc: qAq)
        !$omp target teams distribute parallel do private(i) shared(qAq, Aq, q)
        do i = 1, size(q,2)
                qAq(i) = dot_product(q(:,i),Aq(:,i)) 
        end do
        !$omp end target teams distribute parallel do
        !$omp target exit data map(delete: Aq, A)
        !$omp target update from(qAq)
#else       
        !$omp parallel shared(qAq, Aq, q) private(i)
        !$omp do
        do i = 1, size(q,2)
                qAq(i) = dot_product(q(:,i),Aq(:,i))
        end do
        !$omp end do
        !$omp end parallel
#endif

        deallocate(A)

    end subroutine compute_qAq


end module idiel_linalg
