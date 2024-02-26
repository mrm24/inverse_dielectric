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
!> This file contain subroutines containing linear algebra using MAGMA (GPU) or LAPACK (CPU) like library
!> the routines are specifically written down case by case so the GPU sending and retrieving 
!> is minimized

!> This module contains the subroutines that build up vectors and blocks using linear algebra
module idiel_linalg

    use idiel_constants, only: i32, i32, aip, zzero, zone, fourpi
    use iso_c_binding
#ifdef USE_GPU
    use omp_lib
    use magma2
    use idiel_gpu_magma_t, only: linalg_world_t
#else
    use idiel_cpu_magma_t, only: linalg_world_t
#endif


! Set some macros so the proper linear algebra are called
#ifdef USE_SINGLE_PRECISION
#define _getrf             cgetrf
#define _getri             cgetri
#define _gemm              cgemm
#define _getrf_gpu         magma_cgetrf_gpu
#define _getri_gpu         magma_cgetri_gpu
#define _get_getri_nb_gpu  magma_get_cgetri_nb
#define _gemm_gpu          magma_cgemm
#else
#define _getrf             zgetrf
#define _getri             zgetri
#define _gemm              zgemm
#define _getrf_gpu         magma_zgetrf_gpu
#define _getri_gpu         magma_zgetri_gpu
#define _get_getri_nb_gpu  magma_get_zgetri_nb
#define _gemm_gpu          magma_zgemm
#endif

    
    implicit none

contains

    !> Inverts the A matrix and save the result to inverse_A
    !> @param[in] A - the complex matrix to invert
    !> @param[in] inverse_A - the inverse of A
    !> @param[in] world - the linalg_world_t handler
    subroutine inverse_complex_LU(A, inverse_A, world)

        complex(aip), contiguous, intent(in)            :: A(:,:)
        complex(aip), target, allocatable, intent(out)  :: inverse_A(:,:)
        type(linalg_world_t), intent(inout)             :: world

        ! LAPACK
        integer(i32) :: info, lwork, nb, n
        integer(i32), allocatable :: ipiv(:)
        complex(aip), target, allocatable :: work(:)

#ifdef USE_GPU
        type(C_ptr) :: dA_ptr, dwork_ptr
#else
        external :: _getri, _getrf
#endif
        ! Some constants
        n = size(A,1)
        ! Allocate
        allocate(inverse_A, source = A)

#ifdef USE_GPU

        nb = _get_getri_nb_gpu(n)
        lwork = nb * size(A,1)
        allocate(ipiv(size(A,1)))

        !Allocate in device inverse_A, associate to host inverse_A, and copy from host to device
        call world%register%alloc("inverse_A", n*n*c_sizeof(zzero), world%get_device())
        call world%register%assoc("inverse_A", C_loc(inverse_A))
        call world%register%to_device("inverse_A")
        ! Allocate in device work array
        call world%register%alloc("work", lwork*c_sizeof(zzero), world%get_device())

        ! Obtain device pointers
        dA_ptr    = world%register%device_ptr("inverse_A")
        dwork_ptr = world%register%device_ptr("work")

        call _getrf_gpu(n, n, dA_ptr, n, ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_Xgetrf"
        call _getri_gpu(n, dA_ptr, n, ipiv, dwork_ptr, lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling magma_Xgetri"

        call world%register%from_device("inverse_A")
        call world%register%remove("inverse_A")
        call world%register%remove("work")

        deallocate(ipiv)
#else

        nb = 64
        lwork = nb * size(A,1)
        allocate(ipiv(size(A,1)))
        allocate(work(lwork))

        call _getrf(n, n, inverse_A, n, ipiv, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetrf"
        call _getri(n, inverse_A, n, ipiv, work, lwork, info)
        if (info /= 0) error stop "inverse_complex_LU: error calling zgetri"

        deallocate(work, ipiv)
#endif

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

        complex(aip), target, intent(inout)                        :: Binv(:,:)
        complex(aip), target, intent(inout)                        :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)             :: S(:,:)
        complex(aip), target, intent(inout)                        :: wingU(:,:)
        complex(aip), target, allocatable, intent(out)             :: T(:,:)
        complex(aip), target, intent(inout)                        :: L(:,:)
        type(linalg_world_t), intent(inout)                        :: world

        integer(i32) :: nb

#ifdef USE_GPU
        type(C_ptr) :: L_dptr, S_dptr, T_dptr, Binv_dptr, wingL_dptr, wingU_dptr
#else
        ! External
        external :: _gemm
#endif

        nb = size(Binv,1)
        allocate(S(nb,3))
        allocate(T(3,nb))
#ifdef USE_GPU
        
        ! GPU allocation
        call world%register%alloc("Binv", nb*nb*c_sizeof(zzero), world%get_device())
        call world%register%alloc("wingL", 3*nb*c_sizeof(zzero), world%get_device())
        call world%register%alloc("wingU", 3*nb*c_sizeof(zzero), world%get_device())
        call world%register%alloc("L", 9*c_sizeof(zzero), world%get_device())
        call world%register%alloc("S", 3*nb*c_sizeof(zzero), world%get_device())
        call world%register%alloc("T", 3*nb*c_sizeof(zzero), world%get_device())

        ! GPU-CPU association
        call world%register%assoc("Binv",  C_loc(Binv))
        call world%register%assoc("wingL", C_loc(wingL))
        call world%register%assoc("wingU", C_loc(wingU))
        call world%register%assoc("L",     C_loc(L))
        call world%register%assoc("S",     C_loc(S))
        call world%register%assoc("T",     C_loc(T))

        ! CPU -> GPU copy
        call world%register%to_device("Binv")
        call world%register%to_device("wingL")
        call world%register%to_device("wingU")
        call world%register%to_device("L")

        ! Get device pointers
        L_dptr     = world%register%device_ptr("L")
        S_dptr     = world%register%device_ptr("S")
        T_dptr     = world%register%device_ptr("T")
        Binv_dptr  = world%register%device_ptr("Binv")
        wingL_dptr = world%register%device_ptr("wingL")
        wingU_dptr = world%register%device_ptr("wingU")

        call _gemm_gpu(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, zone, Binv_dptr, &
                    nb, wingL_dptr, nb, zzero, S_dptr, nb, world%get_queue())
        call _gemm_gpu(MagmaTrans, MagmaNoTrans, 3, nb, nb, zone, wingU_dptr, &
                    nb, Binv_dptr, nb, zzero, T_dptr, 3, world%get_queue())
        call _gemm_gpu(MagmaTrans, MagmaNoTrans, 3, 3, nb, -zone, wingU_dptr, &
                    nb, S_dptr, nb, zone, L_dptr, 3, world%get_queue())
        call world%syncronize()
        
        ! Update L 
        call world%register%from_device("L")

        ! Deassociate from device but keep them there
        call world%register%deassoc("L")
        call world%register%deassoc("S")
        call world%register%deassoc("T")

        ! Delete from device
        call world%register%remove("Binv")
        call world%register%remove("wingL")
        call world%register%remove("wingU")

#else
        call _gemm('n', 'n', nb,  3, nb,  zone, Binv,  nb, wingL, nb, zzero, S, nb)
        call _gemm('t', 'n',  3, nb, nb,  zone, wingU, nb, Binv, nb, zzero, T, 3)
        call _gemm('t', 'n',  3,  3, nb, -zone, wingU, nb, S,    nb, zone,  L, 3)
#endif

    end subroutine compute_auxiliary_and_macroscopic_3d_general


    !> It computes the auxiliary S auxiliary vector and the macroscopic dielectric matrix
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] S       - the auxiliary vector related to lower wing
    !> @param[inout] L     - the macroscopic dielectric ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager
    subroutine compute_auxiliary_and_macroscopic_3d_hermitian(Binv, wingL, S, L, world)

        complex(aip), target, intent(inout)                        :: Binv(:,:)
        complex(aip), target, intent(inout)                        :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)             :: S(:,:)
        complex(aip), target, intent(inout)                        :: L(:,:)
        type(linalg_world_t), intent(inout)                        :: world

        integer(i32) :: nb

#ifdef USE_GPU
        type(C_ptr) :: L_dptr, S_dptr, T_dptr, Binv_dptr, wingL_dptr
#else
        ! External
        external :: _gemm
#endif
        nb = size(Binv,1)
        allocate(S(nb,3))

#ifdef USE_GPU
        ! GPU allocation
        call world%register%alloc("Binv", nb*nb*c_sizeof(zzero), world%get_device())
        call world%register%alloc("wingL", 3*nb*c_sizeof(zzero), world%get_device())
        call world%register%alloc("L", 9*c_sizeof(zzero), world%get_device())
        call world%register%alloc("S", 3*nb*c_sizeof(zzero), world%get_device())

        ! GPU-CPU association
        call world%register%assoc("Binv",  C_loc(Binv))
        call world%register%assoc("wingL", C_loc(wingL))
        call world%register%assoc("L",     C_loc(L))
        call world%register%assoc("S",     C_loc(S))

        ! CPU -> GPU copy
        call world%register%to_device("Binv")
        call world%register%to_device("wingL")
        call world%register%to_device("L")

        ! Get device pointers
        L_dptr     = world%register%device_ptr("L")
        S_dptr     = world%register%device_ptr("S")
        Binv_dptr  = world%register%device_ptr("Binv")
        wingL_dptr = world%register%device_ptr("wingL")

        call _gemm_gpu(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, zone, Binv_dptr, &
                        nb, wingL_dptr, nb, zzero, S_dptr, nb, world%get_queue())
        call _gemm_gpu(MagmaConjTrans, MagmaNoTrans, 3, 3, nb, -zone, wingL_dptr, &
                        nb, S_dptr, nb, zone, L_dptr, 3, world%get_queue())
        call world%syncronize()
        
        ! Update L 
        call world%register%from_device("L")

        ! Deassociate from device but keep them there
        call world%register%deassoc("L")
        call world%register%deassoc("S")

        ! Delete from device
        call world%register%remove("Binv")
        call world%register%remove("wingL")

#else
        call _gemm('n', 'n', nb, 3, nb,  zone, Binv,  nb, wingL, nb, zzero, S, nb)
        call _gemm('c', 'n',  3, 3, nb, -zone, wingL, nb, S,     nb, zone,  L, 3)
#endif


    end subroutine compute_auxiliary_and_macroscopic_3d_hermitian

    !> It computes 1 / q^T \cdot L \cdot q 
    !> @param[in] L      - macroscopic dielectric matrix (if GPU it must be loaded)
    !> @param[in] q      - the object q (if GPU it must be loaded)
    !> @param[in] world  - the linalg_world_t handler 
    !> @returns invqLq (if GPU it will also provide the GPU object)
    subroutine compute_inverse_head(L, q, world, invqLq)

        complex(aip),    target, allocatable, intent(inout) :: L(:,:)
        complex(aip),    target, allocatable, intent(in)    :: q(:,:)
        type(linalg_world_t), intent(inout)                 :: world
        complex(aip),    target, allocatable, intent(out)   :: invqLq(:)

        complex(aip), target, allocatable :: Lq(:,:)
        integer(i32) :: nr 
        integer(i32) :: i

#ifdef USE_GPU
        type(C_ptr) :: q_dptr, Lq_dptr, L_dptr 
        integer     :: device
#else
        ! External
        external :: _gemm
#endif 

        nr = size(q,2)
        allocate(Lq, mold=q)

#ifdef USE_GPU

        ! Create Lq in the device
        call world%register%alloc("Lq", size(q) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("Lq", C_loc(Lq))

        ! Reassociate q and L to host memory
        call world%register%assoc("q", C_loc(q))
        call world%register%assoc("L", C_loc(L))

        ! Update L in device to the symmetrized one
        call world%register%to_device("L")
        
        ! Get device pointers
        L_dptr   = world%register%device_ptr("L")
        Lq_dptr  = world%register%device_ptr("Lq")
        q_dptr   = world%register%device_ptr("q")

        call _gemm_gpu(MagmaNoTrans, MagmaNoTrans, 3, nr, 3, zone, L_dptr, &
                         3, q_dptr, 3, zzero, Lq_dptr, 3, world%get_queue())
        call world%syncronize()
#else   
        call _gemm('n', 'n', 3, nr, 3, zone, L, 3, q, 3, zzero, Lq, 3)
#endif  

        if (allocated(invqLq)) deallocate(invqLq)
        allocate(invqLq(size(q,2)))
#ifdef USE_GPU
        call world%register%alloc("invqLq", size(invqLq) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("invqLq", C_loc(invqLq))
        !$omp target teams distribute parallel do private(i) 
        do i = 1, size(q,2)
                invqLq(i) = 1.0_aip / dot_product(q(:,i),Lq(:,i))
        end do
        !$omp end target teams distribute parallel do
        
        call world%register%from_device("invqLq")
        
        call world%register%deassoc("invqLq")
        
        call world%register%deassoc("q")
        call world%register%remove("Lq")
        call world%register%remove("L")
#else       
        !$omp parallel shared(invqLq, Lq, q) private(i)
        !$omp do
        do i = 1, size(q,2)
                invqLq(i) = 1.0_aip / dot_product(q(:,i),Lq(:,i))
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
        complex(aip), target, allocatable,  intent(inout)   :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)   :: S(:,:)
        type(linalg_world_t),       intent(inout)           :: world
        complex(aip), target, allocatable,  intent(inout)   :: qS(:,:)

        integer(i32) :: nr, nb 

#if USE_GPU
        type(C_ptr) :: qS_dptr, S_dptr, q_dptr
#else
        ! External
        external :: _gemm
#endif 

        nr = size(q,2)
        nb = size(S,1)
        allocate(qS(nr,nb))

#ifdef USE_GPU

        ! Create Lq in the device
        call world%register%alloc("qS", size(qS) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("qS", C_loc(qS))

        ! Reassociate q and L to host memory
        call world%register%assoc("q", C_loc(q))
        call world%register%assoc("S", C_loc(S))

        ! Get device pointers
        q_dptr   = world%register%device_ptr("q")
        qS_dptr  = world%register%device_ptr("qS")
        S_dptr   = world%register%device_ptr("S")

        call _gemm_gpu(MagmaTrans, MagmaTrans, nr, nb, 3, zone, q_dptr, &
                        3, S_dptr, nb, zzero, qS_dptr, nr, world%get_queue())
        call world%syncronize()
        
        call world%register%from_device("qS")
        call world%register%remove("S")
        call world%register%deassoc("qS")
        call world%register%deassoc("q")
#else   
        call _gemm('T', 'T', nr, nb, 3, zone, q, 3, S, nb, zzero, qS, nr)
#endif  
        deallocate(S)

    end subroutine compute_inverse_wingL

    !> It computes Tq term of \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}}
    !> @param[in] ref_q - the linear algebra object containing q
    !> @param[in] ref_S - the linear algebra object containing T
    !> @param[in] world - the linalg_world_t handler 
    !> @param[out] qT   - (in GPU if compiled for)
    subroutine compute_inverse_wingU(q, T, world, qT)
        complex(aip), target, allocatable,  intent(inout)  :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)  :: T(:,:)
        type(linalg_world_t),               intent(inout)  :: world
        complex(aip), target, allocatable,  intent(out)    :: qT(:,:)

        integer(i32) :: nr, nb 

#if USE_GPU
        type(C_ptr) :: qT_dptr, T_dptr, q_dptr
#else
        ! External
        external :: _gemm
#endif 

        nr = size(q,2)
        nb = size(T,2)
        allocate(qT(nr,nb))
#ifdef USE_GPU
        call world%register%alloc("qT", size(qT) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("qT", C_loc(qT))

        ! Reassociate q and L to host memory
        call world%register%assoc("q", C_loc(q))
        call world%register%assoc("T", C_loc(T))

        ! Get device pointers
        q_dptr   = world%register%device_ptr("q")
        qT_dptr  = world%register%device_ptr("qT")
        T_dptr   = world%register%device_ptr("T")

        call _gemm_gpu(MagmaTrans, MagmaNoTrans, nr, nb, 3, zone, q_dptr, &
                        3, T_dptr, 3, zzero, qT_dptr, nr, world%get_queue())
        
        call world%syncronize()

        call world%register%from_device("qT")
        call world%register%remove("T")
        call world%register%deassoc("qT")
        call world%register%deassoc("q")
#else   
        call _gemm('T', 'N', nr, nb, 3, zone, q, 3, T, 3, zzero, qT, nr)
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
        complex(aip), target, intent(inout)              :: Binv(:,:)
        complex(aip), target, intent(inout)              :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)   :: ag(:,:)
        complex(aip), target, intent(inout)              :: wingU(:,:)
        complex(aip), target, allocatable, intent(out)   :: bg(:,:)
        complex(aip), target, intent(inout)              :: A(:,:)
        type(linalg_world_t), intent(inout)              :: world

        integer :: nb, i

#ifdef USE_GPU
        type(C_ptr) :: A_dptr, ag_dptr, bg_dptr, Binv_dptr, wingL_dptr, wingU_dptr
#else
        ! External
        external :: _gemm
#endif

        nb = size(Binv,1)
        allocate(ag(nb,3))

        ! Set arrays to the ones required by the formalism
        do i = 1, 3
            A(i,i) = zone - A(i,i)
        end do

        wingL = -cmplx(1.0_aip/sqrt(fourpi), 0.0_aip, aip) * wingL
        A     = -cmplx(1.0_aip/fourpi, 0.0_aip, aip) * A
        wingU = -cmplx(1.0_aip/sqrt(fourpi), 0.0_aip, aip) * wingU
        allocate(bg(3,nb))

#ifdef USE_GPU
        ! Alloc in the device
        call world%register%alloc("A"    , size(A) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("ag"   , size(ag) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("bg"   , size(bg) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("Binv" , size(Binv) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("wingL", size(wingL) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("wingU", size(wingU) * c_sizeof(zzero), world%get_device())

        ! Associate
        call world%register%assoc("A"    , C_loc(A)    )
        call world%register%assoc("ag"   , C_loc(ag)   )
        call world%register%assoc("bg"   , C_loc(bg)   )
        call world%register%assoc("Binv" , C_loc(Binv) )
        call world%register%assoc("wingL", C_loc(wingL))
        call world%register%assoc("wingU", C_loc(wingU))

        ! Host device movement
        call world%register%to_device("A")
        call world%register%to_device("Binv")
        call world%register%to_device("wingL")
        call world%register%to_device("wingU")

        ! Get device pointers
        A_dptr     = world%register%device_ptr("A"    )
        ag_dptr    = world%register%device_ptr("ag"   )
        bg_dptr    = world%register%device_ptr("bg"   )
        Binv_dptr  = world%register%device_ptr("Binv" )
        wingL_dptr = world%register%device_ptr("wingL")
        wingU_dptr = world%register%device_ptr("wingU")

        call _gemm_gpu(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, -zone, Binv_dptr, &
                    nb, wingL_dptr, nb, zzero, ag_dptr, nb, world%get_queue())
        call _gemm_gpu(MagmaTrans, MagmaNoTrans, 3, nb, nb, -zone, wingU_dptr, &
                    nb, Binv_dptr, nb, zzero, bg_dptr, 3, world%get_queue())
        call _gemm_gpu(MagmaTrans, MagmaNoTrans, 3, 3, nb, -zone, wingU_dptr, &
                    nb, ag_dptr, nb, zone, A_dptr, 3, world%get_queue())
        call world%syncronize()

        ! Update A
        call world%register%from_device("A")

        ! Deassociate from device but keep them there
        call world%register%deassoc("A")
        call world%register%deassoc("ag")
        call world%register%deassoc("bg")

        ! Delete from device
        call world%register%remove("Binv")
        call world%register%remove("wingL")
        call world%register%remove("wingU")

#else
        call _gemm('n', 'n', nb,  3, nb,  -zone, Binv,  nb, wingL, nb, zzero, ag, nb)
        call _gemm('t', 'n',  3, nb, nb,  -zone, wingU, nb, Binv, nb, zzero, bg, 3)
        call _gemm('t', 'n',  3,  3, nb,  -zone, wingU, nb, ag,   nb, zone,  A, 3)
#endif

    end subroutine compute_auxiliary_and_A_2d_general

    !> It computes the auxiliary ag vector and the A tensor
    !> @param[in] Binv     - the inverse of the body
    !> @param[in] wingL    - the lower wing
    !> @param[out] ag      - the auxiliary vector related to lower wing
    !> @param[inout] A     - the A ternsor, on entry contains the head
    !> @param[inout] world - the linear algebra manager
    subroutine compute_auxiliary_and_A_2d_hermitian(Binv, wingL, ag, A, world)
        complex(aip), target, intent(inout)                      :: Binv(:,:)
        complex(aip), target, intent(inout)                      :: wingL(:,:)
        complex(aip), target, allocatable, intent(out)           :: ag(:,:)
        complex(aip), target, intent(inout)                      :: A(:,:)
        type(linalg_world_t), intent(inout)                      :: world

        integer :: nb, i

#ifdef USE_GPU
        type(C_ptr) :: A_dptr, ag_dptr, bg_dptr, Binv_dptr, wingL_dptr
#else
        ! External
        external :: _gemm
#endif

        nb = size(Binv,1)
        allocate(ag(nb,3))

        ! Set arrays to the ones required by the formalism
        do i = 1, 3
            A(i,i) = zone - A(i,i)
        end do

        wingL = -cmplx(1.0_aip/sqrt(fourpi), 0.0_aip, aip) * wingL
        A     = -cmplx(1.0_aip/fourpi, 0.0_aip, aip) * A

#ifdef USE_GPU
        ! Alloc in the device
        call world%register%alloc("A"    , size(A) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("ag"   , size(ag) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("Binv" , size(Binv) * c_sizeof(zzero), world%get_device())
        call world%register%alloc("wingL", size(wingL) * c_sizeof(zzero), world%get_device())

        ! Associate
        call world%register%assoc("A"    , C_loc(A)    )
        call world%register%assoc("ag"   , C_loc(ag)   )
        call world%register%assoc("Binv" , C_loc(Binv) )
        call world%register%assoc("wingL", C_loc(wingL))

        ! Host device movement
        call world%register%to_device("A")
        call world%register%to_device("Binv")
        call world%register%to_device("wingL")

        ! Get device pointers
        A_dptr     = world%register%device_ptr("A"    )
        ag_dptr    = world%register%device_ptr("ag"   )
        Binv_dptr  = world%register%device_ptr("Binv" )
        wingL_dptr = world%register%device_ptr("wingL")

        call _gemm_gpu(MagmaNoTrans, MagmaNoTrans, nb, 3, nb, -zone, Binv_dptr, &
                        nb, wingL_dptr, nb, zzero, ag_dptr, nb, world%get_queue())
        call _gemm_gpu(MagmaConjTrans, MagmaNoTrans, 3, 3, nb, -zone, wingL_dptr, &
                        nb, ag_dptr, nb, zone, A_dptr, 3, world%get_queue())
        call world%syncronize()

        ! Update A
        call world%register%from_device("A")

        ! Deassociate from device but keep them there
        call world%register%deassoc("A")
        call world%register%deassoc("ag")

        ! Delete from device
        call world%register%remove("Binv")
        call world%register%remove("wingL")

#else
        call _gemm('n', 'n', nb, 3, nb, -zone, Binv,  nb, wingL, nb, zzero, ag, nb)
        call _gemm('c', 'n',  3, 3, nb, -zone, wingL, nb, ag,     nb, zone,  A, 3)
#endif

    end subroutine compute_auxiliary_and_A_2d_hermitian


    !> It computes q^T \cdot A \cdot q 
    !> @param[in] A      - A tensor (if GPU it must be loaded)
    !> @param[in] q      - the object q (if GPU it must be loaded)
    !> @param[in] world  - the linalg_world_t handler 
    !> @param[out] qAq   - (if GPU it will also provide the GPU object)
    subroutine compute_qAq(A, q, world, qAq)

        complex(aip), target, allocatable, intent(inout)   :: A(:,:)
        complex(aip), target, allocatable, intent(in)      :: q(:,:)
        type(linalg_world_t), intent(inout)                :: world
        complex(aip), target, allocatable, intent(out)     :: qAq(:)

        complex(aip), target, allocatable :: Aq(:,:)
        integer(i32) :: nr 
        integer(i32) :: i

#ifdef USE_GPU
        type(C_ptr) :: q_dptr, Aq_dptr, A_dptr 
#else
        ! External
        external :: _gemm
#endif 

        nr = size(q,2)
        allocate(Aq, mold=q)

#ifdef USE_GPU
        
        ! Create Aq in the device
        call world%register%alloc("Aq", size(q) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("Aq", C_loc(Aq))

        ! Reassociate q and L to host memory
        call world%register%assoc("q", C_loc(q))
        call world%register%assoc("A", C_loc(A))

        ! Update L in device to the symmetrized one
        call world%register%to_device("A")

        ! Get device pointers
        A_dptr   = world%register%device_ptr("A")
        Aq_dptr  = world%register%device_ptr("Aq")
        q_dptr   = world%register%device_ptr("q")

        call _gemm_gpu(MagmaNoTrans, MagmaNoTrans, 3, nr, 3, zone, A_dptr, &
                         3, q_dptr, 3, zzero, Aq_dptr, 3, world%get_queue())
        call world%syncronize()
#else   
        call _gemm('n', 'n', 3, nr, 3, zone, A, 3, q, 3, zzero, Aq, 3)
#endif  

        allocate(qAq(size(q,2)))
#ifdef USE_GPU 
        call world%register%alloc("qAq", size(qAq) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("qAq", C_loc(qAq))
        !$omp target teams distribute parallel do private(i) 
        do i = 1, size(q,2)
                qAq(i) = dot_product(q(:,i),Aq(:,i))
        end do
        !$omp end target teams distribute parallel do
        
        call world%register%from_device("qAq")
        call world%register%remove("qAq")
        
        call world%register%deassoc("q")
        call world%register%remove("Aq")
        call world%register%remove("A")
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

    !> It computes ag term for the correction of the body
    !> @param[in] q      - the linear algebra object containing q
    !> @param[in] ag     - ag vector
    !> @param[in] world  - the linalg_world_t handler 
    !> @param[out] qag   - (in GPU if compiled for)
    subroutine compute_inverse_ag(q, ag, world, qag)
        complex(aip), target, allocatable,  intent(inout)   :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)   :: ag(:,:)
        type(linalg_world_t),       intent(inout)           :: world
        complex(aip), target, allocatable,  intent(inout)   :: qag(:,:)

        integer(i32) :: nr, nb 

#if USE_GPU
        type(C_ptr) :: qag_dptr, ag_dptr, q_dptr
#else
        ! External
        external :: _gemm
#endif 

        nr = size(q,2)
        nb = size(ag,1)
        allocate(qag(nr,nb))

#ifdef USE_GPU

        ! Create Lq in the device
        call world%register%alloc("qag", size(qag) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("qag", C_loc(qag))

        ! Reassociate q and L to host memory
        call world%register%assoc("q",  C_loc(q))
        call world%register%assoc("ag", C_loc(ag))

        ! Get device pointers
        q_dptr    = world%register%device_ptr("q")
        qag_dptr  = world%register%device_ptr("qag")
        ag_dptr   = world%register%device_ptr("ag")

        call _gemm_gpu(MagmaTrans, MagmaTrans, nr, nb, 3, zone, q_dptr, &
                        3, ag_dptr, nb, zzero, qag_dptr, nr, world%get_queue())
        call world%syncronize()
        
        call world%register%from_device("qag")
        call world%register%remove("ag")
        call world%register%remove("qag")
        call world%register%deassoc("q")
#else   
        call _gemm('T', 'T', nr, nb, 3, zone, q, 3, ag, nb, zzero, qag, nr)
#endif  
        deallocate(ag)

    end subroutine compute_inverse_ag

    !> It computes bg term for the correction of the body
    !> @param[in] q -  qs
    !> @param[in] bg - bg vector
    !> @param[in] world - the linalg_world_t handler 
    !> @param[out] qbg   - (in GPU if compiled for)
    subroutine compute_inverse_bg(q, bg, world, qbg)
        complex(aip), target, allocatable,  intent(inout)  :: q(:,:)
        complex(aip), target, allocatable,  intent(inout)  :: bg(:,:)
        type(linalg_world_t),               intent(inout)  :: world
        complex(aip), target, allocatable,  intent(out)    :: qbg(:,:)

        integer(i32) :: nr, nb 

#if USE_GPU
        type(C_ptr) :: qbg_dptr, bg_dptr, q_dptr
#else
        ! External
        external :: _gemm
#endif 

        nr = size(q,2)
        nb = size(bg,2)
        allocate(qbg(nr,nb))
#ifdef USE_GPU
        call world%register%alloc("qbg", size(qbg) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("qbg", C_loc(qbg))

        ! Reassociate q and L to host memory
        call world%register%assoc("q", C_loc(q))
        call world%register%assoc("bg", C_loc(bg))

        ! Get device pointers
        q_dptr    = world%register%device_ptr("q")
        qbg_dptr  = world%register%device_ptr("qbg")
        bg_dptr   = world%register%device_ptr("bg")

        call _gemm_gpu(MagmaTrans, MagmaNoTrans, nr, nb, 3, zone, q_dptr, &
                        3, bg_dptr, 3, zzero, qbg_dptr, nr, world%get_queue())
        
        call world%syncronize()

        call world%register%from_device("qbg")
        call world%register%remove("bg")
        call world%register%remove("qbg")
        call world%register%deassoc("q")
#else   
        call _gemm('T', 'N', nr, nb, 3, zone, q, 3, T, 3, zzero, qT, nr)
#endif  
        deallocate(bg)

    end subroutine compute_inverse_bg

#undef _getrf
#undef _getri
#undef _gemm
#undef _getrf_gpu
#undef _getri_gpu
#undef _gemm_gpu

end module idiel_linalg
