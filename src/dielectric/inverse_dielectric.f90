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
!> Contains procedures for averages of the dielectric matrix

!> This module contains the procedures for the computation of the dielectric matrix averages
module idiel_inverse_dielectric
    
    use idiel_constants, only: i64, r64, twopi, pi, zzero, zone, iunit
    use idiel_crystal_cell, only: cell_t
    use idiel_crystal_symmetry, only: symmetry_t
    use idiel_sph_quadrature, only: compute_angular_mesh_lebedev_131, compute_angular_mesh_lebedev_21
    use idiel_spherical_harmonics, only: sph_harm, sph_harm_expansion
#ifdef USE_GPU
    use idiel_gpu_magma_t, only: linalg_world_t, linalg_obj_t
#else
    use idiel_cpu_magma_t, only: linalg_world_t, linalg_obj_t
#endif

    implicit none

    private
    public inverse_dielectric_t

    !> This type is used to perform all the operations of averaging
    type inverse_dielectric_t
        !> Cell information
        type(cell_t), private :: cell
        !> Symmetry information
        type(symmetry_t), private :: symmetry
        !> The BZ mesh size
        integer(i64), private :: nq(3) 
        !> Angular mesh (Cartesian)
        real(r64), allocatable, private :: ang(:,:) 
        !> Mesh points (Cartesian)
        complex(r64), allocatable, private :: xyz(:,:) 
        !> Linear algebra handler of xyz
        type(linalg_obj_t), private :: ref_xyz
        !> Weights for the integrals
        real(r64), allocatable, private :: weights(:)
        !> Spherical harmonics
        complex(r64), allocatable, private :: ylm(:,:)
        !> angular integrals]
        complex(r64), allocatable, private :: angular_integrals(:)
        !> Pointer to the head of the dielectric matrix
        complex(r64), pointer, private :: head(:,:)    => null()
        !> Pointer to the lower wing, i.e. in PW G0
        complex(r64), pointer, private :: wingL(:,:)   => null()
        !> Pointer to the upper wing, i.e. in PW 0G
        complex(r64), pointer, private :: wingU(:,:)   => null()
        !> Pointer to the inverse of the body block of the dielectric matrix 
        complex(r64), pointer, private :: Binv(:,:)    => null()
        !> Actual data in case that it has been computed using this library
        complex(r64), allocatable :: Binv_data(:,:)
        !> The head of the inverse
        complex(r64) :: inverse_dielectric_head
        !> The lower wing of the inverse
        complex(r64), allocatable :: inverse_dielectric_wingL(:)
        !> The upper wing of the inverse
        complex(r64), allocatable :: inverse_dielectric_wingU(:)
        !> The body of the inverse
        complex(r64), allocatable :: inverse_dielectric_body(:,:)
        !> Number of points of the quadrature
        integer(i64), private :: quadrature_npoints 
        !> The handler of linear algebra queues
        type(linalg_world_t), private :: world
    contains
        procedure, public  :: init_common, set_dielectric_blocks, invert_body, get_n_basis, compute_anisotropic_avg
        procedure, private :: compute_anisotropic_avg_hermitian, compute_anisotropic_avg_general
        final :: clean
    end type inverse_dielectric_t

    !> Order of harmonic spherics expansion
    integer(i64), parameter :: lmax = 10_i64
    ! Number of spherical harmonics
    integer(i64), parameter :: nsph = (lmax + 1)**2

contains

    !> This initializes everything common (i.e. without frequency dependency) to all computations
    !> @param[in,out] this - inverse_dielectric_t object
    !> @param[in]     lattice  - lattice vectors given in rows [nm]
    !> @param[in]     redpos   - reduced positions (3,natoms)
    !> @param[in]     elements - list of elements
    !> @param[in]     nq       - the BZ mesh
    subroutine init_common(this, lattice, redpos, elements, nq)

        class(inverse_dielectric_t), intent(inout) :: this
        real(r64), intent(in)         :: lattice(3,3)
        real(r64),  intent(in)        :: redpos(:,:) 
        integer(r64),  intent(in)     :: elements(:)
        integer(i64), intent(in)      :: nq(3)

        ! Locals
        integer(i64) :: ii
        real(r64) :: v_bz
        real(r64) :: rel_error
        real(r64), parameter :: onethird = 1.0_r64 / 3.0_r64
        real(r64), allocatable :: xyz(:,:)

        ! The volume in which the integral is performed
        real(r64) :: v_integral
        ! The distance to subcell surface with Gamma in the center
        real(r64), allocatable :: kmax(:)
        ! Geometric part of the integral at given point, note that is multiplied by the weight
        real(r64), allocatable :: kmax_f(:)
        ! Spherical harmonics
        complex(r64), allocatable :: ylm(:,:)

        ! Initialize the crystal structure
        call this%cell%initialize(lattice, redpos, elements)

        ! Init reciprocal mesh size
        this%nq = nq

        ! Initialize the symmetry
        call this%symmetry%initialize(this%cell)

        ! Compute the volume in which the integral will be performed
        v_bz = 8.0 * pi**3 / this%cell%vuc 
        v_integral = product(this%nq) / v_bz

        ! Initalize the big angular mesh and the weights for the angular part of the integral
        call compute_angular_mesh_lebedev_131(this%ang, this%weights, xyz)

        ! Get the number of points used in the angular quadrature
        this%quadrature_npoints = size(xyz, 1)

        ! Compute kmax (the boundary)
        call this%cell%get_kmax_subcell_bz(this%nq, xyz, kmax)

        ! Construct the geometric function
        allocate(kmax_f(this%quadrature_npoints))

        ! Compute the geometric part of the integral
        ! Notice that this needs to be done in a dense mesh since
        ! we are trying to integrate parallelepiped using spherical coordinates 
        ! and that is an object of infinite degree as for its spherical harmonic expansion
        ! K function :  \frac{q_{max}^{3}}{3V_{\Gamma}} 
        kmax_f(:) =  this%weights * v_integral * onethird * kmax(:)**3

        ! Precompute the angular part, we first need the spherical harmonics in the big mesh
        call sph_harm(lmax, this%ang, ylm)

        allocate(this%angular_integrals(nsph))

        !$omp parallel shared(ylm, this, kmax_f) private(ii)
        !$omp do schedule(dynamic)
        do ii = 1, nsph
            this%angular_integrals(ii) = sum(ylm(:, ii) * kmax_f(:))
        end do 
        !$omp end do
        !$omp end parallel

        ! Set now to smaller mesh, as dielectric terms converge quite fast with l in comparison
        ! to the geometric part, and thus a smaller Lebedev order can be used
        deallocate(this%ang, this%weights, xyz)
        ! Recompute things in the small mesh
        call compute_angular_mesh_lebedev_21(this%ang, this%weights, xyz)
        allocate(this%xyz,source=transpose(cmplx(xyz,0.0,r64)))
        this%quadrature_npoints = size(xyz, 1)

        ! Compute the spherical harmonics in the smaller mesh
        call sph_harm(lmax, this%ang, this%ylm)

        ! Init algebra world
        call this%world%init()
        call this%ref_xyz%allocate_gpu(this%xyz)
        call this%ref_xyz%transfer_cpu_gpu(this%xyz, this%world)

    end subroutine init_common

    !> This nullify and deallocates the objects
    !> @param[in] this - inverse_dielectric_t object
    subroutine clean(this)

        type(inverse_dielectric_t), intent(inout) :: this

        if (associated(this%head))    nullify(this%head)
        if (associated(this%wingL))   nullify(this%wingL)
        if (associated(this%wingU))   nullify(this%wingU)
        if (associated(this%Binv)) nullify(this%Binv)
        if (allocated(this%Binv_data)) deallocate(this%Binv_data)
        if (allocated(this%weights)) deallocate(this%weights)
        if (allocated(this%angular_integrals)) deallocate(this%angular_integrals)
        if (allocated(this%xyz)) deallocate(this%xyz)
        if (allocated(this%ang)) deallocate(this%ang)
        if (allocated(this%ylm)) deallocate(this%ylm)
        if (allocated(this%inverse_dielectric_wingL)) deallocate(this%inverse_dielectric_wingL)
        if (allocated(this%inverse_dielectric_wingU)) deallocate(this%inverse_dielectric_wingU)
        if (allocated(this%inverse_dielectric_body))  deallocate(this%inverse_dielectric_body)
        call this%ref_xyz%destroy()
        if (this%world%is_queue_set()) call this%world%finish()

    end subroutine clean

    !> This subroutine inits pointers to a given value
    !> @param[in] this - inverse_dielectric_t object
    !> @param[in] h    - head of the dielectric matrix
    !> @param[in] wl   - lower wing of the dielectric matrix
    !> @param[in] wu   - upper wing of the dielectric matrix
    !> @param[in] ib   - inverse of the body of the dielectric matrix (optional)
    subroutine set_dielectric_blocks(this, h, wl, wu, ib)
        
        class(inverse_dielectric_t), intent(inout) :: this
        complex(r64), target, intent(in)      :: h(:,:)
        complex(r64), target, intent(in)      :: wl(:,:)
        complex(r64), target, intent(in)      :: wu(:,:)
        complex(r64), target, optional, intent(in) :: ib(:,:)

        ! Associating to internal objects
        this%head    =>  h
        this%wingL   =>  wl
        this%wingU   =>  wu
        if (present(ib)) this%Binv    =>  ib

    end subroutine set_dielectric_blocks

    !> This computes the block inverse at Gamma 
    !> @param[in] this       - the current inverse_dielectric_t object for which to compute the average
    !> @param[in] hermitian  - is the dielectric matrix hermitian
    subroutine compute_anisotropic_avg(this, hermitian)
        
        class(inverse_dielectric_t), intent(inout) :: this
        logical, intent(in)                        :: hermitian

        if (hermitian) then
            call this%compute_anisotropic_avg_hermitian()
        else
            call this%compute_anisotropic_avg_general()
        end if

    end subroutine compute_anisotropic_avg
    
    !> This computes the block inverse at Gamma when the dielectric matrix is Hermitic
    !> @param[in] this       - the current inverse_dielectric_t object for which to compute the average
    subroutine compute_anisotropic_avg_hermitian(this)

        use idiel_linalg

        class(inverse_dielectric_t), intent(inout) :: this

        ! Auxiliary vectors
        complex(r64), allocatable :: S(:,:)
        type(linalg_obj_t) :: ref_S
        
        ! Macroscopic dielectric matrix
        complex(r64), allocatable :: L(:,:)
        type(linalg_obj_t) :: ref_L

        ! Function from which to compute the integral 
        complex(r64), allocatable :: head_f(:) 
        complex(r64), allocatable :: wingL_f(:,:)
        complex(r64), allocatable :: body_f(:)

        ! Harmonic expansion coefficients
        complex(r64) :: clm_head(nsph)
        complex(r64) :: clm_body(nsph) 
        real(r64)    :: error 

        ! Basis size
        integer(i64) :: nbasis
        
        ! Dummy indexes
        integer(i64) :: ii, jj

        ! Get the basis size
        nbasis = size(this%Binv, 1)

        ! Free final containers
        if (allocated(this%inverse_dielectric_wingL)) deallocate(this%inverse_dielectric_wingL)
        if (allocated(this%inverse_dielectric_wingU)) deallocate(this%inverse_dielectric_wingU)
        if (allocated(this%inverse_dielectric_body))  deallocate(this%inverse_dielectric_body)

        ! Allocate space
        allocate(head_f(this%quadrature_npoints))
        allocate(wingL_f(this%quadrature_npoints, nbasis))
        allocate(body_f(this%quadrature_npoints))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectorws and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute S_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(L, source=this%head)
        call compute_auxiliary_and_macroscopic(this%Binv, this%wingL, S, ref_S, L, ref_L, this%world)
        
        ! Symmetrize the elements of the macroscopic dielectric matrix and update it
        L   = this%symmetry%symmetryze_complex_tensor(L) 
        call ref_L%transfer_cpu_gpu(L, this%world)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               Compute the functions                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! HEAD: Compute \frac{\omega(\Omega)}{\mathbf{\hat{q} L \hat{q}}}  for the space grid
        head_f(:)  =  compute_inverse_head(ref_L, this%xyz, this%ref_xyz, this%world)
        call ref_L%destroy()
        deallocate(L)

        ! WINGS : Compute \frac{ \mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} and
        ! \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} 
        ! Note that there is head missing term to allow for the efficient computation of the body term
        wingL_f =  compute_inverse_wingL(this%ref_xyz, ref_S, this%world)
        call ref_S%destroy()
        deallocate(S)
        
        ! The body is directly computed as saving it to RAM is too intensive
        ! Note that wings functions are odd functions and thus their integral is zero, so no more work is done regarding those.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Do the anisotropic averaging !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        allocate(this%inverse_dielectric_wingL(nbasis), source=zzero)
        allocate(this%inverse_dielectric_wingU(nbasis), source=zzero)
        allocate(this%inverse_dielectric_body, source = this%Binv)
        
        call sph_harm_expansion(lmax, head_f, this%weights, this%ylm, clm_head)
        this%inverse_dielectric_head  = sum(clm_head(:) * this%angular_integrals(:))

        error = maxval(abs(clm_head(100:121)))
        if (error > 1.0e-8_r64) then
            write(*,*) "Warning (compute_anisotropic_avg_hermitian) the expansion coefficient can be not large enough"
        end if
        
        ! Here we compute the body 
        ! \frac{[\mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})] [\mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})]}{\mathbf{\hat{q} L \hat{q}}} 
        !$omp parallel shared(this, head_f, wingL_f, nbasis) private(ii, jj, body_f)
        !$omp do schedule(dynamic) collapse(2)
        do ii = 1, nbasis
            do jj = 1, nbasis
                body_f(:) = head_f(:) * wingL_f(:, jj) * conjg(wingL_f(:, ii))
                call sph_harm_expansion(lmax, body_f, this%weights, this%ylm, clm_body)
                this%inverse_dielectric_body(jj,ii) = this%inverse_dielectric_body(jj,ii) + &
                    sum(clm_body(:) * this%angular_integrals(:))
            end do
        end do 
        !$omp end do
        !$omp end parallel

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(head_f, wingL_f, body_f)
    
    end subroutine compute_anisotropic_avg_hermitian


        !> This computes the block inverse at Gamma for a general dielectric matrix 
    !> @param[in] this       - the current inverse_dielectric_t object for which to compute the average
    subroutine compute_anisotropic_avg_general(this)

        use idiel_linalg

        class(inverse_dielectric_t), intent(inout) :: this

        ! Auxiliary vectors
        complex(r64), allocatable :: S(:,:), T(:,:)
        type(linalg_obj_t) :: ref_S, ref_T
        
        ! Macroscopic dielectric matrix
        complex(r64), allocatable :: L(:,:)
        type(linalg_obj_t) :: ref_L

        ! Function from which to compute the integral 
        complex(r64), allocatable :: head_f(:) 
        complex(r64), allocatable :: wingL_f(:,:)
        complex(r64), allocatable :: wingU_f(:,:)
        complex(r64), allocatable :: body_f(:)

        ! Harmonic expansion coefficients
        complex(r64) :: clm_head(nsph)
        complex(r64) :: clm_body(nsph) 
        real(r64) :: error

        ! Basis size
        integer(i64) :: nbasis
        
        ! Dummy indexes
        integer(i64) :: ii, jj

        ! Get the basis size
        nbasis = size(this%Binv, 1)

        ! Free final containers
        if (allocated(this%inverse_dielectric_wingL)) deallocate(this%inverse_dielectric_wingL)
        if (allocated(this%inverse_dielectric_wingU)) deallocate(this%inverse_dielectric_wingU)
        if (allocated(this%inverse_dielectric_body))  deallocate(this%inverse_dielectric_body)

        ! Allocate space
        allocate(head_f(this%quadrature_npoints))
        allocate(wingL_f(this%quadrature_npoints, nbasis))
        allocate(body_f(this%quadrature_npoints))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectorws and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute S_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! it also computes in this case the term for the upper wing
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(L, source=this%head)
        call compute_auxiliary_and_macroscopic(this%Binv, this%wingL, S, ref_S, L, ref_L, this%world, this%wingU, T, ref_T)
        
        ! Symmetrize the elements of the macroscopic dielectric matrix and update it
        L   = this%symmetry%symmetryze_complex_tensor(L) 
        call ref_L%transfer_cpu_gpu(L, this%world)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               Compute the functions                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! HEAD: Compute \frac{\omega(\Omega)}{\mathbf{\hat{q} L \hat{q}}}  for the space grid
        head_f(:)  =  compute_inverse_head(ref_L, this%xyz, this%ref_xyz, this%world)
        call ref_L%destroy()
        deallocate(L)

        ! WINGS : Compute \frac{ \mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} and
        ! \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} 
        ! Note that there is head missing term to allow for the efficient computation of the body term
        wingL_f =  compute_inverse_wingL(this%ref_xyz, ref_S, this%world)
        call ref_S%destroy()
        deallocate(S)
        wingU_f = compute_inverse_wingU(this%ref_xyz, ref_T, this%world)
        call ref_T%destroy()
        deallocate(T)
        
        ! The body is directly computed as saving it to RAM is too intensive

        ! Note that wings functions are odd functions and thus their integral is zero, so no more is done regarding those.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Do the anisotropic averaging !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        allocate(this%inverse_dielectric_wingL(nbasis), source=zzero)
        allocate(this%inverse_dielectric_wingU(nbasis), source=zzero)
        allocate(this%inverse_dielectric_body, source = this%Binv)
        
        call sph_harm_expansion(lmax, head_f, this%weights, this%ylm, clm_head)
        this%inverse_dielectric_head  = sum(clm_head(:) * this%angular_integrals(:))

        error = maxval(abs(clm_head(100:121)))
        if (error > 1.0e-8_r64) then
            write(*,*) "Warning (compute_anisotropic_avg_general) the expansion coefficient can be not large enough"
        end if
        
        ! Here we compute the body 
        ! \frac{[\mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})] [\mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})]}{\mathbf{\hat{q} L \hat{q}}} 
        !$omp parallel shared(this, head_f, wingL_f, wingU_f, nbasis) private(ii, jj, body_f)
        !$omp do schedule(dynamic) collapse(2)
        do ii = 1, nbasis
            do jj = 1, nbasis
                body_f(:) = head_f(:) * wingL_f(:, jj) * wingU_f(:, ii)
                call sph_harm_expansion(lmax, body_f, this%weights, this%ylm, clm_body)
                this%inverse_dielectric_body(jj,ii) = this%inverse_dielectric_body(jj,ii) + &
                    sum(clm_body(:) * this%angular_integrals(:))
            end do
        end do 
        !$omp end do
        !$omp end parallel

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(head_f, wingL_f, wingU_f, body_f)
    
    end subroutine compute_anisotropic_avg_general

    !> Inverts the body and stores it (GPU or CPU depending on the compilation)
    !> @param[in] this - the inverse_dielectric_t in which to store the inverse of body
    !> @param[in] body - the body to invert
    subroutine invert_body(this, body)
    
        use idiel_linalg, only: inverse_complex_LU

        class(inverse_dielectric_t), intent(inout), target :: this
        complex(r64), intent(in) :: body(:,:)

        if (.not. this%world%is_queue_set()) call this%world%init()

        call inverse_complex_LU(body, this%Binv_data, this%world)        
        this%Binv => this%Binv_data

    end subroutine invert_body

    !> This function returns the basis size (for CXX and Python compatibility)
    !> @param[in] this - the inverse_dielectric_t object to check the size of the used basis
    function get_n_basis(this) result(nbasis)
        class(inverse_dielectric_t), intent(inout), target :: this
        integer(i64) :: nbasis
        if (.not. associated(this%Binv)) error stop "inverse_dielectric_t%get_n_basis: Error set inverse dielectric matrix for this" 
        nbasis = size(this%Binv, 1)
    end function

end module idiel_inverse_dielectric
