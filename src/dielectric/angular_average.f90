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
module m_dielectric_average_gamma
    
    use m_constants, only: i64, r64, twopi, pi, zzero, zone, iunit
    use m_crystal_cell, only: cell_t
    use m_crystal_symmetry, only: symmetry_t
    use m_sph_quadrature, only: compute_angular_mesh_lebedev

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
        !> The volume in which the integral is performed
        real(r64), private :: v_integral
        !> Angular mesh (Cartesian)
        real(r64), allocatable, private :: ang(:,:) 
        !> Mesh points (Cartesian)
        real(r64), allocatable, private :: xyz(:,:) 
        !> Weights for the integrals
        real(r64), allocatable, private :: weights(:)
        !> The distance to subcell surface with Gamma in the center
        real(r64), allocatable, private :: kmax(:)
        !> Geometric part of the integral at given point, note that is multiplied by the weight
        real(r64), allocatable, private :: kmax_f(:)
        !> Pointer to the head of the dielectric matrix
        complex(r64), pointer, private :: head(:,:)    => null()
        !> Pointer to the lower wing, i.e. in PW G0
        complex(r64), pointer, private :: wingL(:,:)   => null()
        !> Pointer to the upper wing, i.e. in PW 0G
        complex(r64), pointer, private :: wingU(:,:)   => null()
        !> Pointer to the inverse of the body block of the dielectric matrix 
        complex(r64), pointer, private :: Binv(:,:) => null()
        !> The head of the inverse
        complex(r64) :: inverse_dielectric_head
        !> The lower wing of the inverse
        complex(r64), allocatable :: inverse_dielectric_wingL(:)
        !> The upper wing of the inverse
        complex(r64), allocatable :: inverse_dielectric_wingU(:)
        !> The body of the inverse
        complex(r64), allocatable :: inverse_dielectric_body(:,:)
        !> Kind of quadrature used
        character(len=:), allocatable, private :: quadrature_scheme
        !> Number of points of the quadrature
        integer(i64), private :: quadrature_npoints 
    contains
        procedure, public  :: init_common, set_dielectric_blocks, compute_anisotropic_avg
        final :: clean
    end type inverse_dielectric_t

contains

    !> This initializes everything common (i.e. without frequency dependency) to all computations
    !> @param[in,out] this - inverse_dielectric_t object
    !> @param[in]     lattice  - lattice vectors given in rows [nm]
    !> @param[in]     redpos   - reduced positions (3,natoms)
    !> @param[in]     elements - list of elements
    !> @param[in]     nq       - the BZ mesh
    !> @param[]
    subroutine init_common(this, lattice, redpos, elements, nq, report)

        class(inverse_dielectric_t), intent(inout) :: this
        real(r64), intent(in) :: lattice(3,3)
        real(r64), allocatable, intent(in) :: redpos(:,:) 
        integer(r64), allocatable, intent(in) :: elements(:)
        integer(i64), intent(in) :: nq(3)
        integer, optional, intent(in) :: report

        ! Locals
        integer(i64) :: ir, il, ii
        real(r64) :: v_bz
        real(r64) :: rel_error
        real(r64), parameter :: onethird = 1.0_r64 / 3.0_r64

        ! Initialize the crystal structure
        call this%cell%initialize(lattice, redpos, elements)

        ! Init reciprocal mesh size
        this%nq = nq

        ! Initialize the symmetry
        call this%symmetry%initialize(this%cell)

        ! Compute the volume in which the integral will be performed
        v_bz = 8.0 * pi**3 / this%cell%vuc 
        this%v_integral = product(this%nq) / v_bz

        ! Initalize the angular mesh and the weights for the angular part of the integral
        call compute_angular_mesh_lebedev(this%ang, this%weights, this%xyz)

        ! Get the number of points used in the angular quadrature
        this%quadrature_npoints = size(this%xyz, 1)

        ! Compute kmax (the boundary)
        call this%cell%get_kmax_subcell_bz(this%nq, this%xyz, this%kmax)

        ! Construct the geometric function
        allocate(this%kmax_f(this%quadrature_npoints))

        ! K function :  \frac{q_{max}^{3}}{3V_{\Gamma}} 
        this%kmax_f(:) =  this%weights * this%v_integral * onethird * this%kmax(:)**3
        
        if (present(report)) then
            rel_error = abs(sum(this%kmax_f) - 1.0_r64)
            write(report,*) 'inverse_dielectric_t%init_common: ' 
            write(report,*) '  Relative error of the volume integral : ', rel_error
        end if

    end subroutine init_common

    !> This nullify and deallocates the objects
    !> @param[in] this - inverse_dielectric_t object
    subroutine clean(this)

        type(inverse_dielectric_t), intent(inout) :: this

        if (associated(this%head))    nullify(this%head)
        if (associated(this%wingL))   nullify(this%wingL)
        if (associated(this%wingU))   nullify(this%wingU)
        if (associated(this%Binv)) nullify(this%Binv)
        if (allocated(this%weights)) deallocate(this%weights)
        if (allocated(this%xyz)) deallocate(this%xyz)
        if (allocated(this%ang)) deallocate(this%ang)
        if (allocated(this%kmax)) deallocate(this%kmax)
        if (allocated(this%kmax_f)) deallocate(this%kmax_f)
        if (allocated(this%inverse_dielectric_wingL)) deallocate(this%inverse_dielectric_wingL)
        if (allocated(this%inverse_dielectric_wingU)) deallocate(this%inverse_dielectric_wingU)
        if (allocated(this%inverse_dielectric_body))  deallocate(this%inverse_dielectric_body)

    end subroutine clean

    !> This subroutine inits pointers to a given value
    !> @param[in] this - inverse_dielectric_t object
    !> @param[in] h    - head of the dielectric matrix
    !> @param[in] wl   - lower wing of the dielectric matrix
    !> @param[in] wu   - upper wing of the dielectric matrix
    !> @param[in] ib   - inverse of the body of the dielectric matrix
    subroutine set_dielectric_blocks(this, h, wl, wu, ib)
        
        class(inverse_dielectric_t), intent(inout) :: this
        complex(r64), target, intent(in)      :: h(:,:)
        complex(r64), target, intent(in)      :: wl(:,:)
        complex(r64), target, intent(in)      :: wu(:,:)
        complex(r64), target, intent(in)      :: ib(:,:)

        ! Associating to internal objects
        this%head    =>  h
        this%wingL   =>  wl
        this%wingU   =>  wu
        this%Binv    =>  ib

    end subroutine set_dielectric_blocks
    
    !> This computes the block inverse at Gamma
    !> @param[in] this   - the current inverse_dielectric_t object for which to compute the average
    subroutine compute_anisotropic_avg(this)

        class(inverse_dielectric_t), intent(inout) :: this

        ! Auxiliary vectors
        complex(r64), allocatable :: U(:,:), V(:,:)
        complex(r64), allocatable :: S(:,:), T(:,:)
        
        ! local field effects 
        complex(r64) :: lfe(3,3)

        ! Macroscopic dielectric matrix
        complex(r64) :: L(3,3)

        ! Function from which to compute the integral 
        complex(r64), allocatable :: kmax_f(:)
        complex(r64), allocatable :: head_f(:) 
        complex(r64), allocatable :: wingL_f(:,:)
        complex(r64), allocatable :: wingU_f(:,:)
        complex(r64), allocatable :: body_f(:,:,:)

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
        allocate(wingU_f(this%quadrature_npoints, nbasis), source=zzero)
        allocate(body_f(this%quadrature_npoints, nbasis, nbasis), source=zzero)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           AUXILIARY VECTORS         ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Init the U_{G\alpha} auxiliary vector (Eq. B.7) in 10.1016/j.cpc.2006.07.018
        allocate(U, source=this%wingL)

        ! Compute S_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        S = matmul(this%Binv,U)

        ! Now for the upper wing
        allocate(V, source=transpose(this%wingU))
        T = matmul(V, this%Binv)

        ! Clean
        deallocate(U,V)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !      Compute macroscopic dielectric matrix      ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute the macroscopic dielectric matrix (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        ! First compute the local field effects
        lfe = matmul(transpose(this%wingU),S)
        ! Now compute the macroscopic dielectric matrix
        L   = this%head - lfe 
        ! Symmetrize the elements
        L   = this%symmetry%symmetryze_complex_tensor(L)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               Compute the functions                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! HEAD: Compute \frac{1.0}{\mathbf{\hat{q} L \hat{q}}}  for the space grid
        head_f(:)  =  this%kmax_f(:) / sum(transpose(this%xyz) * matmul(L,transpose(this%xyz)), dim=1)

        ! WINGS : Compute \frac{ \mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} and
        ! \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} 
        ! Note that there is head missing term to allow for the efficient computation of the body term
        wingL_f =  matmul(this%xyz, transpose(S))
        wingU_f =  matmul(this%xyz, T)
        
        ! Compute the part the body
        ! \frac{[\mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})] [\mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})]}{\mathbf{\hat{q} L \hat{q}}} 
        !$omp parallel shared(this, wingL_f, wingU_f, body_f) private(ii, jj)
        !$omp do 
        do ii = 1, nbasis
            do jj = 1, nbasis
                body_f(:, ii, jj) = head_f(:) * wingL_f(:, ii) * wingU_f(:, jj)
            end do 
        end do
        !$omp end do
        !$omp end parallel

        ! Here finish the Wings object and build the body 
        !$omp parallel shared(wingL_f, wingU_f, head_f) private(ii)
        !$omp do 
        do ii = 1, nbasis
            wingL_f(:, ii) =  - head_f(:) * wingL_f(:, ii)
            wingU_f(:, ii) =  - head_f(:) * wingU_f(:, ii)
        end do
        !$omp end do
        !$omp end parallel

        deallocate(S,T)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Do the anisotropic averaging !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(this%inverse_dielectric_wingL(nbasis), source=zzero)
        allocate(this%inverse_dielectric_wingU(nbasis), source=zzero)
        allocate(this%inverse_dielectric_body, source = this%Binv)

        this%inverse_dielectric_head  = sum(head_f)
        
        !$omp parallel shared(this, wingL_f, wingU_f, body_f) private(ii, jj)
        !$omp do 
        do ii = 1, nbasis
            this%inverse_dielectric_wingL(ii) = sum(wingL_f(:,ii))
            this%inverse_dielectric_wingU(ii) = sum(wingU_f(:,ii))
            do jj = 1, nbasis
                this%inverse_dielectric_body(ii,jj) = this%inverse_dielectric_body(ii,jj) + &
                    sum(body_f(:,ii,jj) )
            end do
        end do 
        !$omp end do
        !$omp end parallel

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(head_f, wingL_f, wingU_f, body_f)
    
    end subroutine compute_anisotropic_avg

end module m_dielectric_average_gamma
