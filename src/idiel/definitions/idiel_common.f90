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
submodule (idiel) idiel_common

    implicit none

contains

    module subroutine init_common(this, lattice, redpos, elements, nq, dim)

        class(idiel_t), intent(inout) :: this
        real(r64), intent(in)         :: lattice(3,3)
        real(r64),  intent(in)        :: redpos(:,:) 
        integer(r64),  intent(in)     :: elements(:)
        integer(i64), intent(in)      :: nq(3)
        integer(i64), intent(in), optional :: dim

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
        ! Some cutoff values to create the radial mesh for the 2D case
        real(r64) :: xlim, xcut, rmax, dx1, dx2

        ! Initialize the crystal structure
        call this%cell%initialize(lattice, redpos, elements)

        ! Init reciprocal mesh size
        this%nq = nq

        ! Initialize the symmetry
        call this%symmetry%initialize(this%cell)

        ! Process dimensionality
        if (present(dim)) this%dim = dim

        select case(this%dim)
        
        case(2) ! 2D case

            ! rcut
            this%rcut = 0.5_r64 * this%cell%lattice(3,3) 

            ! Initalize the big angular mesh and the weights for the angular part of the integral
            call compute_angular_mesh_gauss_legendre(size_mesh_2d_fine , this%ang, this%weights_fine, xyz)

            ! Add the area factor to the weights 
            this%weights_fine = this%weights_fine / this%cell%area_r_12

            ! Save the large phi mesh
            allocate(this%phi, source=this%ang(:,2))
            
            ! Compute circular basis in the fine mesh
            call circ_harm(lmax, this%ang(:,2), this%blm_fine)

            ! Get the number of points used in the angular quadrature
            this%quadrature_npoints = size(xyz, 1)

            ! Compute kmax (the boundary)
            call this%cell%get_kmax_subcell_bz(this%nq, xyz, this%rmax2d)

            ! Init the radial mesh (this mesh is intended to provide good treatment of the queue)
            rmax = maxval(this%rmax2d)
            xlim = 0.005_r64
            xcut = -log(xlim) / this%rcut
            if (xcut <= 0.5_r64 * rmax ) then  
                dx1  = (1.0_r64 - exp(-xcut)) / (nr/2 - 1)
                dx2  = (rmax - xcut) / (nr / 2)
                do ii = 1, nr / 2
                    this%radii(ii) = -log(-(ii - 1) * dx1 + 1)
                    this%radii(ii + nr/2) = xcut + ii * dx2
                end do 
            else
                xcut = rmax
                dx1  = (1.0_r64 - exp(-xcut)) / (nr - 1)
                do ii = 1, nr
                    this%radii(ii) = -log(-(ii - 1) * dx1 + 1)
                end do
            end if 

            do ii = 1, nr
                this%vr(:,ii) = fourpi * (1.0_r64 - exp(-this%rcut * this%radii(ii)))
            end do

            ! Compute the small mesh
            deallocate(this%ang, xyz)
            call compute_angular_mesh_gauss_legendre(size_mesh_2d_coarse, this%ang, this%weights, xyz)
            allocate(this%xyz, source=transpose(cmplx(xyz,0.0,r64)))
            
            ! Compute circular basis in the coarse mesh
            call circ_harm(lmax, this%ang(:,2), this%blm_coarse)

        case(3) ! 3D case
            ! Compute the volume in which the integral will be performed
            v_bz = 8.0 * pi**3 / this%cell%vuc 
            v_integral = product(this%nq) / v_bz

            ! Initalize the big angular mesh and the weights for the angular part of the integral
            call compute_angular_mesh_lebedev_131(this%ang, this%weights_fine, xyz)

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
            kmax_f(:) =  this%weights_fine * v_integral * onethird * kmax(:)**3

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
            deallocate(this%ang, this%weights_fine, xyz)
            ! Recompute things in the small mesh
            call compute_angular_mesh_lebedev_21(this%ang, this%weights, xyz)
            allocate(this%xyz,source=transpose(cmplx(xyz,0.0,r64)))
            this%quadrature_npoints = size(xyz, 1)

            ! Compute the spherical harmonics in the smaller mesh
            call sph_harm(lmax, this%ang, this%ylm)
        case default
            error stop "Error(idiel_t%init_common): Error dimension should be either 3 or 2"
        end select 

        ! Init algebra world
        call this%world%init()
        call this%ref_xyz%allocate_gpu(this%xyz)
        call this%ref_xyz%transfer_cpu_gpu(this%xyz, this%world)

    end subroutine init_common

    !> This nullify and deallocates the objects
    !> @param[in] this - idiel_t object
    module subroutine clean(this)

        type(idiel_t), intent(inout) :: this

        if (associated(this%head))    nullify(this%head)
        if (associated(this%wingL))   nullify(this%wingL)
        if (associated(this%wingU))   nullify(this%wingU)
        if (associated(this%Binv)) nullify(this%Binv)
        if (allocated(this%Binv_data)) deallocate(this%Binv_data)
        if (allocated(this%weights_fine)) deallocate(this%weights_fine)
        if (allocated(this%weights)) deallocate(this%weights)
        if (allocated(this%xyz)) deallocate(this%xyz)
        if (allocated(this%ang)) deallocate(this%ang)
        if (allocated(this%idiel_wingL)) deallocate(this%idiel_wingL)
        if (allocated(this%idiel_wingU)) deallocate(this%idiel_wingU)
        if (allocated(this%idiel_body))  deallocate(this%idiel_body)
        call this%ref_xyz%destroy()
        if (this%world%is_queue_set()) call this%world%finish()

        ! 2D stuff
        if (allocated(this%phi)) deallocate(this%phi)
        if (allocated(this%rmax2d))  deallocate(this%rmax2d)
        if (allocated(this%blm_coarse))  deallocate(this%blm_coarse)
        if (allocated(this%blm_fine))  deallocate(this%blm_fine)
        
        ! 3D stuff
        if (allocated(this%angular_integrals)) deallocate(this%angular_integrals)
        if (allocated(this%ylm)) deallocate(this%ylm)
        
    end subroutine clean

    !> This subroutine inits pointers to a given value
    !> @param[in] this - idiel_t object
    !> @param[in] h    - head of the dielectric matrix
    !> @param[in] wl   - lower wing of the dielectric matrix
    !> @param[in] wu   - upper wing of the dielectric matrix
    !> @param[in] ib   - inverse of the body of the dielectric matrix (optional)
    module subroutine set_dielectric_blocks(this, h, wl, wu, ib)
        
        class(idiel_t), intent(inout) :: this
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


    module subroutine invert_body(this, body)
        
        use idiel_linalg, only: inverse_complex_LU

        class(idiel_t), intent(inout), target :: this
        complex(r64), intent(in) :: body(:,:)

        if (.not. this%world%is_queue_set()) call this%world%init()

        call inverse_complex_LU(body, this%Binv_data, this%world)        
        this%Binv => this%Binv_data

    end subroutine invert_body

    module function get_n_basis(this) result(nbasis)
        class(idiel_t), intent(inout), target :: this
        integer(i64) :: nbasis
        if (.not. associated(this%Binv)) error stop "idiel_t%get_n_basis: Error set inverse dielectric matrix for this" 
        nbasis = size(this%Binv, 1)
    end function

end submodule idiel_common
