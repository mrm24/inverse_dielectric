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
!> Contains procedures for averages of the dielectric matrix

!> This module contains the procedures for the computation of the dielectric matrix averages
submodule (idiel) idiel_2D

    implicit none

contains

    module subroutine compute_anisotropic_avg_scrcoulomb_2d(this, hermitian)
            
        class(idiel_t), intent(inout) :: this
        logical, intent(in)           :: hermitian

        if (this%dim /= 2) then
            error stop "Error(idiel_t%compute_anisotropic_avg_scrcoulomb_2d): this procedure can only be used with 2D materials"
        end if

        if (hermitian) then
            call this%compute_anisotropic_avg_scrcoulomb_2d_hermitian()
        else
            call this%compute_anisotropic_avg_scrcoulomb_2d_general()
        end if

    end subroutine compute_anisotropic_avg_scrcoulomb_2d
    
    module subroutine compute_anisotropic_avg_scrcoulomb_2d_general(this)

        use idiel_linalg

        class(idiel_t), target, intent(inout) :: this

        ! Auxiliary vectors
        complex(aip), allocatable :: ag(:,:), bg(:,:)
        
        ! Macroscopic dielectric matrix
        complex(aip), allocatable :: A(:,:)

        ! Function from which to compute the integral 
        complex(aip), target, allocatable :: qAq(:) 
        complex(aip), target, allocatable :: wingL_f(:,:)
        complex(aip), target, allocatable :: wingU_f(:,:)
        complex(aip), target, allocatable :: wLwU_f(:)

        ! Basis size
        integer(i32) :: nbasis
        
        ! Dummy indexes
        integer(i32) :: ii, jj

        ! Get the basis size
        nbasis = size(this%Binv, 1)

        ! Free final containers
        if (allocated(this%idiel_wingL)) deallocate(this%idiel_wingL)
        if (allocated(this%idiel_wingU)) deallocate(this%idiel_wingU)
        if (allocated(this%idiel_body))  deallocate(this%idiel_body)

        ! Allocate space
        allocate(qAq(this%quadrature_npoints))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectors and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute ag_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! it also computes in this case the term for the upper wing
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(A, source=this%head)
        call compute_auxiliary_and_A_2d_general(this%Binv, this%wingL, ag, this%wingU, bg, A, this%world)
        
        ! agymmetrize the elements of the macroscopic dielectric matrix and update it
        A   = this%symmetry%symmetryze_complex_tensor(A)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !         Compute the functions  (angular part)         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! HEAD: Compute \mathbf{\hat{q} A \hat{q}}  for the space grid
        call compute_qAq(A, this%xyz, this%world, qAq)

        ! WINGag : Compute \frac{ \mathbf{\hat{q}} \cdot ag_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} A \hat{q}}} and
        ! \frac{ \mathbf{\hat{q}} \cdot bg_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} A \hat{q}}} 
        ! Note that there is head missing term to allow for the efficient computation of the body term
        call compute_inverse_ag(this%xyz, ag, this%world, wingL_f)
        call compute_inverse_bg(this%xyz, bg, this%world, wingU_f)
        
        ! The body is directly computed as saving it to RAM is too intensive

        ! Note that wings functions are odd functions and thus their integral is zero, so no more is done regarding those.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Do the anisotropic averaging !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        allocate(this%idiel_wingL(nbasis), source=zzero)
        allocate(this%idiel_wingU(nbasis), source=zzero)
        allocate(this%idiel_body, source = this%Binv)
        
        this%idiel_head = head_2d(this%blm_coarse, this%weights, this%blm_fine, this%weights_fine, &
                                  this%radii, this%rmax2d, this%rcut, qAq, this%vr)

        ! Here we compute the body 
#if defined(USE_GPU) && defined(HAVEOMP5)
        allocate(wLwU_f(this%quadrature_npoints))
        !$omp target enter data map(to: this%idiel_body) map(alloc: wLwU_f)
        !$omp target teams distribute parallel do shared(this, qAq, wingL_f, wingU_f, nbasis) private(ii, jj, wLwU_f)
        do ii = 1, nbasis
            this%idiel_body(ii, ii) = this%idiel_body(ii, ii) - zone
            do jj = 1, nbasis
                wLwU_f(:) = wingL_f(:, jj) *  wingU_f(:, ii) 
                this%idiel_body(jj, ii) = this%idiel_body(jj, ii) + &
                                          body_2d(this%blm_coarse, this%weights, this%blm_fine, this%weights_fine, &
                                              this%radii, this%rmax2d, this%rcut, qAq, wLwU_f, this%vr)
            end do
        end do
        !$omp end target teams distribute parallel do
        !$omp target update from(this%idiel_body)
        !$omp target exit data map(delete: this%idiel_body, wLwU_f)
        deallocate(wLwU_f)
#else
        !$omp parallel shared(this, qAq, wingL_f, wingU_f, nbasis) private(ii, jj, wLwU_f)
        allocate(wLwU_f(this%quadrature_npoints))
        !$omp do schedule(dynamic)
        do ii = 1, nbasis
            this%idiel_body(ii, ii) = this%idiel_body(ii, ii) - zone
            do jj = 1, nbasis
                wLwU_f(:) = wingL_f(:, jj) *  wingU_f(:, ii) 
                this%idiel_body(jj, ii) = this%idiel_body(jj, ii) + &
                                        body_2d(this%blm_coarse, this%weights, this%blm_fine, this%weights_fine, &
                                            this%radii, this%rmax2d, this%rcut, qAq, wLwU_f, this%vr)
            end do
        end do 
        !$omp end do
        deallocate(wLwU_f)
        !$omp end parallel
#endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
        call this%world%register%remove("qAq")
        call this%world%register%remove("qag")
        call this%world%register%remove("qbg")
#endif

        deallocate(qAq, wingL_f, wingU_f)

    end subroutine compute_anisotropic_avg_scrcoulomb_2d_general


    module subroutine compute_anisotropic_avg_scrcoulomb_2d_hermitian(this)

        use idiel_linalg

        class(idiel_t), target, intent(inout) :: this

        ! Auxiliary vectors
        complex(aip), allocatable :: ag(:,:)
        
        ! Macroscopic dielectric matrix
        complex(aip), allocatable :: A(:,:)

        ! Function from which to compute the integral 
        complex(aip), target, allocatable :: qAq(:) 
        complex(aip), target, allocatable :: wingL_f(:,:)
        complex(aip), target, allocatable :: wLwU_f(:)

        ! Basis size
        integer(i32) :: nbasis
        
        ! Dummy indexes
        integer(i32) :: ii, jj

        ! Get the basis size
        nbasis = size(this%Binv, 1)

        ! Free final containers
        if (allocated(this%idiel_wingL)) deallocate(this%idiel_wingL)
        if (allocated(this%idiel_wingU)) deallocate(this%idiel_wingU)
        if (allocated(this%idiel_body))  deallocate(this%idiel_body)

        ! Allocate space
        allocate(qAq(this%quadrature_npoints))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectors and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute ag_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! it also computes in this case the term for the upper wing
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(A, source=this%head)
        call compute_auxiliary_and_A_2d_hermitian(this%Binv, this%wingL, ag, A, this%world)
        
        ! agymmetrize the elements of the macroscopic dielectric matrix and update it
        A   = this%symmetry%symmetryze_complex_tensor(A)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !         Compute the functions  (angular part)         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! HEAD: Compute \mathbf{\hat{q} A \hat{q}}  for the space grid
        call compute_qAq(A, this%xyz, this%world, qAq)

        ! WINGag : Compute \frac{ \mathbf{\hat{q}} \cdot ag_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} A \hat{q}}} and
        ! \frac{ \mathbf{\hat{q}} \cdot bg_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} A \hat{q}}} 
        ! Note that there is head missing term to allow for the efficient computation of the body term
        call compute_inverse_ag(this%xyz, ag, this%world, wingL_f)
        
        ! bghe body is directly computed as saving it to RAM is too intensive

        ! Note that wings functions are odd functions and thus their integral is zero, so no more is done regarding those.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Do the anisotropic averaging !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        allocate(this%idiel_wingL(nbasis), source=zzero)
        allocate(this%idiel_wingU(nbasis), source=zzero)
        allocate(this%idiel_body, source = this%Binv)
        
        this%idiel_head = head_2d(this%blm_coarse, this%weights, this%blm_fine, this%weights_fine, &
                                  this%radii, this%rmax2d, this%rcut, qAq, this%vr)
        
#if defined(USE_GPU) && defined(HAVEOMP5)
        allocate(wLwU_f(this%quadrature_npoints))
        !$omp target enter data map(to: this%idiel_body) map(alloc: wLwU_f)
        !$omp target teams distribute parallel do shared(this, qAq, wingL_f, nbasis) private(ii, jj, wLwU_f)
        do ii = 1, nbasis
            this%idiel_body(ii, ii) = this%idiel_body(ii, ii) - zone
            do jj = 1, nbasis
                wLwU_f(:) = wingL_f(:, jj) * conjg(wingL_f(:, ii))
                this%idiel_body(jj, ii) = this%idiel_body(jj, ii) + &
                                          body_2d(this%blm_coarse, this%weights, this%blm_fine, this%weights_fine, &
                                              this%radii, this%rmax2d, this%rcut, qAq, wLwU_f, this%vr)
            end do
        end do
        !$omp end target teams distribute parallel do
        !$omp target update from(this%idiel_body)
        !$omp target exit data map(delete: this%idiel_body, wLwU_f)
        deallocate(wLwU_f)
#else
        !$omp parallel shared(this, qAq, wingL_f, nbasis) private(ii, jj, wLwU_f)
        allocate(wLwU_f(this%quadrature_npoints))
        !$omp do schedule(dynamic)
        do ii = 1, nbasis
            this%idiel_body(ii, ii) = this%idiel_body(ii, ii) - zone
            do jj = 1, nbasis
                wLwU_f(:) = wingL_f(:, jj) * conjg(wingL_f(:, ii))
                this%idiel_body(jj, ii) = this%idiel_body(jj, ii) + &
                                        body_2d(this%blm_coarse, this%weights, this%blm_fine, this%weights_fine, &
                                            this%radii, this%rmax2d, this%rcut, qAq, wLwU_f, this%vr)
            end do
        end do 
        !$omp end do
        deallocate(wLwU_f)
        !$omp end parallel
#endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_GPU
        call this%world%register%remove("qAq")
        call this%world%register%remove("qag")
#endif
        deallocate(qAq, wingL_f)
        
    end subroutine compute_anisotropic_avg_scrcoulomb_2d_hermitian

    !> This function computes the screened Coulomb head contribution
    !> @param[in] blm_coarse - circular harmonics in the coarse angular mesh
    !> @param[in] wc - quadrature weights in the coarse angular mesh
    !> @param[in] blm_fine - circular harmonics in the fine angular mesh
    !> @param[in] wf - quadrature weights in the fine angular mesh
    !> @param[in] r  -  the radial mesh
    !> @param[in] rmax - the upper limit of the geometric integral
    !> @param[in] rcut - the cutoff radius
    !> @param[in] qAq  - the qAq product for all the angles in the coarse mesh
    !> @param[in] vr   - the 1-exp(-rcut*r) for the radial mesh
    !> @returns   w_head - the head component of the 2D screened Coulomb potential
    function head_2d(blm_coarse, wc, blm_fine, wf, r, rmax, rcut, qAq, vr) result(w_head)
        
        use idiel_cubic_spline

        complex(aip), intent(in) :: blm_coarse(:,:)
        real(aip), intent(in)    :: wc(:)
        complex(aip), intent(in) :: blm_fine(:,:)
        real(aip), intent(in)    :: wf(:)
        real(aip), intent(in)    :: r(:)
        real(aip), intent(in)    :: rmax(:)
        real(aip), intent(in)    :: rcut
        complex(aip), intent(in) :: qAq(:)
        complex(aip), intent(in) :: vr(:,:)

        complex(aip), allocatable         :: w_head_f(:,:)
        complex(aip)                      :: clm_head(ncir, nr)
        complex(aip)                      :: r_integral(size_mesh_2d_fine, ncir)
        integer(i32)                      :: ii, jj

        !! Splines
        real(aip), allocatable    :: x(:)
        complex(aip), allocatable :: splines(:,:), integrals(:)

        complex(aip) :: w_head

        allocate(w_head_f, source=vr)

        ! Compute the head component
        w_head_f(:,1) =  - (fourpi * rcut)**2 * qAq(:)

        do ii = 2, size(w_head_f,2)
            w_head_f(:, ii) =  - (w_head_f(:,ii) / r(ii))**2 * qAq(:) / ( 1.0_aip + w_head_f(:,ii) *  qAq(:))
        end do

        ! Expand the head component into Fourier
        do ii = 1, size(w_head_f, 2)
            call circ_harm_expansion(lmax, w_head_f(:,ii), wc, blm_coarse, clm_head(:,ii))
        end do

        ! Create interpolator for the Fourier components
        ! and integrate those as function of r multiplied by r 
        ! to obtain the rmax(phi) integral
        r_integral = zzero
        do ii = 1, ncir
            if ( maxval(abs(clm_head(ii,:))) < 1.0e-7_aip) cycle
             call init_cubic_spline_t(x, splines, integrals, r, r * clm_head(ii,:))
            do jj = 1, size_mesh_2d_fine
                 r_integral(jj, ii) = integrate(x, splines, integrals, 0.0_aip, rmax(jj))
            end do 
            deallocate(x, splines, integrals)
        end do

        ! Compute the integral
        w_head = zzero
        do ii = 1, ncir
            if ( maxval(abs(clm_head(ii,:))) < 1.0e-7_aip) cycle
            w_head = w_head + sum(wf(:) * blm_fine(:,ii) * r_integral(:,ii))
        end do

    end function head_2d

    !> This function computes the screened Coulomb head contribution coming from wings and head
    !> @param[in] blm_coarse - circular harmonics in the coarse angular mesh
    !> @param[in] wc - quadrature weights in the coarse angular mesh
    !> @param[in] blm_fine - circular harmonics in the fine angular mesh
    !> @param[in] wf - quadrature weights in the fine angular mesh
    !> @param[in] r  -  the radial mesh
    !> @param[in] rmax - the upper limit of the geometric integral
    !> @param[in] rcut - the cutoff radius
    !> @param[in] qAq  - the qAq product for all the angles in the coarse mesh
    !> @param[in] wLwU - the terms comming from wings at all the angles in the coarse mesh
    !> @param[in] vr   - the 1-exp(-rcut*r) for the radial mesh
    !> @returns   w_body - the correction body component of 2D screened Coulomb potential
    function body_2d(blm_coarse, wc, blm_fine, wf, r, rmax, rcut, qAq, wLwU, vr) result(w_body)

        use idiel_cubic_spline

        complex(aip), intent(in) :: blm_coarse(:,:)
        real(aip), intent(in)    :: wc(:)
        complex(aip), intent(in) :: blm_fine(:,:)
        real(aip), intent(in)    :: wf(:)
        real(aip), intent(in)    :: r(:)
        real(aip), intent(in)    :: rmax(:)
        real(aip), intent(in)    :: rcut
        complex(aip), intent(in) :: wLwU(:)
        complex(aip), intent(in) :: qAq(:)
        complex(aip), intent(in) :: vr(:,:)

#if defined(USE_GPU) && defined(HAVEOMP5)
        !$omp declare target
#endif  

        complex(aip), allocatable         :: w_body_f(:,:)
        complex(aip)                      :: clm_body(ncir, nr)
        complex(aip)                      :: r_integral(size_mesh_2d_fine, ncir)
        integer(i32)                      :: ii, jj

        !! Splines
        real(aip), allocatable    :: x(:)
        complex(aip), allocatable :: splines(:,:), integrals(:)

        complex(aip) :: w_body

        allocate(w_body_f, source=vr)

        ! Compute the body correction
        do ii = 1, size(w_body_f,2)
            w_body_f(:, ii) =  w_body_f(:,ii) * wLwU(:) / ( 1.0_aip + w_body_f(:,ii) *  qAq(:))
        end do 

        ! Expand the body component into Fourier
        do ii = 1, size(w_body_f,2)
            call circ_harm_expansion(lmax, w_body_f(:,ii), wc, blm_coarse, clm_body(:,ii))
        end do

        ! Create interpolator for the Fourier components
        ! and integrate those as function of r multiplied by r 
        ! to obtain the rmax(phi) integral
        r_integral(:,:) = zzero
        do ii = 1, ncir
            if (maxval(abs(clm_body(ii,:))) < 1.0e-7_aip) cycle
             call init_cubic_spline_t(x, splines, integrals, r, r * clm_body(ii,:))
            do jj = 1, size_mesh_2d_fine
                 r_integral(jj, ii) = integrate(x, splines, integrals, 0.0_aip, rmax(jj))
            end do 
             deallocate(x, splines, integrals)
        end do

        ! Compute the integral
        w_body = zzero
        do ii = 1, ncir
            if ( maxval(abs(clm_body(ii,:))) < 1.0e-7_aip) cycle
            w_body = w_body + sum(wf(:) * blm_fine(:,ii) * r_integral(:,ii))
        end do

    end function body_2d

! #if defined(USE_GPU) && defined(HAVEOMP5)
!     !> This function computes the body term inside a GPU
!     !> @param[in] world      - GPU-CPU register
!     !> @param[in] blm_coarse - circular harmonics in the coarse angular mesh
!     !> @param[in] wc - quadrature weights in the coarse angular mesh
!     !> @param[in] blm_fine - circular harmonics in the fine angular mesh
!     !> @param[in] wf - quadrature weights in the fine angular mesh
!     !> @param[in] r  -  the radial mesh
!     !> @param[in] rmax - the upper limit of the geometric integral
!     !> @param[in] rcut - the cutoff radius
!     !> @param[in] qAq  - the qAq product for all the angles in the coarse mesh
!     !> @param[in] vr   - the 1-exp(-rcut*r) for the radial mesh
!     !> @param[inout] body - the body term
!     !> @param[in] wL    - the lower wing term
!     !> @param[in] wU    - the upper wing term
!     subroutine body_2d_kernel(world, blm_coarse, wc, blm_fine, wf, radii, rmax2d, rcut, qAq, vr, body, wL, wU) 

!         use idiel_cubic_spline
!         use iso_c_binding
!         use omp_lib

!         type(gpu_world_t),    intent(inout)        :: world
!         complex(aip), target, intent(in)           :: blm_coarse(:,:)
!         real(aip),    target, intent(in)           :: wc(:)
!         complex(aip), target, intent(in)           :: blm_fine(:,:)
!         real(aip),    target, intent(in)           :: wf(:)
!         real(aip),    target, intent(in)           :: radii(:)
!         real(aip),    target, intent(in)           :: rmax2d(:)
!         real(aip),    target, intent(in)           :: rcut
!         complex(aip), target, intent(in)           :: qAq(:)
!         complex(aip), target, intent(in)           :: vr(:,:)
!         complex(aip), target, intent(inout)        :: body(:,:)
!         complex(aip), target, intent(in)           :: wL(:,:)
!         complex(aip), target, optional, intent(in) :: wU(:,:)

!         integer(i32) :: nbasis 
!         integer(i32) :: nquadrature_points
!         integer(i32) :: nteams, id_team 
!         integer(i32) :: ii, jj, kk, ll

!         complex(aip), allocatable :: wLwU(:,:)
!         complex(aip), allocatable :: body_f(:,:,:)
!         ! complex(aip), allocatable :: clm(:,:,:)

!         ! Get sizes
!         nbasis             = size(body, 1)
!         nquadrature_points = size(qAq)
!         nteams             = max(omp_get_max_teams(), 24)

!         ! Allocate
!         allocate(wLwU(nquadrature_points, nteams))
!         allocate(body_f(nquadrature_points, nr, nteams))

!         ! Upload the body 
!         call world%register%alloc("body", size(body) * c_sizeof(zzero), world%get_device())
!         call world%register%assoc("body", C_loc(body))
!         call world%register%to_device("body")

!         ! Associate
!         call world%register%assoc("blm_coarse"   , C_loc(blm_coarse)  )
!         call world%register%assoc("weights"      , C_loc(wc)     )
!         call world%register%assoc("blm_fine"     , C_loc(blm_fine)    )
!         call world%register%assoc("weights_fine" , C_loc(wf))
!         call world%register%assoc("radii"        , C_loc(radii)       )
!         call world%register%assoc("rmax2d"       , C_loc(rmax2d)      )
!         call world%register%assoc("vr"           , C_loc(vr)          )

!         call world%register%assoc("qag", C_loc(wL))
!         call world%register%assoc("qAq", C_loc(qAq))
!         if (present(wU)) then
!             call world%register%assoc("qbg", C_loc(wU))
!         end if


!         ! Main loop in the GPU
!         if (present(wU)) then
!             error stop "Not implemented"
!         else 
!             !$omp target data map(always, to: wLwU, body_f)
!             !$omp target teams distribute private(id_team) collapse(2) num_teams(nteams)
!             do ii = 1, nbasis
!                 do jj = 1, nbasis
!                     ! Get team id (add 1 because is 0 based)
!                     id_team = omp_get_team_num() + 1
!                     ! Update diagonal
!                     if (ii == jj) then 
!                         body(jj, ii) = body(jj, ii) - zone
!                     end if
!                     !$omp workshare
!                     wLwU(:, id_team) = wL(:, jj) * conjg(wL(:, ii))
!                     !$omp end workshare
                    
!                     ! Compute the body correction
!                     !$omp parallel do simd collapse(2) private(kk, ll)
!                     do kk = 1, nr
!                         do ll = 1, nquadrature_points
!                             body_f(ll, kk, id_team) =  body_f(ll, kk, id_team) * wLwU(ll, id_team) / ( 1.0_aip + body_f(ll, kk, id_team) *  qAq(ll))
!                         end do
!                     end do 
!                     !$omp end parallel do simd

!                 end do
!             end do
!             !$omp end target teams distribute
!             !$omp end target data
!         end if

!         ! Update
!         call world%register%from_device("body")
!         write(*,*) body(100,100)
!         ! Clean
!         call world%register%deassoc("blm_coarse"  )
!         call world%register%deassoc("weights"     )
!         call world%register%deassoc("blm_fine"    )
!         call world%register%deassoc("weights_fine")
!         call world%register%deassoc("radii"       )
!         call world%register%deassoc("rmax2d"      )
!         call world%register%deassoc("vr"          )

!         call world%register%remove("body")

!         return

!     end subroutine body_2d_kernel

! #endif 

end submodule idiel_2D
