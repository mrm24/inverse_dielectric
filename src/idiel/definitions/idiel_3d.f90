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
submodule (idiel) idiel_3D

    implicit none

contains

    module subroutine compute_anisotropic_avg_inversedielectric_3d(this, hermitian)
            
        class(idiel_t), intent(inout) :: this
        logical, intent(in)           :: hermitian

        if (this%dim /= 3) then
            error stop "Error(idiel_t%compute_anisotropic_avg_inversedielectic_3d: this procedure can only be used with 3D materials"
        end if

        if (hermitian) then
            call this%compute_anisotropic_avg_hermitian()
        else
            call this%compute_anisotropic_avg_general()
        end if

    end subroutine compute_anisotropic_avg_inversedielectric_3d
    
    module subroutine compute_anisotropic_avg_hermitian(this)

        use idiel_linalg

        class(idiel_t), target, intent(inout) :: this

        ! Auxiliary vectors
        complex(aip), allocatable :: S(:,:)
        
        ! Macroscopic dielectric matrix
        complex(aip), allocatable :: L(:,:)

        ! Function from which to compute the integral 
        complex(aip), target, allocatable :: head_f(:) 
        complex(aip), target, allocatable :: wingL_f(:,:)
        complex(aip), target, allocatable :: body_f(:)

        ! Harmonic expansion coefficients
        complex(aip), target, allocatable :: clm_head(:)
        complex(aip), target, allocatable :: clm_body(:)

        ! Basis size
        integer(i32) :: nbasis
        
        ! Dummy indexes
        integer(i32) :: ii, jj, kk

        ! For device
        integer(i32)              :: nr
        complex(aip), allocatable :: ylm(:,:)
        real(aip), allocatable    :: weights(:)
        complex(aip), allocatable :: angular_integrals(:)
        complex(aip), allocatable :: body(:,:)

        ! Get the basis size
        nbasis = size(this%Binv, 1)

        ! Free final containers
        if (allocated(this%idiel_wingL)) deallocate(this%idiel_wingL)
        if (allocated(this%idiel_wingU)) deallocate(this%idiel_wingU)
        if (allocated(this%idiel_body))  deallocate(this%idiel_body)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectors and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute S_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(L, source=this%head)
        call compute_auxiliary_and_macroscopic_3d_hermitian(this%Binv, this%wingL, S, L, this%world)

        ! Symmetrize the elements of the macroscopic dielectric matrix and update it
        L   = this%symmetry%symmetryze_complex_tensor(L)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               Compute the functions                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! HEAD: Compute \frac{\omega(\Omega)}{\mathbf{\hat{q} L \hat{q}}}  for the space grid
        call compute_inverse_head(L, this%xyz, this%world, head_f)

        ! WINGS : Compute \frac{ \mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} and
        ! \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} 
        ! Note that there is head missing term to allow for the efficient computation of the body term
        call compute_inverse_wingL(this%xyz, S, this%world, wingL_f)
        
        ! The body is directly computed as saving it to RAM is too intensive
        ! Note that wings functions are odd functions and thus their integral is zero, so no more work is done regarding those.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Do the anisotropic averaging !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(this%idiel_wingL(nbasis), source=zzero)
        allocate(this%idiel_wingU(nbasis), source=zzero)
        allocate(this%idiel_body, source = this%Binv)

        allocate(clm_head(nsph_pair))
        call sph_harm_expansion(nsph_pair, head_f, this%weights, this%ylm, clm_head)
        this%idiel_head  = sum(clm_head(:) * this%angular_integrals(:))
        deallocate(clm_head)
        ! Here we compute the body 
        ! \frac{[\mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})] [\mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})]}{\mathbf{\hat{q} L \hat{q}}} 
#if defined(USE_GPU) && defined(HAVEOMP5)
        ! We need to use associations for OMP to understand
        associate(quadrature_npoints => this%quadrature_npoints, body => this%idiel_body, &
                  ylm => this%ylm, weights => this%weights, angular_integrals => this%angular_integrals, &
                  world => this%world)

            allocate(body_f(quadrature_npoints), clm_body(nsph_pair))

            call world%register%alloc("body", size(body) * c_sizeof(zzero), world%get_device())
            call world%register%assoc("body", C_loc(body))
            call world%register%to_device("body")

            call world%register%alloc("body_f", size(body_f) * c_sizeof(zzero), world%get_device())
            call world%register%assoc("body_f", C_loc(body_f))

            call world%register%alloc("clm_body", size(clm_body) * c_sizeof(zzero), world%get_device())
            call world%register%assoc("clm_body", C_loc(clm_body))

            call world%register%assoc("ylm", C_loc(ylm))
            call world%register%assoc("weights", C_loc(weights)) 
            call world%register%assoc("angular_integrals", C_loc(angular_integrals))

            call world%register%assoc("qS", C_loc(wingL_f))
            call world%register%assoc("invqLq", C_loc(head_f))

            !$omp target teams distribute private(ii, jj, body_f, clm_body)
            do ii = 1, nbasis
                do jj = 1, nbasis
                    body_f(:) = head_f(:) * wingL_f(:, jj) * conjg(wingL_f(:, ii))
                    ! This call itself is parallelized so each team can use its threads
                    ! to solve it
                    call sph_harm_expansion(nsph_pair, body_f, weights, ylm, clm_body)
                    body(jj,ii) = body(jj,ii) + sum(clm_body(:) * angular_integrals(:))
                end do
            end do
            !$omp end target teams distribute
            
            call world%register%from_device("body")

            call world%register%deassoc("ylm")
            call world%register%deassoc("weights")
            call world%register%deassoc("angular_integrals")
            
            call world%register%remove("body")
            call world%register%remove("body_f")
            call world%register%remove("clm_body")

            deallocate(body_f, clm_body)

        end associate
#else
        !$omp parallel shared(this, head_f, wingL_f, nbasis) private(ii, jj, body_f, clm_body)
        allocate(body_f(this%quadrature_npoints), clm_body(nsph_pair))
        !$omp do schedule(dynamic) collapse(2)
        do ii = 1, nbasis
            do jj = 1, nbasis
                body_f(:) = head_f(:) * wingL_f(:, jj) * conjg(wingL_f(:, ii))
                call sph_harm_expansion(nsph_pair, body_f, this%weights, this%ylm, clm_body)
                this%idiel_body(jj,ii) = this%idiel_body(jj,ii) + &
                    sum(clm_body(:) * this%angular_integrals(:))
            end do
        end do 
        !$omp end do
        deallocate(body_f, clm_body)
        !$omp end parallel
#endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
        call this%world%register%remove("invqLq")
        call this%world%register%remove("qS")
#endif

        deallocate(head_f, wingL_f)

    end subroutine compute_anisotropic_avg_hermitian

    module subroutine compute_anisotropic_avg_general(this)

        use idiel_linalg

        class(idiel_t), target, intent(inout) :: this

        ! Auxiliary vectors
        complex(aip), allocatable :: S(:,:), T(:,:)
        
        ! Macroscopic dielectric matrix
        complex(aip), allocatable :: L(:,:)

        ! Function from which to compute the integral 
        complex(aip), target, allocatable :: head_f(:) 
        complex(aip), target, allocatable :: wingL_f(:,:)
        complex(aip), target, allocatable :: wingU_f(:,:)
        complex(aip), target, allocatable :: body_f(:)

        ! Harmonic expansion coefficients
        complex(aip), target, allocatable :: clm_head(:)
        complex(aip), target, allocatable :: clm_body(:)

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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectorws and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute S_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! it also computes in this case the term for the upper wing
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(L, source=this%head)
        call compute_auxiliary_and_macroscopic_3d_general(this%Binv, this%wingL, S, this%wingU, T, L, this%world)
        
        ! Symmetrize the elements of the macroscopic dielectric matrix and update it
        L   = this%symmetry%symmetryze_complex_tensor(L) 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               Compute the functions                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! HEAD: Compute \frac{\omega(\Omega)}{\mathbf{\hat{q} L \hat{q}}}  for the space grid
        call  compute_inverse_head(L, this%xyz, this%world, head_f)

        ! WINGS : Compute \frac{ \mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} and
        ! \frac{ \mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})}{\mathbf{\hat{q} L \hat{q}}} 
        ! Note that there is head missing term to allow for the efficient computation of the body term
        call compute_inverse_wingL(this%xyz, S, this%world, wingL_f)
        call compute_inverse_wingU(this%xyz, T, this%world, wingU_f)
        
        ! The body is directly computed as saving it to RAM is too intensive

        ! Note that wings functions are odd functions and thus their integral is zero, so no more is done regarding those.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Do the anisotropic averaging !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        allocate(this%idiel_wingL(nbasis), source=zzero)
        allocate(this%idiel_wingU(nbasis), source=zzero)
        allocate(this%idiel_body, source = this%Binv)

        allocate(clm_head(nsph_pair))
        call sph_harm_expansion(nsph_pair, head_f, this%weights, this%ylm, clm_head)
        this%idiel_head  = sum(clm_head(:) * this%angular_integrals(:))
        deallocate(clm_head)

        ! Here we compute the body 
        ! \frac{[\mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})] [\mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})]}{\mathbf{\hat{q} L \hat{q}}} 
#if defined(USE_GPU) && defined(HAVEOMP5)
        ! We need to use associations for OMP to understand
        associate(quadrature_npoints => this%quadrature_npoints, body => this%idiel_body, &
            ylm => this%ylm, weights => this%weights, angular_integrals => this%angular_integrals, &
            world => this%world)

        allocate(body_f(quadrature_npoints), clm_body(nsph_pair))

        call world%register%alloc("body", size(body) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("body", C_loc(body))
        call world%register%to_device("body")

        call world%register%alloc("body_f", size(body_f) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("body_f", C_loc(body_f))

        call world%register%alloc("clm_body", size(clm_body) * c_sizeof(zzero), world%get_device())
        call world%register%assoc("clm_body", C_loc(clm_body))

        call world%register%assoc("ylm", C_loc(ylm))
        call world%register%assoc("weights", C_loc(weights)) 
        call world%register%assoc("angular_integrals", C_loc(angular_integrals))

        call world%register%assoc("qS", C_loc(wingL_f))
        call world%register%assoc("qT", C_loc(wingU_f))
        call world%register%assoc("invqLq", C_loc(head_f))

        !$omp target teams distribute private(ii, jj, body_f, clm_body)
        do ii = 1, nbasis
            do jj = 1, nbasis
                body_f(:) = head_f(:) * wingL_f(:, jj) * wingU_f(:, ii)
                ! This call itself is parallelized so each team can use its threads
                ! to solve it
                call sph_harm_expansion(nsph_pair, body_f, weights, ylm, clm_body)
                body(jj,ii) = body(jj,ii) + sum(clm_body(:) * angular_integrals(:))
            end do
        end do
        !$omp end target teams distribute

        call world%register%from_device("body")

        call world%register%deassoc("ylm")
        call world%register%deassoc("weights")
        call world%register%deassoc("angular_integrals")

        call world%register%remove("body")
        call world%register%remove("body_f")
        call world%register%remove("clm_body")

        deallocate(body_f, clm_body)

        end associate
#else
        !$omp parallel shared(this, head_f, wingL_f, wingU_f, nbasis) private(ii, jj, body_f, clm_body)
        allocate(body_f(this%quadrature_npoints),  clm_body(nsph_pair))
        !$omp do schedule(dynamic) collapse(2)
        do ii = 1, nbasis
            do jj = 1, nbasis
                body_f(:) = head_f(:) * wingL_f(:, jj) * wingU_f(:, ii)
                call sph_harm_expansion(nsph_pair, body_f, this%weights, this%ylm, clm_body)
                this%idiel_body(jj,ii) = this%idiel_body(jj,ii) + &
                    sum(clm_body(:) * this%angular_integrals(:))
            end do
        end do 
        !$omp end do
        deallocate(body_f,  clm_body)
        !$omp end parallel
#endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
        call this%world%register%remove("invqLq")
        call this%world%register%remove("qS")
        call this%world%register%remove("qT")
#endif

        deallocate(head_f, wingL_f, wingU_f)

    end subroutine compute_anisotropic_avg_general

end submodule idiel_3D
