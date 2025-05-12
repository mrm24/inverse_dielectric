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

#include "offload.fpp"



#ifdef DEVICEOFFLOAD
    use omp_lib
#define _omp_get_team_num omp_get_team_num()
#else 
#define _omp_get_team_num 1
#endif

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
        complex(aip), target, allocatable :: body_f(:,:)

        ! Harmonic expansion coefficients
        complex(aip), target, allocatable :: clm_head(:)
        complex(aip), target, allocatable :: clm_body(:,:)

        ! Basis size
        integer(i32) :: nbasis
        
        ! Dummy indexes
        integer(i32) :: ii, jj, kk, my_team

        ! For device aware compilation the team id
        integer(i32) :: my_team

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
        
        allocate(body_f(this%quadrature_npoints, this%world%get_num_teams()))
        allocate(clm_body(nsph_pair, this%world%get_num_teams()))
        
        OMP_OFFLOAD target enter data map(to: this%idiel_body, clm_body, body_f, head_f, wingL_f)

        OMP_OFFLOAD target
        OMP_OFFLOAD teams distribute num_teams(this%world%get_num_teams()) &
        OMP_NO_OFFLOAD parallel do &
        !$omp collapse(2) default(none) shared(this, head_f, wingL_f, nbasis, body_f, clm_body) private(ii, jj, kk, my_team)
        do ii = 1, nbasis
            do jj = 1, nbasis
                my_team = _omp_get_team_num

                OMP_OFFLOAD parallel do default(none) shared(body_f, head_f, wingL_f, ii, jj, my_team) private(kk)
                do kk = 1, this%quadrature_npoints
                    body_f(kk, my_team) = head_f(kk) * wingL_f(kk, jj) * conjg(wingL_f(kk, ii))
                end do
                OMP_OFFLOAD end parallel do
                
                call sph_harm_expansion(nsph_pair, body_f(:,my_team), this%weights, this%ylm, clm_body(:,my_team))
                
                this%idiel_body(jj,ii) = this%idiel_body(jj,ii) + &
                    sum(clm_body(:, my_team) * this%angular_integrals(:))
            end do
        end do 
        OMP_NO_OFFLOAD end parallel do
        OMP_OFFLOAD end teams distribute
        OMP_OFFLOAD end target

        OMP_OFFLOAD target update from(this%idiel_body)
        OMP_OFFLOAD target exit data map(delete: this%idiel_body, clm_body, body_f, head_f, wingL_f)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(body_f, clm_body)
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
        complex(aip), target, allocatable :: body_f(:,:)

        ! Harmonic expansion coefficients
        complex(aip), target, allocatable :: clm_head(:)
        complex(aip), target, allocatable :: clm_body(:,:)

        ! Basis size
        integer(i32) :: nbasis
        
        ! Dummy indexes
        integer(i32) :: ii, jj, kk

        ! For device aware compilation the team id
        integer(i32) :: my_team

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

        allocate(body_f(this%quadrature_npoints, this%world%get_num_teams()))
        allocate(clm_body(nsph_pair, this%world%get_num_teams()))

        OMP_OFFLOAD target enter data map(to: this%idiel_body, clm_body, body_f, head_f, wingL_f, wingU_f)

        OMP_OFFLOAD target
        OMP_OFFLOAD teams distribute num_teams(this%world%get_num_teams()) &
        OMP_NO_OFFLOAD parallel do &
        !$omp collapse(2) default(none) shared(this, head_f, wingL_f, wingU_f, nbasis, body_f, clm_body) private(ii, jj, kk, my_team)
        do ii = 1, nbasis
            do jj = 1, nbasis
                my_team = _omp_get_team_num

                OMP_OFFLOAD parallel do default(none) shared(body_f, head_f, wingL_f, wingU_f, ii, jj, my_team) private(kk)
                do kk = 1, this%quadrature_npoints
                    body_f(kk, my_team) = head_f(kk) * wingL_f(kk, jj) * wingU_f(kk, ii)
                end do
                OMP_OFFLOAD end parallel do

                call sph_harm_expansion(nsph_pair, body_f(:,my_team), this%weights, this%ylm, clm_body(:,my_team))

                this%idiel_body(jj,ii) = this%idiel_body(jj,ii) + &
                    sum(clm_body(:, my_team) * this%angular_integrals(:))
            end do
        end do
        OMP_NO_OFFLOAD end parallel do
        OMP_OFFLOAD end teams distribute
        OMP_OFFLOAD end target

        OMP_OFFLOAD target update from(this%idiel_body)
        OMP_OFFLOAD target exit data map(delete: this%idiel_body, clm_body, body_f, head_f, wingL_f, wingU_f)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(body_f, clm_body)
        deallocate(head_f, wingL_f, wingU_f)


    end subroutine compute_anisotropic_avg_general

end submodule idiel_3D
