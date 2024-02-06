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
submodule (idiel) idiel_3D

    implicit none

contains

    module subroutine compute_anisotropic_avg_inversedielectic_3d(this, hermitian)
            
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

    end subroutine compute_anisotropic_avg_inversedielectic_3d
    
    module subroutine compute_anisotropic_avg_hermitian(this)

        use idiel_linalg

        class(idiel_t), intent(inout) :: this

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
        complex(r64) :: clm_head(nsph_pair)
        complex(r64) :: clm_body(nsph_pair) 
        real(r64)    :: error 

        ! Basis size
        integer(i64) :: nbasis
        
        ! Dummy indexes
        integer(i64) :: ii, jj, kk

        ! Get the basis size
        nbasis = size(this%Binv, 1)

        ! Free final containers
        if (allocated(this%idiel_wingL)) deallocate(this%idiel_wingL)
        if (allocated(this%idiel_wingU)) deallocate(this%idiel_wingU)
        if (allocated(this%idiel_body))  deallocate(this%idiel_body)

        ! Allocate space
        allocate(head_f(this%quadrature_npoints))
        allocate(wingL_f(this%quadrature_npoints, nbasis))
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectorws and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute S_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(L, source=this%head)
        call compute_auxiliary_and_macroscopic_3d(this%Binv, this%wingL, S, ref_S, L, ref_L, this%world)
        
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
        
        allocate(this%idiel_wingL(nbasis), source=zzero)
        allocate(this%idiel_wingU(nbasis), source=zzero)
        allocate(this%idiel_body, source = this%Binv)
        
        call sph_harm_expansion(nsph_pair, head_f, this%weights, this%ylm, clm_head)
        this%idiel_head  = sum(clm_head(:) * this%angular_integrals(:))

        error = maxval(abs(clm_head(nsph_pair-41:nsph_pair)))/abs(clm_head(1)) 
        if (error > 1.0e-4_r64) then
            write(*,*) "Warning (compute_anisotropic_avg_hermitian) the expansion coefficient can be not large enough"
        end if
        
        ! Here we compute the body 
        ! \frac{[\mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})] [\mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})]}{\mathbf{\hat{q} L \hat{q}}} 
        !$omp parallel shared(this, head_f, wingL_f, nbasis) private(ii, jj, body_f, clm_body)
        allocate(body_f(this%quadrature_npoints))
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
        deallocate(body_f)
        !$omp end parallel

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(head_f, wingL_f)

    end subroutine compute_anisotropic_avg_hermitian

    module subroutine compute_anisotropic_avg_general(this)

        use idiel_linalg

        class(idiel_t), intent(inout) :: this

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
        complex(r64) :: clm_head(nsph_pair)
        complex(r64) :: clm_body(nsph_pair) 
        real(r64) :: error

        ! Basis size
        integer(i64) :: nbasis
        
        ! Dummy indexes
        integer(i64) :: ii, jj

        ! Get the basis size
        nbasis = size(this%Binv, 1)

        ! Free final containers
        if (allocated(this%idiel_wingL)) deallocate(this%idiel_wingL)
        if (allocated(this%idiel_wingU)) deallocate(this%idiel_wingU)
        if (allocated(this%idiel_body))  deallocate(this%idiel_body)

        ! Allocate space
        allocate(head_f(this%quadrature_npoints))
        allocate(wingL_f(this%quadrature_npoints, nbasis))
        allocate(wingU_f(this%quadrature_npoints, nbasis))
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Auxiliary vectorws and macroscopic dielectric matrix !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute S_{\alpha}(\mathbf{G}) = \sum_{\mathbf{G'\neq 0}} B^{-1}_{\mathbf{GG'}} U_{\alpha}(\mathbf{G}) (Eq. B.13)
        ! it also computes in this case the term for the upper wing
        ! and the local field effects of the head (Eq. B.14) in 10.1016/j.cpc.2006.07.018
        allocate(L, source=this%head)
        call compute_auxiliary_and_macroscopic_3d(this%Binv, this%wingL, S, ref_S, L, ref_L, this%world, this%wingU, T, ref_T)
        
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
        
        allocate(this%idiel_wingL(nbasis), source=zzero)
        allocate(this%idiel_wingU(nbasis), source=zzero)
        allocate(this%idiel_body, source = this%Binv)
        
        call sph_harm_expansion(nsph_pair, head_f, this%weights, this%ylm, clm_head)
        this%idiel_head  = sum(clm_head(:) * this%angular_integrals(:))

        error = maxval(abs(clm_head(nsph_pair-41:nsph_pair)))/abs(clm_head(1)) 
        if (error > 1.0e-4_r64) then
            write(*,*) "Warning (compute_anisotropic_avg_general) the expansion coefficient can be not large enough"
        end if
        
        ! Here we compute the body 
        ! \frac{[\mathbf{\hat{q}} \cdot T_{\alpha}(\mathbf{G})] [\mathbf{\hat{q}} \cdot S_{\alpha}(\mathbf{G})]}{\mathbf{\hat{q} L \hat{q}}} 
        !$omp parallel shared(this, head_f, wingL_f, wingU_f, nbasis) private(ii, jj, body_f, clm_body)
        allocate(body_f(this%quadrature_npoints))
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
        deallocate(body_f)
        !$omp end parallel

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            Clean              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(head_f, wingL_f, wingU_f)

    end subroutine compute_anisotropic_avg_general

end submodule idiel_3D
