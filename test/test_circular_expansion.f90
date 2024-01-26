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
!> Tests the circular expansion of a function similar to the head component
program test_circular_expansion

    use idiel_constants, only: r64, i64, twopi, pi, fourpi, zzero
    use idiel_circular_harmonics
    use idiel_circle_quadrature 

    implicit none

    integer(i64), parameter :: msize = 55_i64
    integer(i64), parameter :: lmax  = 10_i64
    real(r64),    parameter :: tolerance = 1.0e-12_r64
    real(r64), allocatable :: rphi(:,:), w(:), xyz(:,:)
    complex(r64), allocatable :: blm(:,:)
    real(r64) :: integral
    integer(i64) :: i


    write(*,*) '[TEST : idiel_circular_harmonics and idiel_circle_quadrature]' 
    write(*,*)  '  * Check exact integration of circular basis'
    ! We compute the angular mesh
    call compute_angular_mesh_gauss_legendre(msize, rphi, w, xyz)

    ! We compute the circular harmonics up to 22
    call circ_harm(2*msize, rphi(:,2), blm)

    ! Check that the quadrature integrates the circular harmonics
    integral = real(sum(w * blm(:,1)))
    if ( abs(integral-2.50662827463100_r64)/2.50662827463100_r64 <= tolerance) then
        write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: PASSED]'
    else   
        write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: FAILED]'
        stop 1
    end if

    do i = 2, 41
        integral = real(sum(w * blm(:,i)))
        if ( integral <= tolerance) then
            write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: PASSED]'
        else   
            write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: FAILED]'
            stop 1
        end if
    end do

    stop 0

end program test_circular_expansion

! contains 
! Check the quadrature coefficients
! do i = 1, nr/2
!     radii(i) = -log(-(i - 1) * dr1 + 1)
!     radii(i + nr/2) = xcut + i * dr2
! end do
! head = calculate_headf(radii,rphi(:,2))
! ! Test the expansion
! allocate(clm(2*lmax+1,nr))
! do i = 1, nr
!     call circ_harm_expansion(lmax, head(:,i), w, blm(:,1:2*lmax+1), clm(:,i))
!     write(*,*) i, radii(i), real(head(1,i)), real(clm(1,i)), real(clm(5,i))
! end do
! function calculate_headf(r,ang) result(h)
!     !> Radius
!     real(r64), intent(in) :: r(:)
!     !> Angles
!     real(r64), intent(in) :: ang(:)

!     complex(r64) :: Axx, Ayy, vrii
!     integer(i64) :: ii, jj, nrad, nang
!     complex(r64), allocatable :: qAq(:), vr(:), h(:,:)

!     nrad = size(r)
!     nang = size(ang) 

!     Axx = cmplx(1.0_r64,0.0_r64,r64)
!     Ayy = cmplx(3.0_r64,0.0_r64,r64)

!     allocate(vr(nrad))
!     vr = fourpi * (1.0_r64 - exp(-rcut*r))

!     allocate(qAq(nang))
        
!     qAq = cos(ang)**2 * Axx + sin(ang)**2 * Ayy

!     allocate(h(nang,nrad))
        
!     h(:,1) =  - (fourpi * rcut)**2 * qAq(:)

!     do ii = 2, nrad
!         do jj = 1, nang
!             h(jj,ii) =  - (vr(ii) / r(ii))**2 * qAq(jj) / ( 1 + vr(ii) *  qAq(jj))
!         end do 
!     end do 
        
! end function calculate_headf