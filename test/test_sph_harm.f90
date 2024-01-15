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
!> Tests the spherical harmonics generation
program test_sph_harm

    use idiel_constants, only: r64, i64
    use idiel_sph_quadrature, only: compute_angular_mesh_lebedev_21
    use idiel_spherical_harmonics,  only: sph_harm

    implicit none 

    ! Spherical harmonics
    integer(i64), parameter :: lmax = 5_i64
    integer(i64), parameter :: nsph = (lmax+1)**2
    complex(r64), allocatable :: ylm(:,:)
    integer :: i
    real(r64) :: real_ylm(nsph), imag_ylm(nsph)
    complex(r64) :: ylm_ref(170, nsph)

    ! Tolerance
    real(r64), parameter :: tolerance = 1.0e-12_r64
    real(r64) :: rdiff

    ! Weights and points
    integer(i64), parameter :: np = 170_i64
    real(r64), allocatable  :: ang_mesh(:,:)
    real(r64), allocatable  :: xyz_mesh(:,:)
    real(r64), allocatable  :: weights(:)

    ! Print info
    write(*,*) '[TEST : sph_harm]' 
    write(*,*) ' * kind : regression test against precomputed data'

    ! Compute the integral points and weights
    call compute_angular_mesh_lebedev_21(ang_mesh, weights, xyz_mesh)

    ! Compute the spherical harmonics
    call sph_harm(lmax, ang_mesh, ylm)

    ! Read data to compare to
    open(unit=101, file='real_part.txt', status='old', action='read')
    open(unit=102, file='imag_part.txt', status='old', action='read')
    do i = 1, np
        read(101, *) real_ylm(:)
        read(102, *) imag_ylm(:)
        ylm_ref(i,:) = cmplx( real_ylm, imag_ylm, r64)
    end do
    close(101)
    close(102)

    rdiff = sum(abs(ylm - ylm_ref))/sum(abs(ylm_ref))
    write(*,'(A, e20.13)')  '  * Regression spherical harmonics result (relative difference): ', rdiff
    if ( rdiff < tolerance ) then
        write(*,*)  '[TEST : sph_harm: PASSED]'
    else
        write(*,*)  '[TEST : sph_harm: FAILED]'
        stop 1
    end if

    stop 0

end program test_sph_harm
