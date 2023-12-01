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
!> Tests the spherical quadrature using Lebedev-mesh works by comparing to objective functions

program test_quadrature_lebedev
    
    use idiel_constants, only: r64, i64, pi, twopi
    use idiel_sph_quadrature, only: compute_angular_mesh_lebedev

    implicit none

    ! Tolerance
    real(r64), parameter :: tolerance = 1.0e-10_r64
    ! The alpha parameter
    real(r64), parameter :: alpha = 12.0_r64
    ! The exact integrals of test functions 
    real(r64), parameter :: I1 = 216_r64 * pi / 35.0_r64
    real(r64), parameter :: I2 = 6.6961822200736179523_r64
    real(r64), parameter :: I3 = 4.0_r64 * pi / alpha
    real(r64), parameter :: I4 = 4.0_r64 * pi / alpha
    real(r64), parameter :: I5 = 4.0_r64 * pi / alpha

    ! Weights and points
    real(r64), allocatable :: ang_mesh(:,:)
    real(r64), allocatable :: xyz_mesh(:,:)
    real(r64), allocatable :: weights(:)

    ! Storage
    real(r64), allocatable :: f(:)
    real(r64) :: nuidiel_integral, rdiff

    write(*,*) '[TEST : compute_ang_lebedev]' 
    write(*,*) ' * kind : regression test against precomputed data'

    ! Compute the integral points and weights
    call compute_angular_mesh_lebedev(ang_mesh, weights, xyz_mesh)

    ! Compute integral 1 and check
    f = f1(xyz_mesh)
    nuidiel_integral = sum(weights * f)
    rdiff = abs( nuidiel_integral - I1 ) / I1
    write(*,'(A, E11.6)')  '  * Regression (f1) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f1: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f1: FAILED]'
        stop 1
    end if

    !Compute integral 2 and check (this one has lower convergence for all quadratures, thus we reduce convergence criteria)
    f = f2(xyz_mesh)
    nuidiel_integral = sum(weights * f)
    rdiff = abs( nuidiel_integral - I2 ) / I2
    write(*,'(A, E11.6)')  '  * Regression (f2) result (relative difference): ', rdiff
    if (rdiff .lt. 1.0e+3_r64 * tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f2: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f2: FAILED]'
        stop 1
    end if

    ! Compute integral 3 and check
    f = f3(xyz_mesh)
    nuidiel_integral = sum(weights * f)
    rdiff = abs( nuidiel_integral - I3 ) / I3
    write(*,'(A, E11.6)')  '  * Regression (f3) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f3: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f3: FAILED]'
        stop 1
    end if

    ! Compute integral 4 and check
    f = f4(xyz_mesh)
    nuidiel_integral = sum(weights * f)
    rdiff = abs( nuidiel_integral - I4 ) / I4
    write(*,'(A, E11.6)')  '  * Regression (f4) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f4: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f4: FAILED]'
        stop 1
    end if

    ! Compute integral 5 and check
    f = f5(xyz_mesh)
    nuidiel_integral = sum(weights * f)
    rdiff = abs( nuidiel_integral - I5 ) / I5
    write(*,'(A, E11.6)')  '  * Regression (f5) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f5: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f5: FAILED]'
        stop 1
    end if

    stop 0
contains

    !> Test function 1
    pure function f1(r)
        !> Mesh in which to compute
        real(r64), intent(in)  :: r(:,:)
        !> Result
        real(r64), allocatable :: f1(:)
        !> Elements
        real(r64), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f1 = 1.0_r64 + x + y**2 + y * x**2 + x**4 + y**5 + (x*y*z)**2 

    end function f1

    !> Test function 2
    pure function f2(r)
        !> Mesh in which to compute
        real(r64), intent(in)  :: r(:,:)
        !> Result
        real(r64), allocatable :: f2(:)
        !> Elements
        real(r64), allocatable  :: x(:), y(:), z(:)
        real(r64), allocatable, dimension(:) :: exp1, exp2, exp3, exp4

        x = 9.0_r64 * r(:,1)
        y = 9.0_r64 * r(:,2)
        z = 9.0_r64 * r(:,3)

        exp1 = 0.75 * exp( -0.25 * (x - 2.0)**2)     * exp(-0.25*(y - 2.0)**2)  * exp(-0.25*(z - 2.0)**2)
        exp2 = 0.75 * exp( -1.0/49.0 * (x + 1.0)**2) * exp(-0.10*(y + 1.0))     * exp(-0.10*(z + 1.0))
        exp3 = 0.50 * exp( -0.25 * (x - 7.0)**2)     * exp(-0.25*(y - 3.0)**2)  * exp(-0.25*(z - 5.0)**2)
        exp4 = 0.20 * exp( -(x - 4.0)**2)            * exp(-(y - 7.0)**2)       * exp(-(z - 5.0)**2)

        f2 = exp1 + exp2 + exp3 - exp4

    end function f2

    !> Test function 3
    pure function f3(r)
        !> Mesh in which to compute
        real(r64), intent(in) :: r(:,:)
        !> Result
        real(r64), allocatable :: f3(:)
        !> Elements
        real(r64), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f3 = (1.0_r64 + tanh(-alpha*(x+y-z)))/alpha

    end function f3

    !> Sgn function
    pure function sgn(value)
        
        real(r64), intent(in) :: value(:)
        real(r64), allocatable :: sgn(:)
        integer(i64) :: i
        
        allocate(sgn,source=value)

        do i = 1, size(value,1)
            if (sgn(i) .lt. 0.0_r64) then
                sgn(i) = -1.0_r64
            else if (sgn(i) .eq. 0.0_r64) then
                sgn(i) =  0.0_r64
            else 
                sgn(i) =  1.0_r64
            end if
        end do

    end function sgn

    !> Test function 4
    pure function f4(r)
        !> Mesh in which to compute
        real(r64), intent(in) :: r(:,:)
        !> Result
        real(r64), allocatable :: f4(:)
        !> Elements
        real(r64), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f4 = (1.0_r64 - sgn(x+y-z))/alpha
        
    end function f4

    !> Test function 4
    pure function f5(r)
        !> Mesh in which to compute
        real(r64), intent(in) :: r(:,:)
        !> Result
        real(r64), allocatable :: f5(:)
        !> Elements
        real(r64), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f5 = (1.0_r64 - sgn(pi*x + y))/alpha
        
    end function f5
    

end program test_quadrature_lebedev