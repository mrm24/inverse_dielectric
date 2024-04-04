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
!> Tests the computation of 1D splines
program test_spline_1d

    use idiel_constants, only: aip, i32, twopi, pi
    use idiel_cubic_spline

    implicit none
    
    integer(i32), parameter :: n = 60_i32
    real(aip), parameter :: tolerance  = 5.0e-3_aip
    real(aip), parameter :: rtolerance = 5.0e-5_aip
    real(aip), parameter :: ana_integral = 2.35718_aip ! from 0.5 to pi
    real(aip), parameter :: ana_derivative = -0.550858_aip ! for 0.3*pi
    complex(aip) :: y(n), yref(10*n), yint
    real(aip)    :: x(n), xdenser(10*n)
    integer(i32) ::  i
    real(aip) :: dx, dx2, diff, rdiff, ics, dcs
    !! Splines
    real(aip), allocatable    :: xlim(:)
    complex(aip), allocatable :: splines(:,:), integrals(:)

    ! Print info
    write(*,*) '[TEST : cubic_spline_t]' 
    write(*,*) ' * kind : regression test against precomputed data'
    write(*,*) ' * Test function : y(x) = cos(x)**2 * x '  

    dx  = twopi / (n-1)
    dx2 = twopi / (10*n-1) 

    do i = 1, n
        x(i) = (i - 1) * dx
        y(i) = cmplx(cos(x(i))**2 * x(i) , 0.0, aip) 
    end do

    do i = 1, 10*n
        xdenser(i) = (i - 1) * dx2
        yref(i) = cmplx(cos(xdenser(i))**2 * xdenser(i) , 0.0, aip) 
    end do

    call init_cubic_spline_t(xlim, splines, integrals, x, y)
    
    write(*,'(A, e20.13)')  '  * Compare vs exact function:  cubic_spline_t%interpolate'
    do i = 1, 9*n
        yint = interpolate(xlim, splines, xdenser(i))
        if (real(yref(i)) == 0) cycle
        diff = abs(real(yint)-real(yref(i)))
        write(*,*)  xdenser(i), real(yref(i)), real(yint)
        if ( diff > tolerance) then
           write(*,*)  '[TEST : cubic_spline_t%interpolate: FAILED] , diff : ', diff
           stop 1
        else
           write(*,*)  '[TEST : cubic_spline_t%interpolate: PASSED] , diff : ', diff
        end if
    end do

    write(*,'(A, e20.13)')  '  * Compare vs exact result:  cubic_spline_t%derivative'
    dcs = real(derivative(xlim, splines,0.3_aip*real(pi,aip)))
    rdiff = abs(dcs - ana_derivative)/abs(ana_derivative)
    if ( rdiff > 100*rtolerance) then
        write(*,*)  '[TEST : cubic_spline_t%derivative: FAILED] , rdiff : ', rdiff
        stop 1
     else
        write(*,*)  '[TEST : cubic_spline_t%derivative: PASSED] , rdiff : ', rdiff
     end if

    write(*,'(A, e20.13)')  '  * Compare vs exact result:  cubic_spline_t%integrate'
    ics = real(integrate(xlim, splines, integrals,0.5_aip,real(pi,aip)))
    rdiff = abs(ics - ana_integral)/abs(ana_integral)
    if ( rdiff > rtolerance) then
        write(*,*)  '[TEST : cubic_spline_t%integrate: FAILED] , rdiff : ', rdiff
        stop 1
     else 
        write(*,*)  '[TEST : cubic_spline_t%integrate: PASSED] , rdiff : ', rdiff
     end if 

     stop 0

end program test_spline_1d
