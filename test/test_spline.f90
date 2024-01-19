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
!> Tests the computation of circular 
program test_spline_1d

    use idiel_constants, only: r64, i64, twopi, pi
    use idiel_cubic_spline

    implicit none
    
    integer(i64), parameter :: n = 50_i64
    complex(r64) :: y(n), yref(10*n), yint
    real(r64)    :: x(n), xdenser(10*n)
    integer(i64) ::  i
    real(r64) :: dx, dx2
    type(cubic_spline_t) :: cs

    dx  = twopi / (n-1)
    dx2 = twopi / (10*n-1) 

    do i = 1, n
        x(i) = (i - 1) * dx
        y(i) = cmplx(cos(x(i))**2 * x(i) , 0.0, r64) 
    end do

    do i = 1, 10*n
        xdenser(i) = (i - 1) * dx2
        yref(i) = cmplx(cos(xdenser(i))**2 * xdenser(i) , 0.0, r64) 
    end do

    call cs%initialize(x, y)
    
    do i = 1, 10*n
        yint = cs%interpolate(xdenser(i))
        write(*,*) i, xdenser(i), real(yref(i)), real(yint)
    end do

    write(*,*) 2.35718, cs%integrate(0.5_r64,pi)

end program test_spline_1d