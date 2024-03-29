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
!> Tests the circular expansion and quadrature
program test_circular_expansion

    use idiel_constants, only: r64, i64, twopi, pi, fourpi, zzero
    use idiel_circular_harmonics
    use idiel_circle_quadrature 

    implicit none

    integer(i64), parameter :: msize = 75_i64
    integer(i64), parameter :: lmax  = 20_i64
    real(r64),    parameter :: tolerance = 1.0e-12_r64
    real(r64), allocatable :: rphi(:,:), w(:), xyz(:,:)
    complex(r64), allocatable :: blm(:,:)
    complex(r64), allocatable :: f2expand(:), fexact(:), fint(:)
    complex(r64) :: clm(2*lmax+1)
    real(r64) :: integral, rdiff
    integer(i64) :: i


    write(*,*) '[TEST : idiel_circular_harmonics and idiel_circle_quadrature]' 
    write(*,*)  '  * Check exact integration of circular basis'
    ! We compute the angular mesh
    call compute_angular_mesh_gauss_legendre(msize, rphi, w, xyz)

    ! We compute the circular harmonics up to 22
    call circ_harm(msize, rphi(:,2), blm)

    ! Check that the quadrature integrates the circular harmonics
    integral = real(sum(w * blm(:,1)))
    if ( abs(integral-2.50662827463100_r64)/2.50662827463100_r64 <= tolerance) then
        write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: PASSED]'
    else   
        write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: FAILED]'
        stop 1
    end if

    do i = 2, 2*lmax+1
        integral = real(sum(w * blm(:,i)))
        if ( integral <= tolerance) then
            write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: PASSED]'
        else   
            write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: FAILED]'
            stop 1
        end if
    end do

    write(*,*)  '[TEST : circ_harm_expansion]'
    f2expand = f1(xyz)
    call circ_harm_expansion(10_i64, f2expand, w, blm, clm)
    deallocate(rphi, w, xyz, blm)
    
    call compute_angular_mesh_gauss_legendre(10*msize, rphi, w, xyz)
    call circ_harm(msize, rphi(:,2), blm)
    fexact = f1(xyz)
    allocate(fint, mold=fexact)
    fint = zzero
    
    ! Compute the expansion
    do i = 1, 2*lmax+1
        fint(:) = fint(:) + clm(i) * blm(:,i) 
    end do

    ! Compute the function
    do i = 1, 10*msize
        rdiff = abs(real(fint(i))-real(fexact(i))) / abs(real(fexact(i)))
        if ( rdiff <= tolerance) then
            write(*,*)  '[TEST : circ_harm_expansion: PASSED]'
        else   
            write(*,*)  '[TEST : circ_harm_expansion: FAILED]'
            stop 1
        end if
    end do 

    stop 0

contains

    pure function f1(r)
        !> Mesh in which to compute
        real(r64), intent(in)  :: r(:,:)
        !> Result
        complex(r64), allocatable :: f1(:)
        !> Elements
        real(r64), allocatable :: x(:), y(:)

        x = r(:,1)
        y = r(:,2)
        
        allocate(f1(size(r,1)))

        f1 = cmplx(1.0_r64 + x + y**2 + y * x**2 + x**4 + y**5 + (x*y)**3, 0.0_r64, r64)

    end function f1

end program test_circular_expansion
