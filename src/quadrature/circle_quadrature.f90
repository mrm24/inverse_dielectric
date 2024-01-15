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
!> Contains elements to compute Gauss-Legendre grid of 2500th order for a unitary
!> circle
module idiel_circle_quadrature
        
    use idiel_constants,   only: i64, r64, pi, twopi
    use iso_c_binding

    implicit none
    
    private 
    public :: compute_angular_mesh_gauss_legendre

    !> Interface to CXX boost implementation of gauss legendre quadrature
    subroutine compute_gauss_legendre_cxx(n, x, w) bind(C,name="double_compute_gauss_legendre")
        import
        integer(C_int), value  :: n
        type(C_ptr), value     :: x
        type(C_ptr), value     :: w
    end subroutine compute_gauss_legendre_cxx

contains 

    !> Polar to cartesian coordinates
    !> @param[in] rphi - polar mesh (r,phi)
    !> @result    xyz - the Cartesian coordinates 
    pure function polar_2_cartesian(rphi) result(xyz)
        !> The angular mesh (theta,phi)
        real(r64), intent(in) :: rphi(:,:)
        !> The Cartesian mesh
        real(r64), allocatable :: xyz(:,:)

        allocate(xyz(size(rphi,1),3))

        xyz(:,1) = rphi(:,1) * cos(rphi(:,2))
        xyz(:,2) = rphi(:,1) * sin(rphi(:,2))
        xyz(:,3) = 0.0_r64

    end function polar_2_cartesian

    !> Cartesian to polar coordinates
    !> @param[in] xyz - the Cartesian coordinates
    !> @result    rphi - polar mesh (r,phi) 
    pure function cartesian_2_polar(xyz) result(rphi)
        !> The Cartesian mesh
        real(r64), intent(in)  :: xyz(:,:)
        !> The angular mesh (theta,phi)
        real(r64), allocatable :: rphi(:,:)

        allocate(rphi(size(xyz,1),2), sourece=0.0_r64)

        rphi(:,1) = hypot(xyz(:,1),xyz(:,2))
        rphi(:,2) = atan2(xyz(:,2),xyz(:,1))

    end function cartesian_2_polar


    !> Compute the gauss legendre mesh
    !> @param[in]   n - the number of points 
    !> @param[in]   a - the lower bound of the integral
    !> @param[in]   b - the upper bound of the integral
    !> @param[out]  x - the abscisa points
    !> @param[out]  w - weights
    subroutine compute_gauss_legendre_f90(n, a, b, x, w)

        integer(i64), intent(in) :: n
        real(r64), intent(in) :: a
        real(r64), intent(in) :: b
        real(r64), allocatable, target, intent(out) :: x(:)
        real(r64), allocatable, target, intent(out) :: w(:)

        real(r64) :: dx, shift
         
        allocate(x(n),w(n))

        call compute_gauss_legendre_cxx(int(n,C_int),C_loc(x),C_loc(w))

        ! Remap
        dx    = 0.5_r64 * ( b - a )
        shift = 0.5_r64 * ( a + b )
        x(:)  = dx * x(:) + shift
        w(:)  = dx * w(:)

    end subroutine compute_gauss_legendre_f90


    !> Compute an angular mesh Gauss-Legendre for the angle
    !> @param[out]  rphi - the angles (r,phi)
    !> @param[out]  w    - weights
    !> @param[out]  xyz  - the mesh in cartesian coordinates
    subroutine compute_angular_mesh_gauss_legendre(rphi, w, xyz)
        
        real(r64), allocatable, intent(out) :: rphi(:,:)
        real(r64), allocatable, intent(out) :: w(:)
        real(r64), allocatable, intent(out) :: xyz(:,:)

        ! Phi mesh
        real(r64), allocatable :: w_phi(:), x_phi(:)

        ! Mesh size
        integer(i64), parameter :: mesh_size = 2500_i64 
        ! Idx
        integer(i64) :: i, idx 

        ! Build theta mesh
        call compute_gauss_legendre_f90(mesh_size, 0.0_r64, twopi, x_theta, w_theta)

        allocate(w(mesh_size))
        allocate(ang(mesh_size,2), source=1.0_r64)

        ! Save in proper format
        rphi(:,2) = x_phi(:)
        rphi(:,1) = w_phi(:)

        ! Go to cartesian
        xyz =  polar_2_cartesian(rphi)

    end subroutine compute_angular_mesh_gauss_legendre


    

end module idiel_circle_quadrature