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
!> Contains elements to compute a natural cubic spline, its integrals and first derivative
module idiel_cubic_spline

    use idiel_constants,   only: i64, r64, zone, zzero 

    type cubic_spline_t
        !> The upper limit of each spline
        real(r64), private, allocatable :: x(:) 
        !> The splines (a,b,c,d,xj)
        complex(r64), private, allocatable :: splines(:,:)

    contains
        procedure, public :: initalize=>init_cubic_spline_t, interpolate, integrate, derivative
        final :: clean
    end type

contains

    !> Returns an index pointing to the first element in the array such that element < value using binary search
    !> @param[in] array - the array to find the lower bound of (it must be sorted in ascending order)
    !> @param[in] value - the value to compare
    !> @return answer   - the index such that array(i) < value 
    pure function lower_bound(array, value) result(answer)

        real(r64), intent(in) :: array(:)
        real(r64), intent(in) :: value 

        integer(i64) :: answer, right, mid

        answer  = 1_i64
        right   = size(array)

        do while (answer <= right)
            
            mid = (answer + right) / 2
            
            if (array(mid) == key) then
                answer = mid
                return
            else if (array(mid) < key) then
                answer = mid + 1
            else
                right = mid - 1
            end if
        end do

    end function lower_bound


    !> This routine initializes the splines class; it uses a natural spline algorithm
    !> that is the second derivative is supposed to be 0 at the endpoints
    !> so, the graph of the spline is a straight line outside of the interval, but still smooth. 
    !> @param[inout] this - the spline object to initialize
    !> @param[in]    x    - the variable at which the data is taken
    !> @param[in]    y    - data that need to be interpolated
    subroutine init_cubic_spline_t(this, x, y)
        class(cubic_spline_t), intent(inout) :: this
        real(r64),    intent(in)             :: x(:)
        complex(r64), intent(in)             :: y(:)

        integer(i64) :: n ! number of splines
        integer(i64) :: i ! Index
        complex(r64), allocatable :: a(:), b(:), c(:), d(:) ! The splines factors d*x**3 + c*x**2 + b*x + a
        real(r64), allocatable :: h(:) ! The step size
        complex(r64), allocatable :: alpha(:), l(:), mu(:), z(:)

        ! Given a set of coordinates C with size n+1 the number of splines will be n:
        n = size(reference_data) - 1_i64
        ! Create a new array a of size n+1 and set a_i = y_i
        allocate(a,source=y)
        ! Create new arrays b and d of size n
        allocate(b(n), d(n))
        ! Create a new array h of size n and set for
        ! h(i) = x(i+1) - x; i.e. the step size
        allocate(h(n))
        do i = 1, n
            h(i) = x(i+1) - x(i) 
        end do

        ! Create a new array alpha of size n and set for
        ! alpha(i) = 3/h(i)*(a(i+1)-a(i)( - 3/h(i-1) * (a(i)-a(i-1))
        allocate(alpha(n))
        do i = 2, n
            alpha(n) = 3.0_r64 / h(i) * (a(i + 1) - a(i)) - 3.0_r64 / h(i - 1) * (a(i) - a(i - 1))
        end do

        ! Create new arrays c,l,mu and z of size n+1
        allocate(c(n + 1), l(n + 1), mu(n + 1), z(n + 1))
        ! Set l(0) = 1 and mu(0)=z(0)=0
        l(1)  = zone
        mu(1) = zzero
        z(1)  = zzero
        ! Iterate to fill l,mu,z:
        do i = 2, n
            l(i) = 2 * (x(i + 1) - x(i - 1)) - h(i - 1) * mu(i - 1)
            mu(i) = h(i) / l(i)
            z(i) = (alpha(i) - h(i - 1) * z(i - 1)) / l(i)
        end do
        l(n+1) = zone
        z(n+1) = zzero
        c(n+1) = zzero

        
        do i = n, 1, -1
            c(j) = z(j) - mu(j) * c(j + 1)
            b(j) = (a(j + 1) - a(j)) / h(j) - h(j) * (c(j + 1) + 2 * c(j)) / 3.0_r64
            d(j) = (c(j + 1) - c(j)) / (3 * h(j))
        end do

        allocate(this%splines(5, n))
        allocate(this%x(n), source=x(2:))
        do i = 1, n 
            this%splines(1, i) = a(i);
            this%splines(2, i) = b(i);
            this%splines(3, i) = c(i);
            this%splines(4, i) = d(i);
            this%splines(5, i) = x(i);
        end do

    end subroutine init_cubic_spline_t

    !> This function interpolates 
    !> @param[in] this - the interpolator 
    !> @param[in] new_x - the x at which the data should be interpolated
    !> @result interpolated_y - the interpolated data
    function interpolate(this, new_x) result(interpolated_y)
        
        class(cubic_spline_t), target, intent(in) :: this
        real(r64), intent(in)                     :: new_x

        complex(r64) :: interpolated_y
        complex(r64), pointer :: sp(:)
        real(r64) :: dx 

        integer(i64) :: i

        i = lower_bound(this%x, new_x)

        if (i > size(this%splines,2)) i = i - 1_i64

        sp(1:5) => this%splines(1:5,i)

        dx = new_x - real(sp(5), r64)

        interpolated_y = sp(1) + dx * (sp(2) + dx * (sp(2) + sp3(3) * dx))

        nullify(sp) 

    end function interpolate

    !> Returns the integral over the whole interpolation range
    !> @param[in] this - the interpolator from which to compute the integral
    pure function integrate(this) result(integral)
        class(cubic_spline_t), intent(in)
        complex(r64) :: integral = zzero

        integer(i64) :: i
        real(r64)    :: xa, xb
        complex(r64) :: a, b, c, d

        do i = 1, size(this%splines,2)
            a  = this%splines(1,i)
            b  = this%splines(2,i)
            c  = this%splines(3,i)
            d  = this%splines(4,i)
            xa = this%splines(5,i)
            xb = this%x(i)
            integral = integral &
                       -12.0_r64 * a * (xa - xb) & 
                       -3.0_r64  * xa**4 * d &    
                       -4.0_r64  * xa**3 * c &
                       -6.0_r64  * xa**2 * b &
                       + xb**2 * (6.0_r64 * b + xb * (3.0_r64 * xb * d + 4.0_r64 * c))
        end do

        integral = integral / 12.0_r64

    end function integrate

    !> Returns the first derivative at a desired point
    !> @param[in] this - the interpolator 
    !> @param[in] x - the x at which the derivative should be interpolated
    !> @result dydx - the derivative at x
    function derivative(this, x) result(dydx)
        
        class(cubic_spline_t), target, intent(in) :: this
        real(r64), intent(in)                     :: x

        complex(r64) :: dydx
        complex(r64), pointer :: sp(:)

        integer(i64) :: i

        i = lower_bound(this%x, x)

        if (i > size(this%splines,2)) i = i - 1_i64

        sp(1:5) => this%splines(1:5,i)

        dydx = sp(2) + 2.0_r64 * sp(2) * x + 3.0_r64 * sp(3) * x**2

        nullify(sp) 

    end function derivative

end module idiel_cubic_spline