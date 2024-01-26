! Copyright 2023 EXCITING developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   htang://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied. See the License for the specific language governing
! permissions and limitations under the License.

!> @file
!> Contains elements to compute circular harmonics and related elements

!> This module contains procedures to compute circular harmonics and functions expansion in terms of such a basis
module idiel_circular_harmonics
   
   use idiel_constants,   only: i64, r64, pi, twopi, zone, zzero

   implicit none 

   private
   public :: circ_harm, circ_harm_expansion


contains 
    
   !>   Generates a sequence of circular harmonics
   !>   @param[in] lmax - the maximum momentum for which the circular harmonic are computed
   !>   @param[in] ang  - the angular points for which the circular harmonics are computed 
   !>   @param[out] blm - the spherical harmonics at ang points
   subroutine circ_harm(lmax, ang, blm)
   
      integer(i64), intent(in)  :: lmax
      real(r64),    intent(in)  :: ang(:)
      complex(r64), allocatable, intent(out) :: blm(:,:)
      
      ! local variables
      integer(i64) :: l, m, npoints
      
      npoints = size(ang)
      
      allocate(blm(npoints,(2*lmax+1)), source = cmplx(1.0/sqrt(pi), 0.0, r64))
      
      blm(:,1) = cmplx(1.0/sqrt(twopi), 0.0, r64)

      !$omp parallel default(shared) private(l, m)
      !$omp do schedule(dynamic)
      do l = 1, lmax
         blm(:, 2*l    ) = blm(:, 2*l    ) * cmplx(sin(l * ang(:)), 0.0, r64)
         blm(:, 2*l + 1) = blm(:, 2*l + 1) * cmplx(cos(l * ang(:)), 0.0, r64)
      end do
      !$omp end do
      !$omp end parallel 
   
   end subroutine circ_harm
    
   !> It expands the function f in terms of spherical harmonics (up to l = 10)
   !> Notice that this function is intended to be used with a Gauss-Legendre grid of size 55 which integrates all harmonics up to 10
   !> expansion coefficients should be accurate up to l = 10
   !> @param[in] lmax     - the spherical expansion order
   !> @param[in] f        - the function to expand (the values for different radius are provided in the radii)
   !> @param[in] weights  - the mesh points weights
   !> @param[in] blm      - the spherical harmonics in the mesh
   !> @param[out] clm     - the expansion coefficients
   subroutine circ_harm_expansion(lmax, f, weights, blm, clm)

      integer(i64), intent(in)   :: lmax
      complex(r64), intent(in)   :: f(:)
      real(r64),    intent(in)   :: weights(:)
      complex(r64), intent(in)   :: blm(:,:)
      complex(r64), intent(out)  :: clm(2*lmax+1)

      integer(i64)              :: i
      integer(i64), parameter   :: nr = 55_i64
      complex(r64), allocatable :: fw(:)
      external                  :: zgemv

      ! Checks
      if (lmax > 10) then
         error stop "Error(circ_harm_expansion): maximum order of expansion is 10"
      end if

      if (size(f,1) /= nr .or. size(blm,1) /= nr .or. size(weights) /= nr) then
         error stop "Error(circ_harm_expansion): the expected mesh is for Gauss-Legendre grid of order 21"
      end if

      ! Multiply the function with the weights so the matrix-vector products solves the integral for all the (l,m)
      allocate(fw, source=f)
      fw(:) = fw(:) * weights(:)
      ! This is a small calculation so it is done in the CPU (Notice that we use only T not H since blm imaginary part is zero)
      ! clm = blm**T \cdot f
      call zgemv('T', size(blm,1), size(blm,2), zone, blm, size(blm,1), fw, 1, zzero, clm, 1)

   end subroutine circ_harm_expansion

end module idiel_circular_harmonics
