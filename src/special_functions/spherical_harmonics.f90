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
!> Contains elements to compute harmonic spherics and related elements

!> This module contains procedures to compute spherical harmonics and functions expansion in terms of such a basis
module idiel_spherical_harmonics
   
   use idiel_constants,   only: i64, r64, fourpi, y00, zzero, zone

   implicit none 

   private
   public :: sph_harm, sph_harm_expansion

contains 
    
   !>   Generates a sequence of spherical harmonics, including the Condon-Shortley
   !>   phase, evaluated at angles $(\theta,\phi)$ for $0<l<l_{\rm max}$. The values
   !>   are returned in a packed array {\tt ylm} indexed with $j=l(l+1)+m+1$. The
   !>   algorithm of Masters and Richards-Dinger is used, {\it Geophys. J. Int.}
   !>   {\bf 135}, 307 (1998). This routine is numerically stable and accurate to
   !>   near machine precision for $l\le 50$.
   !>   @param[in] lmax - the maximum momentum for which the spherical harmonic are computed
   !>   @param[in] ang  - the angular points for which the spherical harmonics are computed 
   !>   @param[out] ylm - the spherical harmonics at ang points
   subroutine sph_harm(lmax, ang, ylm)
   
      integer(i64), intent(in)  :: lmax
      real(r64),    intent(in)  :: ang(:,:)
      complex(r64), allocatable, intent(out) :: ylm(:,:)
      
      ! local variables
      integer(i64) :: l, m, lm1, lm2, npoints
      real(r64), allocatable :: sn(:), cs(:), dx(:), cumul(:), t1(:)
      complex(r64), allocatable :: z(:,:)
      real(r64), allocatable :: x(:,:)
      
      npoints = size(ang, 1)
   
      allocate(sn(npoints), cs(npoints), dx(npoints), cumul(npoints), t1(npoints))
      allocate(x(npoints,0:lmax))
      allocate(z(npoints,lmax))
      
      if ((lmax < 0) .or. (lmax > 50)) then
          error stop "sph_harm: lmax out of range, max value for 0 =< lmax <= 50"
      end if
      
      allocate(ylm(npoints,(lmax+1)**2), source = y00)
      
      !$omp parallel default(shared) private(l, m, lm1, lm2, sn, cs, dx, cumul, t1, z, x)
      !$omp do schedule(dynamic)
      do l = 1, lmax
         x = 0.0_r64
         if (mod(l, 2) == 0) then
            x(:,0:l) = 1.0_r64
         else
            x(:,0:l) = -1.0_r64
         endif
         
         sn = sin(ang(:, 1))
         cs = cos(ang(:, 1))
         
         ! phase factors exp(i*m*phi)
         do m = 1, l
            t1 = real(m, kind=r64) * ang(:, 2)
            z(:,m) = cmplx(cos(t1), sin(t1), kind=r64)
         end do
         
         dx = 0.0_r64
         do m = l, 1, -1
            t1 = sqrt(real((l+m)*(l-m+1), kind=r64))
            x(:,m-1) = -(sn*dx+real(2*m, kind=r64) * cs * x(:,m)) / t1
            dx = sn * x(:,m) * t1
         end do
         
         t1 = sn
         cumul = 0.0_r64
         do m = 1, l
            x(:,m) = t1 * x(:,m)
            cumul = cumul + x(:,m)**2
            t1 = t1 * sn
         end do
         
         cumul = 2.0_r64 * cumul + x(:,0)**2
         t1 = sqrt(real(2*l+1, kind=r64)/(fourpi*cumul))
         lm1 = l * (l+1) + 1
         lm2 = lm1
         
         ylm(:, lm1) = t1 * x(:,0)
         do m = 1, l
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            ylm(:, lm1) = t1 * x(:,m) * z(:,m)
            ylm(:, lm2) = conjg(ylm(:, lm1))
            if (mod(m, 2) /= 0) ylm(:, lm2) = -ylm(:, lm2)
         end do
      end do
      !$omp end do
      !$omp end parallel 
   
      deallocate(sn, cs, dx, cumul, t1, x)
   
   end subroutine sph_harm
   
   !> It expands the function f in terms of spherical harmonics 
   !> Notice that this function is intended to be used with a Lebedev grid of order N so the
   !> expansion coefficients should be accurate up to half of the order of the grid
   !> @param[in] n    - number of spherical harmonics
   !> @param[in] f    - the function to expand
   !> @param[in] weights - the mesh points weights
   !> @param[in] ylm     - the spherical harmonics in the mesh
   !> @param[out] clm     - the expansion coefficients
   subroutine sph_harm_expansion(n, f, weights, ylm, clm)
#ifdef USE_GPU
      !$omp declare target
#endif 
      integer(i64), intent(in)  :: n
      complex(r64), intent(in)  :: f(:)
      real(r64),    intent(in)  :: weights(:)
      complex(r64), intent(in)  :: ylm(:,:)
      complex(r64), intent(out) :: clm(n)

      complex(r64), allocatable :: fw(:)
      integer(i64) :: i, j

      ! Multiply the function with the weights so the matrix-vector products solves the integral for all the (l,m)
      allocate(fw, source=f)
      fw(:) = fw(:) * weights(:)

      ! This is a small calculation so it is done with a hand-made definition of matmul so the GPU kernel can be built
      ! clm = Ylm**H \cdot f
      do i = 1, n
          clm(i) = dot_product(ylm(:,i),fw(:))
      end do

      deallocate(fw)

   end subroutine sph_harm_expansion

end module idiel_spherical_harmonics
