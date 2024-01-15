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
!> This module contains parameters and constans
!>
module idiel_constants
    
    use iso_fortran_env, only: real64, int32, int64
    
    implicit none
    
    public 
    
    !> Real double precision
    integer, parameter :: r64 = real64
    !> Integer single precision
    integer, parameter :: i32 = int32
    !> Integer double precision
    integer, parameter :: i64 = int64
    
    !> Pi 
    real(r64), parameter :: pi    = 4.0_r64 * atan(1.0_r64)
    !> Two pi
    real(r64), parameter :: twopi = 2.0_r64 * pi
    !> Four pi
    real(r64), parameter :: fourpi = 4.0_r64 * pi
    !> 1.0/3.0
    real(r64), parameter :: onethird = 1.0_r64 / 3.0_r64
    !> The Y_{0}^{0} spherical harmonic
    complex(r64), parameter :: y00 = cmplx(0.5_r64 * sqrt(1.0/pi),  0.0_r64, kind=r64)
    
    !> Imaginary unit
    complex(r64), parameter :: iunit = cmplx(0.0_r64, 1.0_r64, kind=r64)
    !> Complex one
    complex(r64), parameter :: zone  = cmplx(1.0_r64, 0.0_r64, kind=r64)
    !> Complex zero
    complex(r64), parameter :: zzero = cmplx(0.0_r64, 0.0_r64, kind=r64)
    
end module idiel_constants
