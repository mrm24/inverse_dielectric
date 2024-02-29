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
!> Test the inversion using either GPU or CPU 
program test_inverse
    
    use idiel_constants,         only: i32, aip, zzero
    use idiel_linalg,            only: inverse_complex_LU
    use idiel_gpu_world_t,       only: gpu_world_t
    
    implicit none

    integer(i32), parameter :: n = 3000_i32
    complex(aip), allocatable :: A(:,:), inverse_A(:,:)
    type(gpu_world_t) :: world
    integer(i32) :: i
    real(aip) :: rdiff
    real(aip), parameter    :: tolerance = 1.0e-12_aip

    allocate(A(n,n), source=zzero)

    do i = 1, n
        A(i,i) = cmplx(i, 0, aip)
    end do
    call world%init()
    call inverse_complex_LU(A, inverse_A, world)
    call world%finish()
    write(*,*) '[TEST : inverse_complex_LU]' 
    do i = 1, n
        write(*,*) 1.0_aip / A(i,i), inverse_A(i,i)
        rdiff = abs(1.0_aip / A(i,i) - inverse_A(i,i))/ abs(1.0_aip / A(i,i))
        write(*,'(A,I0,A, E16.7)')  '  * Regression (inverse,', i ,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST : inverse_complex_LU (inverse,', i,'): PASSED]'
        else
            write(*,*)  '[TEST : inverse_complex_LU (inverse,', i,'): FAILED]'
            stop 1
        end if
    end do

    stop 0

end program test_inverse
