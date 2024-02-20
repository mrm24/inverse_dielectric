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
    
    use idiel_constants,         only: i64, r64, zzero
    use idiel_linalg,            only: inverse_complex_LU
#ifdef USE_GPU
    use idiel_gpu_magma_t, only: linalg_world_t
#else
    use idiel_cpu_magma_t, only: linalg_world_t
#endif
    
    implicit none

    integer(i64), parameter :: n = 3000_i64
    complex(r64), allocatable :: A(:,:), inverse_A(:,:)
    type(linalg_world_t) :: world
    integer(i64) :: i
    real(r64) :: rdiff
    real(r64), parameter    :: tolerance = 1.0e-12_r64

    allocate(A(n,n), source=zzero)

    do i = 1, n
        A(i,i) = cmplx(i, 0, r64)
    end do
    call world%init()
    call inverse_complex_LU(A, inverse_A, world)
    call world%finish()
    write(*,*) '[TEST : inverse_complex_LU]' 
    do i = 1, n
        write(*,*) 1.0_r64 / A(i,i), inverse_A(i,i)
        rdiff = abs(1.0_r64 / A(i,i) - inverse_A(i,i))/ abs(1.0_r64 / A(i,i))
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
