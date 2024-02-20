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
!> Contains procedures to test averages in 2D cases (MoS2)

!> This module contains the procedures for the test of averaging of the dielectric matrix in 2D cases
program test_2d

    use idiel_constants, only: r64, i64, pi, twopi, zzero
    use idiel, only: idiel_t

    implicit none

    ! GW data
    integer, parameter :: mbsize = 438 ! Mixed basis size
    integer, parameter :: nomega = 12  ! Number of frequencies
    integer(i64) :: iom

    ! Mesh data
    integer(i64), parameter :: ngrid(3) = [2,2,1]

    ! Silicon crystal data
    integer(i64), parameter :: natoms = 3_i64
    real(r64) :: a(3), b(3), c(3), lattice(3,3)
    real(r64), allocatable :: redpos(:,:)
    integer(i64), allocatable :: types(:)

    ! Dielectric data
    complex(r64), allocatable :: head(:,:,:), wingL(:,:,:)
    complex(r64), allocatable :: wingU(:,:,:), body(:,:,:)

    ! Reference data to compare
    real(r64) :: head_ref(4) = [ -0.27165343E+03_r64, -0.97552958E+02_r64, -0.49406608E+02_r64, -0.46006706E+01_r64]
    real(r64) :: rdiff
    real(r64), parameter    :: tolerance = 0.05_r64

    ! Computation object
    type(idiel_t) :: inv_diel

    ! In case that no SPG we set only Identity (i.e. no symmetry)
    integer(i64) :: nsym = 1
    real(r64)    :: crot(3,3,1)

    ! Lattice vectors
    a(:) = [0.00000000_r64, 6.020289060_r64, 0.00000000_r64]
    b(:) = [5.213723260_r64, 3.010144530_r64, 0.00000000_r64]
    c(:) = [0.0_r64,  0.0_r64, 20.0_r64]

    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [1_i64, 2_i64, 2_i64]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.00_r64, 0.00_r64, 0.00_r64]
    redpos(:,2) = [2.0_r64/3.0_r64 , 2.0_r64/3.0_r64,  0.14764992_r64]
    redpos(:,3) = [2.0_r64/3.0_r64 , 2.0_r64/3.0_r64, -0.14764992_r64]

    ! Init common objects
    write(*,*) 'Init'
#ifdef USE_SPGLIB
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i64)
#else 
    crot = 0.0_r64
    crot(1,1,1) = 1.0_r64
    crot(2,3,1) = 1.0_r64
    crot(3,3,1) = 1.0_r64
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i64, nsym, crot)
#endif

    ! Read the dielectric data from previous G0W0 run
    call load_from_file('head.2d.dat',  head)
    call load_from_file('wings.2d.dat', wingL, wingU)
    call load_from_file('body.2d.dat',  body)

    ! Do the average
    write(*,*) '[TEST: idiel_t (2d)]'
    do iom = 1, size(head,3)
        ! Invert body
        write(*,*) 'Invert'
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        write(*,*) 'Set blocks'
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        write(*,*) 'Compute'
        call inv_diel%compute_anisotropic_avg_scrcoulomb_2d(.false.)

        ! Check the head
        rdiff = abs(inv_diel%idiel_head - head_ref(iom))/abs(head_ref(iom))
        write(*,'(A,I2,A,e20.13)')  '  * Regression (HEAD,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
             write(*,*)  '[TEST: idiel_t (2d) (HEAD,',iom,'): PASSED]'
        else
             write(*,*)  '[TEST: idiel_t (2d) (HEAD,',iom,'): FAILED]'
             stop 1
        end if

    end do

    stop 0

contains 

    !> Process a file and load to anray
    !> @param[in] fname - the file name   
    !> @param[in] data_shape - the shape of the data to reads
    !> @param[out] data - the data
    subroutine load_from_file(fname, data, data2)

        character(len=*), intent(in) :: fname
        complex(r64), allocatable, intent(out) :: data(:,:,:)
        complex(r64), allocatable, optional, intent(out) :: data2(:,:,:)

        integer :: fin
        real(r64), allocatable  :: dreal(:,:,:), dimag(:,:,:)
        integer(i64) :: data_shape(3)
        
        open(file=fname, newunit=fin, status='old', action='read')

        read(fin, *) data_shape

        if (present(data2)) then 
            allocate(data(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dreal(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dimag, mold=dreal)
            allocate(data2, mold=data)
            read(fin, *) dreal
            read(fin, *) dimag
            data = cmplx(dreal, dimag, r64) 
            read(fin, *) dreal
            read(fin, *) dimag
            data2 = cmplx(dreal, dimag, r64) 
        else
            allocate(data(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dreal(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dimag, mold=dreal)
            read(fin, *) dreal
            read(fin, *) dimag
            data = cmplx(dreal, dimag, r64) 
        end if

        close(fin)

    end subroutine load_from_file

end program test_2d
