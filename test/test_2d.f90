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
    integer(i64), parameter :: ngrid(3) = [6,6,6]

    ! Silicon crystal data
    integer(i64), parameter :: natoms = 3_i64
    real(r64) :: a(3), b(3), c(3), lattice(3,3)
    real(r64), allocatable :: redpos(:,:)
    integer(i64), allocatable :: types(:)

    ! Dielectric data
    complex(r64), allocatable :: head(:,:,:), wingL(:,:,:)
    complex(r64), allocatable :: wingU(:,:,:), Binv(:,:,:)

    ! Reference data to compare
    real(r64) :: head_ref(12) = [ 6.442031104885004E-002, 0.132222402626860, 0.316390106319300, 0.516127938510703, &
                                  0.645769722906067, 0.707300869613546, 0.733254256631726, 0.785227466265463, &
                                  0.864527818748701, 0.941940576487719, 0.987846211386285, 1.00132129472395 ]
    real(r64) :: rdiff
    real(r64), parameter    :: tolerance = 1.0e+7_r64

    ! Computation object
    type(idiel_t) :: inv_diel

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
    redpos(:,2) = [0.8213672_r64 , 0.57735027_r64,  0.14764992_r64]
    redpos(:,3) = [0.8213672_r64 , 0.57735027_r64, -0.14764992_r64]

    ! Init common objects
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i64)

    ! Read the dielectric data from previous G0W0 run
    call load_from_file('head.dat', head)
    call load_from_file('wingL.dat', wingL)
    call load_from_file('wingU.dat', wingU)
    call load_from_file('Binv.dat', Binv)

    ! Do the average
    write(*,*) '[TEST: idiel_t (2d)]'
    do iom = 1, size(head,3)
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom), Binv(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_scrcoulomb_2d(.true.)

        ! Check the head
        rdiff = abs(inv_diel%idiel_head - head_ref(iom))/head_ref(iom)
        write(*,'(A,I2,A,e20.13)')  '  * Regression (HEAD,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST: idiel_t (2d) (HEAD,',iom,'): PASSED]'
        else
            write(*,*)  '[TEST: idiel_t (2d) (HEAD,',iom,'): FAILED]'
            stop 1
        end if

        ! Check that the wings are zzero (though they are by construction)
        rdiff = sum(abs(inv_diel%idiel_wingL))
        write(*,'(A,I3,A, e20.13)')  '  * Regression (WING L,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST: idiel_t (2d) (WING L,',iom,'): PASSED]'
        else
            write(*,*)  '[TEST: idiel_t (2d) (WING L,',iom,'): FAILED]'
            stop 1
        end if

        rdiff = sum(abs(inv_diel%idiel_wingU))
        write(*,'(A,I3,A, e20.13)')  '  * Regression (WING U,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST: idiel_t (2d) (WING U,',iom,'): PASSED]'
        else
            write(*,*)  '[TEST: idiel_t (2d) (WING U,',iom,'): FAILED]'
            stop 1
        end if

    end do

    stop 0

contains 

    !> Process a file and load to anray
    !> @param[in] fname - the file name   
    !> @param[in] data_shape - the shape of the data to reads
    !> @param[out] data - the data
    subroutine load_from_file(fname, data)

        character(len=*), intent(in) :: fname
        complex(r64), allocatable, intent(out) :: data(..)

        integer(i64) :: data_shape(3)
        real(r64) :: dr1, di1
        real(r64), allocatable :: dr2(:), di2(:), dr3(:,:), di3(:,:)
        integer :: fin, iom, niom
        integer :: ii, jj

        open(file=fname, newunit=fin)
        read(fin,*) data_shape

        select rank(data)
        rank(1)
            niom = data_shape(1)
            allocate(data(data_shape(1)))
        rank(2)
            niom = data_shape(2)
            allocate(data(data_shape(1),data_shape(2)))
        rank(3)
            niom = data_shape(3)
            allocate(data(data_shape(1),data_shape(2),data_shape(3)))
        end select

        do iom = 1, niom
            select rank(data)
            rank(1)
                read(fin,*) dr1
                read(fin,*) di1
                data(iom) = cmplx(dr1,di1) 
            rank(2)
                allocate(dr2(data_shape(1)), di2(data_shape(1)))
                read(fin,*) dr2
                read(fin,*) di2
                data(:,iom) = cmplx(dr2,di2)
                deallocate(dr2, di2)
            rank(3)
                allocate(dr3(data_shape(1),data_shape(2)), di3(data_shape(1),data_shape(2)))
                read(fin,*) dr3
                read(fin,*) di3
                data(:,:,iom) = cmplx(dr3,di3)
                deallocate(dr3, di3)
            end select
        end do

        close(fin)

    end subroutine load_from_file

end program test_2d
