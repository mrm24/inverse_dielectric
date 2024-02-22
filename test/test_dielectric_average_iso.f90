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
!> Contains procedures to test averages of the dielectric matrix

!> This module contains the procedures for the test of averaging of the dielectric matrix averages for isotropic case
program test_dielectric_average_iso

    use idiel_constants, only: aip, i32, pi, twopi, zzero
    use idiel, only: idiel_t

    implicit none

    ! GW data
    integer, parameter :: mbsize = 438 ! Mixed basis size
    integer, parameter :: nomega = 12  ! Number of frequencies
    integer(i32) :: iom

    ! Mesh data
    integer(i32), parameter :: ngrid(3) = [6,6,6]

    ! Silicon crystal data
    integer(i32), parameter :: natoms = 2_i32
    real(aip), parameter    :: latpar = 0.513_aip 
    real(aip) :: a(3), b(3), c(3), lattice(3,3)
    real(aip), allocatable :: redpos(:,:)
    integer(i32), allocatable :: types(:)
    ! In case that no SPG we set only Identity (i.e. no symmetry)
    integer(i32) :: nsym = 1
    real(aip)    :: crot(3,3,1)

    ! Dielectric data
    complex(aip), allocatable :: head(:,:,:), wingL(:,:,:)
    complex(aip), allocatable :: wingU(:,:,:), Binv(:,:,:)

    ! Reference data to compare
    real(aip) :: head_ref(12) = [ 6.442031104885004E-002, 0.132222402626860, 0.316390106319300, 0.516127938510703, &
                                  0.645769722906067, 0.707300869613546, 0.733254256631726, 0.785227466265463, &
                                  0.864527818748701, 0.941940576487719, 0.987846211386285, 1.00132129472395 ]
    complex(aip) ::  body_ref(12) =    [(0.99929660640900964_aip,  1.74324304756898034E-021_aip), &
                                        (0.99931353362452746_aip, -4.57572808687006110E-022_aip), &
                                        (0.99934423000483308_aip, -3.15453273667631766E-022_aip), &
                                        (0.99937212492878569_aip, -6.82160838609231033E-022_aip), &
                                        (0.99939388057919998_aip, -4.45690387272341898E-022_aip), &
                                        (0.99940735120538882_aip,  1.44789913499821991E-022_aip), &
                                        (0.99941408651989239_aip,  3.88464023785692835E-022_aip), &
                                        (0.99943077579868456_aip,  4.45218643325182083E-023_aip), &
                                        (0.99947190296600508_aip,  2.63938186362001598E-022_aip), &
                                        (0.99957442290201637_aip, -1.71839395069885518E-022_aip), &
                                        (0.99980527163746624_aip,  1.04080503186989247E-022_aip), &
                                        (0.99998843669895232_aip, -7.63825135003370565E-024_aip) ]
    real(aip) :: rdiff
#ifdef USE_SINGLE_PRECISION
    real(aip), parameter    :: tolerance = 1.0e-4_aip
#else
    real(aip), parameter    :: tolerance = 1.0e-7_aip
#endif
    ! Computation object
    type(idiel_t) :: inv_diel

    ! Lattice vectors
    a(:) = latpar * [0.5_aip,  0.5_aip, 0.0_aip]
    b(:) = latpar * [0.5_aip,  0.0_aip, 0.5_aip]
    c(:) = latpar * [0.0_aip,  0.5_aip, 0.5_aip]


    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [14_i32, 14_i32]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.00_aip, 0.00_aip, 0.00_aip]
    redpos(:,2) = [0.25_aip, 0.25_aip, 0.25_aip]

    ! Init common objects
#ifdef USE_SPGLIB
    call inv_diel%init_common(lattice, redpos, types, ngrid)
#else 
    crot = 0.0_aip
    crot(1,1,1) = 1.0_aip
    crot(2,2,1) = 1.0_aip
    crot(3,3,1) = 1.0_aip
    call inv_diel%init_common(lattice, redpos, types, ngrid, 3_i32, nsym, crot)
#endif
    ! Read the dielectric data from previous G0W0 run
    call load_from_file('head.dat', head)
    call load_from_file('wingL.dat', wingL)
    call load_from_file('wingU.dat', wingU)
    call load_from_file('Binv.dat', Binv)

    ! Do the average
    write(*,*) '[TEST : idiel_t]' 
    do iom = 1, size(head,3)
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom), Binv(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_inversedielectric_3d(.true.)
        
        ! Check the head
        rdiff = abs(inv_diel%idiel_head - head_ref(iom))/head_ref(iom)
        write(*,'(A,I2,A,e20.13)')  '  * Regression (HEAD,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST : idiel_t (HEAD,',iom,'): PASSED]'
        else
            write(*,*)  '[TEST : idiel_t (HEAD,',iom,'): FAILED]'
            stop 1
        end if

        ! Check that the wings are zzero (though they are by construction)
        rdiff = sum(abs(inv_diel%idiel_wingL))
        write(*,'(A,I3,A, e20.13)')  '  * Regression (WING L,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST : idiel_t (WING L,',iom,'): PASSED]'
        else
            write(*,*)  '[TEST : idiel_t (WING L,',iom,'): FAILED]'
            stop 1
        end if

        rdiff = sum(abs(inv_diel%idiel_wingU))
        write(*,'(A,I3,A, e20.13)')  '  * Regression (WING U,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST : idiel_t (WING U,',iom,'): PASSED]'
        else
            write(*,*)  '[TEST : idiel_t (WING U,',iom,'): FAILED]'
            stop 1
        end if

        ! Check the body only one element
        rdiff = abs(inv_diel%idiel_body(1,1) - body_ref(iom))/body_ref(iom)
        write(*,'(A,I3,A, e20.13)')  '  * Regression (BODY,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST : idiel_t (BODY,',iom,'): PASSED]'
        else
            write(*,*)  '[TEST : idiel_t (BODY,',iom,'): FAILED]'
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
        complex(aip), allocatable, intent(out) :: data(:,:,:)

        integer(i32) :: data_shape(3)
        real(aip), allocatable :: dr3(:,:), di3(:,:)
        integer :: fin, iom, niom
        integer :: ii, jj

        open(file=fname, newunit=fin,  status='old', action='read')
        read(fin,*) data_shape

        niom = data_shape(3)
        allocate(data(data_shape(1),data_shape(2),data_shape(3)))

        do iom = 1, niom
            allocate(dr3(data_shape(1),data_shape(2)), di3(data_shape(1),data_shape(2)))
            read(fin,*) dr3
            read(fin,*) di3
            data(:,:,iom) = cmplx(dr3,di3)
            deallocate(dr3, di3)
        end do

        close(fin)

    end subroutine load_from_file

end program test_dielectric_average_iso
