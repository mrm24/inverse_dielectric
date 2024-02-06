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
!> Contains program to see procedure 2d for a real case (phosphorene)

!> Program showcasing th 2D routines
program test_2d

    use idiel_constants, only: r64, i64, pi, twopi, zzero
    use idiel, only: idiel_t

    implicit none

    ! GW data
    integer(i64), parameter :: nomega = 32_i64  ! Number of frequencies
    integer(i64) :: iom

    ! Mesh data
    integer(i64), parameter :: ngrid(3) = [18_i64, 14_i64, 1_i64]

    ! Silicon crystal data
    integer(i64), parameter :: natoms = 4_i64
    real(r64) :: a(3), b(3), c(3), lattice(3,3)
    real(r64), allocatable :: redpos(:,:)
    integer(i64), allocatable :: types(:)

    ! Dielectric data
    complex(r64), allocatable :: head(:,:,:), wingL(:,:,:)
    complex(r64), allocatable :: wingU(:,:,:), body(:,:,:)

    ! Computation object
    type(idiel_t) :: inv_diel

    ! Lattice vectors
    a(:) = [6.26122957273549_r64, 0.00000000000000_r64, 0.000000000000000_r64]
    b(:) = [0.00000000000000_r64, 8.26566207441072_r64, 0.000000000000000_r64]
    c(:) = [0.00000000000000_r64, 0.00000000000000_r64, 19.79110171589015_r64]

    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [15_i64, 15_i64, 15_i64, 15_i64]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.50000000000000_r64, 0.08060000000000_r64, 0.60340000000000_r64]
    redpos(:,2) = [0.00000000000000_r64, 0.58060000000000_r64, 0.39660000000000_r64]
    redpos(:,3) = [0.00000000000000_r64, 0.41940000000000_r64, 0.60340000000000_r64]
    redpos(:,4) = [0.50000000000000_r64, 0.91940000000000_r64, 0.39660000000000_r64]

    ! Reading files
    call load_from_file('head',  head)
    call load_from_file('wingL', wingL)
    call load_from_file('wingU', wingU)
    call load_from_file('body',  body)

    !!!!!!!!!!!! Here it goes the main programm !!!!!!!!!!!!!!!!!!!

    ! Init common objects
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i64)

    ! Do the average (Hermitian)
    do iom = 1, nomega
        ! Invert body
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_scrcoulomb_2d(.true.)
    end do

    ! Do the average (General case)
    do iom = 1, nomega
        ! Invert body
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_scrcoulomb_2d(.false.)
    end do

contains 

    !> Process a file and load to anray
    !> @param[in] fname - the file name   
    !> @param[in] data_shape - the shape of the data to reads
    !> @param[out] data - the data
    subroutine load_from_file(fname, data)

        character(len=*), intent(in) :: fname
        complex(r64), allocatable, intent(out) :: data(:,:,:)

        integer :: fin
        real(r64), allocatable  :: dreal(:,:,:), dimag(:,:,:)
        integer(i64) :: data_shape(3)
        
        open(file=fname, newunit=fin, status='old', action='read')

        read(fin, *) data_shape

        allocate(data(data_shape(1),data_shape(2),data_shape(3)))
        allocate(dreal(data_shape(1),data_shape(2),data_shape(3)))
        allocate(dimag, mold=dreal)
        read(fin, *) dreal
        read(fin, *) dimag
        data = cmplx(dreal, dimag, r64) 
        
        close(fin)

    end subroutine load_from_file

end program test_2d
