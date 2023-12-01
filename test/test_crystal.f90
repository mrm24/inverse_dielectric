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
!> Tests the generated crystal structure as well as the symmetry and compares it with reference results
program test_crystal
    
    use idiel_constants,         only: i64, r64, pi, twopi
    use idiel_crystal_cell,      only: cell_t
    use idiel_crystal_symmetry,  only: symmetry_t

    implicit none

    real(r64), parameter    :: tolerance = 1.0e-12_r64
    real(r64) :: rdiff

    ! Elements
    type(cell_t)     :: my_cell
    type(symmetry_t) :: my_symmetry

    ! cell_t elements
    integer(i64), parameter :: natoms = 4_i64
    real(r64) :: a(3) , b(3), c(3)
    real(r64) :: lattice(3,3) 
    real(r64), allocatable :: redpos(:,:)
    integer(r64), allocatable :: types(:)

    ! ref_values
    real(r64), parameter :: vuc_ref = 85.990737984502_r64
    real(r64) :: invlat(3,3)

    ! symmetry_t reference values
    integer(i64), parameter :: space_group = 186_i64
    integer(i64), parameter :: nsym = 12_i64
    integer(i64) :: rot_ref(3,3,12)
    real(r64)    :: crot_ref(3,3,12)

    ! Locals
    integer(i64) :: ii

    invlat(1,:) = [1.601845996473_r64, 0.924826217264_r64, 0.000000000000_r64]
    invlat(2,:) = [0.000000000000_r64, 1.849652434528_r64, 0.000000000000_r64]
    invlat(3,:) = [0.000000000000_r64, 0.000000000000_r64, 0.973592098042_r64]

    rot_ref(1,:,1) = [1_i64, 0_i64, 0_i64]
    rot_ref(2,:,1) = [0_i64, 1_i64, 0_i64]
    rot_ref(3,:,1) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,2) = [1_i64, -1_i64, 0_i64]
    rot_ref(2,:,2) = [1_i64, 0_i64, 0_i64]
    rot_ref(3,:,2) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,3) = [0_i64, -1_i64, 0_i64]
    rot_ref(2,:,3) = [1_i64, -1_i64, 0_i64]
    rot_ref(3,:,3) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,4) = [-1_i64, 0_i64, 0_i64]
    rot_ref(2,:,4) = [0_i64, -1_i64, 0_i64]
    rot_ref(3,:,4) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,5) = [-1_i64, 1_i64, 0_i64]
    rot_ref(2,:,5) = [-1_i64, 0_i64, 0_i64]
    rot_ref(3,:,5) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,6) = [0_i64, 1_i64, 0_i64]
    rot_ref(2,:,6) = [-1_i64, 1_i64, 0_i64]
    rot_ref(3,:,6) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,7) = [0_i64, 1_i64, 0_i64]
    rot_ref(2,:,7) = [1_i64, 0_i64, 0_i64]
    rot_ref(3,:,7) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,8) = [1_i64, 0_i64, 0_i64]
    rot_ref(2,:,8) = [1_i64, -1_i64, 0_i64]
    rot_ref(3,:,8) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,9) = [1_i64, -1_i64, 0_i64]
    rot_ref(2,:,9) = [0_i64, -1_i64, 0_i64]
    rot_ref(3,:,9) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,10) = [0_i64, -1_i64, 0_i64]
    rot_ref(2,:,10) = [-1_i64, 0_i64, 0_i64]
    rot_ref(3,:,10) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,11) = [-1_i64, 0_i64, 0_i64]
    rot_ref(2,:,11) = [-1_i64, 1_i64, 0_i64]
    rot_ref(3,:,11) = [0_i64, 0_i64, 1_i64]
    rot_ref(1,:,12) = [-1_i64, 1_i64, 0_i64]
    rot_ref(2,:,12) = [0_i64, 1_i64, 0_i64]
    rot_ref(3,:,12) = [0_i64, 0_i64, 1_i64]

    crot_ref(1,:,1) = [0.9999999999999999_r64, -1.2798967289084123e-16_r64, 0.0_r64]
    crot_ref(2,:,1) = [0.0_r64, 0.9999999999999999_r64, 0.0_r64]
    crot_ref(3,:,1) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,2) = [0.49999999999999994_r64, -0.8660254037844387_r64, 0.0_r64]
    crot_ref(2,:,2) = [0.8660254037844385_r64, 0.4999999999999999_r64, 0.0_r64]
    crot_ref(3,:,2) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,3) = [-0.49999999999999994_r64, -0.8660254037844386_r64, 0.0_r64]
    crot_ref(2,:,3) = [0.8660254037844385_r64, -0.5_r64, 0.0_r64]
    crot_ref(3,:,3) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,4) = [-0.9999999999999999_r64, 1.2798967289084123e-16_r64, 0.0_r64]
    crot_ref(2,:,4) = [0.0_r64, -0.9999999999999999_r64, 0.0_r64]
    crot_ref(3,:,4) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,5) = [-0.49999999999999994_r64, 0.8660254037844387_r64, 0.0_r64]
    crot_ref(2,:,5) = [-0.8660254037844385_r64, -0.4999999999999999_r64, 0.0_r64]
    crot_ref(3,:,5) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,6) = [0.49999999999999994_r64, 0.8660254037844386_r64, 0.0_r64]
    crot_ref(2,:,6) = [-0.8660254037844385_r64, 0.5_r64, 0.0_r64]
    crot_ref(3,:,6) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,7) = [-0.49999999999999994_r64, 0.8660254037844387_r64, 0.0_r64]
    crot_ref(2,:,7) = [0.8660254037844385_r64, 0.4999999999999999_r64, 0.0_r64]
    crot_ref(3,:,7) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,8) = [0.49999999999999994_r64, 0.8660254037844386_r64, 0.0_r64]
    crot_ref(2,:,8) = [0.8660254037844385_r64, -0.5_r64, 0.0_r64]
    crot_ref(3,:,8) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,9) = [0.9999999999999999_r64, -1.2798967289084123e-16_r64, 0.0_r64]
    crot_ref(2,:,9) = [0.0_r64, -0.9999999999999999_r64, 0.0_r64]
    crot_ref(3,:,9) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,10) = [0.49999999999999994_r64, -0.8660254037844387_r64, 0.0_r64]
    crot_ref(2,:,10) = [-0.8660254037844385_r64, -0.4999999999999999_r64, 0.0_r64]
    crot_ref(3,:,10) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,11) = [-0.49999999999999994_r64, -0.8660254037844386_r64, 0.0_r64]
    crot_ref(2,:,11) = [-0.8660254037844385_r64, 0.5_r64, 0.0_r64]
    crot_ref(3,:,11) = [0.0_r64, 0.0_r64, 1.0_r64]
    crot_ref(1,:,12) = [-0.9999999999999999_r64, 1.2798967289084123e-16_r64, 0.0_r64]
    crot_ref(2,:,12) = [0.0_r64, 0.9999999999999999_r64, 0.0_r64]
    crot_ref(3,:,12) = [0.0_r64, 0.0_r64, 1.0_r64]


    write(*,*) '[TEST : cell_t and symmetry_t]' 
    write(*,*) ' * kind : regression test against precomputed data'
    write(*,'(A, E11.6)') '  * tolerance : ', tolerance

    ! Data for WZ ZnSe

    ! Lattice vectors
    a(:) = 3.9224652813161520_r64 * [1.0_r64,  0.0_r64,               0.0_r64]
    b(:) = 3.9224652813161520_r64 * [-0.5_r64, sqrt(3.0_r64)/2.0_r64, 0.0_r64]
    c(:) = 3.9224652813161520_r64 * [0.0_r64,  0.0_r64,               1.6452947797075332_r64]

    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [30_i64, 30_i64, 34_i64, 34_i64]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.3333333333333357_r64, 0.6666666666666643_r64, 0.0000000000000000_r64]
    redpos(:,2) = [0.6666666666666643_r64, 0.3333333333333357_r64, 0.5000000000000000_r64]
    redpos(:,3) = [0.3333333333333357_r64, 0.6666666666666643_r64, 0.3736232276372607_r64]
    redpos(:,4) = [0.6666666666666643_r64, 0.3333333333333357_r64, 0.8736232276372607_r64]

    ! Initialize the crystal structure
    call my_cell%initialize(lattice, redpos, types)

    ! Check computed values
    rdiff = abs(my_cell%vuc - vuc_ref)/vuc_ref
    write(*,'(A, E11.6)')  '  * Regression (vuc) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : cell_t%vuc : PASSED]'
    else
        write(*,*)  '[TEST : cell_t%vuc : FAILED]'
        stop 1
    end if

    rdiff = relative_diference_m(my_cell%rlattice,invlat)
    write(*,'(A, E11.6)')  '  * Regression (rlattice) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : cell_t%rlattice : PASSED]'
    else
        write(*,*)  '[TEST : cell_t%rlattice : FAILED]'
        stop 1
    end if

    ! Computing symmetry stuff
    call my_symmetry%initialize(my_cell)

    ! Check the space group
    if ( my_symmetry%spgid .eq. space_group) then
        write(*,*)  '[TEST : symmetry_t%spgid : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%spgid : FAILED]'
        stop 1
    end if

    !Check the rotations
    if ( my_symmetry%nsym .eq. nsym) then
        write(*,*)  '[TEST : symmetry_t%nsym : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%nsym : FAILED]'
        stop 1
    end if

    rdiff = relative_diference_i3(my_symmetry%rot, rot_ref)
    write(*,'(A, E11.6)')  '  * Regression (rotations) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : symmetry_t%rot : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%rot : FAILED]'
        stop 1
    end if

    rdiff = relative_diference_r3(my_symmetry%crot, crot_ref)
    write(*,'(A, E11.6)')  '  * Regression (Cartesian rotations) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : symmetry_t%crot : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%crot : FAILED]'
        stop 1
    end if

    stop 0
contains
        !> Computes the Frobenious norm of a matrix
        pure function relative_diference_m(array,reference) result(rd)
            !> The array for which the Frobenious norm will be computed
            real(r64), intent(in) :: array(:,:)
            real(r64), intent(in) :: reference(:,:)
            real(r64) :: rd
            rd = sqrt(real(sum( (array-reference)**2 ))) / sqrt(real(sum((reference)**2 )))
        end function relative_diference_m
        
        !> Computes the Frobenious norm of a matrix
        pure function relative_diference_i3(array,reference) result(rd)
            !> The array for which the Frobenious norm will be computed
            integer(i64), intent(in) :: array(:,:,:)
            integer(i64), intent(in) :: reference(:,:,:)
            real(r64) :: rd
            rd = sqrt(real(sum( (array-reference)**2 ))) / sqrt(real(sum((reference)**2 )))
        end function relative_diference_i3

        !> Computes the Frobenious norm of a matrix
        pure function relative_diference_r3(array,reference) result(rd)
            !> The array for which the Frobenious norm will be computed
            real(r64), intent(in) :: array(:,:,:)
            real(r64), intent(in) :: reference(:,:,:)
            real(r64) :: rd
            rd = sqrt(real(sum( (array-reference)**2 ))) / sqrt(real(sum((reference)**2 )))
        end function relative_diference_r3

end program test_crystal