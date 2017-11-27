module types
   TYPE :: t_start_end
      REAL(KIND=8) ::xstart,xend,ystart,yend
!< GITTERPUNKTE -POSITION (I,J,K,COORD)
   END TYPE
end module types

PROGRAM GridGen
use types
use mod_gridgen
IMPLICIT NONE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!                                                                                        !!!!!
!!!!!!                              GITTERABMESSUNGEN                                         !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, PARAMETER :: NUM_OF_BLOCKS = 3


!!!!!!!!!!!!!!!! nj values 
integer                 :: nj
!!!!!!!!!!!!!!!! ni values 
integer                 :: ni1
integer                 :: ni2

real(kind=8)            :: x1
real(kind=8), parameter :: x2  =  0.0d0
real(kind=8)            :: x3

real(kind=8), parameter :: y2 = 0.0d0
real(kind=8)            :: y3

TYPE(t_start_end) :: block_start_end(NUM_OF_BLOCKS)

integer :: b, testcase
character(len = 1) :: arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!                                                                                        !!!!!
!!!!!!                    BERECHUNG DER KLEINSTEN ZELLGRÖßE IN RADIALER RICHTUNG              !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!       Annahme: konstante Vergrößerung für bestimmte Anzahl an Zellen, dann konst.      !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*) "======================================================" &
          ,"==   GRIDGEN for Channel with changing Refinement   =="  &
          ,"======================================================"
    testcase = 1
call get_command_argument(1, arg)
if (len_trim(arg) /= 0) then
   read(arg,*) testcase
end if

if (testcase == 1) then
   nj  =  50
   ni1 =  50
   ni2 =  50

   x1 = -5.000d-02
   x3 =  5.000d-02
   y3 =  5.080d-02
else if (testcase == 2) then
   nj  =  50
   ni1 =  10
   ni2 =  50

   x1 = -1.000d-02
   x3 =  5.000d-02
   y3 =  5.080d-02
else
   write(*,*) "case unknown:",testcase
   stop 1
end if

call set_dimension(2)
call add_block(NI1,NJ)
block_start_end(1) % xstart = x1
block_start_end(1) % xend   = x2
block_start_end(1) % ystart = y2
block_start_end(1) % yend   = y3


call add_block( NI2,NJ)
block_start_end(2) % xstart = x2
block_start_end(2) % xend   = x3
block_start_end(2) % ystart = y2
block_start_end(2) % yend   = y3


call allocate_blocks(0)

do B = 1, NUM_OF_BLOCKS
   call make_block_grid(block_start_end(b),blocks(B))
end do

call write_grid()
open(666,file="bc.cfg")
if (testcase == 1) then
   write(666,'(A)') "fixed-wall: 1S     ! Wall"
   write(666,'(A)') "dn = 1E-6    ! Spacing of first Cell"
   write(666,'(A)') "pos-fixed: 2S,1N,2N     ! Position on the Boundary Fixed"
   !write(666,'(A)') "pos-fixed: 2S     ! Position on the Boundary Fixed"
else if (testcase == 2) then
   write(666,'(A)') "fixed-wall: 1S     ! Wall"
   write(666,'(A)') "dn = 1E-6    ! Spacing of first Cell"
   write(666,'(A)') "x-pos-fixed: 2S,1N,2N     ! Position on the Boundary Fixed"
   !write(666,'(A)') "pos-fixed: 2S     ! Position on the Boundary Fixed"
end if
close(666)

WRITE(*,*) "======================================================" &
          ,"==                    FINISHED                      ==" &
          ,"======================================================"

END PROGRAM GridGen

subroutine make_block_grid(se,block)
use types
use mod_gridgen
implicit none
TYPE(t_start_end) :: se
TYPE(tblock) :: block


integer :: i,j
real(kind=8) :: dx,dy

dx = ( se % xend - se % xstart ) / block % nCells(1)
dy = ( se % yend - se % ystart ) / block % nCells(2)
do j = 1, block%nPkts(2)
   do i = 1, block%nPkts(1)
      block % xyzs(i,j,1,1) = se % xstart + dx * (i-1)
      block % xyzs(i,j,1,2) = se % ystart + dy * (j-1)
      block % xyzs(i,j,1,3) = 0.0d0
   end do
end do
end subroutine make_block_grid
