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
integer                 :: nj1
integer                 :: nj2
!!!!!!!!!!!!!!!! ni values 
integer                 :: ni1
integer                 :: ni2

real(kind=8)            :: x1
real(kind=8), parameter :: x2  =  0.0d0
real(kind=8)            :: x3

real(kind=8)            :: y1
real(kind=8), parameter :: y2 = 0.0d0
real(kind=8)            :: y3

TYPE(t_start_end) :: block_start_end(NUM_OF_BLOCKS)

integer :: b, bfs_case
character(len = 1) :: arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!                                                                                        !!!!!
!!!!!!                    BERECHUNG DER KLEINSTEN ZELLGRÖßE IN RADIALER RICHTUNG              !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!       Annahme: konstante Vergrößerung für bestimmte Anzahl an Zellen, dann konst.      !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*) "======================================================" &
          ,"==       GRIDGEN for Backward Facing Step           =="  &
          ,"======================================================"
bfs_case = 1
call get_command_argument(1, arg)
if (len_trim(arg) /= 0) then
   read(arg,*) bfs_case
end if

if (bfs_case == 1) then
   nj1 = 100
   nj2 = 200
   ni1 = 100
   ni2 = 500

   x1 = -5.080d-02
   x3 =  3.810d-01
   y1 = -1.270d-02
   y3 =  1.016d-01
else
   nj1 =  50
   nj2 = 100
   ni1 =  50
   ni2 = 100

   x1 = -2.540d-02
   x3 =  7.620d-02
   y1 = -6.350d-03
   y3 =  5.080d-02
end if

call set_dimension(2)
call add_block(NI1,NJ2)
block_start_end(1) % xstart = x1
block_start_end(1) % xend   = x2
block_start_end(1) % ystart = y2
block_start_end(1) % yend   = y3

call add_block( NI2,NJ1)
block_start_end(2) % xstart = x2
block_start_end(2) % xend   = x3
block_start_end(2) % ystart = y1
block_start_end(2) % yend   = y2

call add_block( NI2,NJ2)
block_start_end(3) % xstart = x2
block_start_end(3) % xend   = x3
block_start_end(3) % ystart = y2
block_start_end(3) % yend   = y3


call allocate_blocks(0)

do B = 1, NUM_OF_BLOCKS
   call make_block_grid(block_start_end(b),blocks(B))
end do

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "!wall: 1N,1S,2S,2W,3N ! Wall"
write(666,'(A)') "!dn = 1E-6    ! Spacing of first Cell"
close(666)

WRITE(*,*) "======================================================" &
          ,"FINISHED"  &
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
