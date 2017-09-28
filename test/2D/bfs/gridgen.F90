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


!!!!!!!!!!!!!!!! NJ VALUES 
INTEGER, PARAMETER :: NJ1 = 100
INTEGER, PARAMETER :: NJ2 = 200
!!!!!!!!!!!!!!!! NI VALUES 
INTEGER, PARAMETER :: NI1 = 100
INTEGER, PARAMETER :: NI2 = 500

REAL(KIND=8), PARAMETER :: x1 = -5.08D-2
REAL(KIND=8), PARAMETER :: x2  =  0.0D0
REAL(KIND=8), PARAMETER :: x3  =  3.81D-1

REAL(KIND=8), PARAMETER :: y1 = -1.27D-2
REAL(KIND=8), PARAMETER :: y2 = 0.0D0
REAL(KIND=8), PARAMETER :: y3 = 1.016D-1

TYPE(t_start_end) :: block_start_end(NUM_OF_BLOCKS)

INTEGER :: B

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
write(666,'(A)') "wall: 1N,1S,2S,2W,3N ! Wall"
write(666,'(A)') "dn = 1E-6    ! Spacing of first Cell"
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
