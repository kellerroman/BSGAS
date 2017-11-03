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
INTEGER, PARAMETER :: NUM_OF_BLOCKS = 6


!!!!!!!!!!!!!!!! NJ VALUES 
INTEGER, PARAMETER :: N1 = 15
INTEGER, PARAMETER :: N2 = 20
INTEGER, PARAMETER :: N3 = 15
INTEGER, PARAMETER :: N4 = 30
!!!!!!!!!!!!!!!! NI VALUES 
INTEGER, PARAMETER :: N5 = 15
INTEGER, PARAMETER :: N6 = 50

REAL(KIND=8), PARAMETER :: xm1 = -1.0D0
REAL(KIND=8), PARAMETER :: x0  =  0.0D0
REAL(KIND=8), PARAMETER :: x5  =  6.0D0

REAL(KIND=8), PARAMETER :: y0 = 0.0D0
REAL(KIND=8), PARAMETER :: y1 = 1.0D0
REAL(KIND=8), PARAMETER :: y2 = 2.0D0
REAL(KIND=8), PARAMETER :: y3 = 2.5D0
REAL(KIND=8), PARAMETER :: y4 = 5.0D0

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
          ,"GRIDGEN"  &
          ,"======================================================"

call set_dimension(2)
call add_block(N5,N1)
block_start_end(1) % xstart = xm1
block_start_end(1) % xend   = x0
block_start_end(1) % ystart = y0
block_start_end(1) % yend   = y1

call add_block( N6, N1)
block_start_end(2) % xstart = x0
block_start_end(2) % xend   = x5
block_start_end(2) % ystart = y0
block_start_end(2) % yend   = y1

call add_block( N6, N2)
block_start_end(3) % xstart = x0
block_start_end(3) % xend   = x5
block_start_end(3) % ystart = y1
block_start_end(3) % yend   = y2

call add_block( N5, N3)
block_start_end(4) % xstart = xm1
block_start_end(4) % xend   = x0
block_start_end(4) % ystart = y2
block_start_end(4) % yend   = y3

call add_block( N6, N3)
block_start_end(5) % xstart = x0
block_start_end(5) % xend   = x5
block_start_end(5) % ystart = y2
block_start_end(5) % yend   = y3

call add_block( N6, N4)
block_start_end(6) % xstart = x0
block_start_end(6) % xend   = x5
block_start_end(6) % ystart = y3
block_start_end(6) % yend   = y4

call allocate_blocks(0)

DO B = 1, NUM_OF_BLOCKS
call make_block_grid(block_start_end(b),blocks(B))
END DO

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
write(666,'(A)') "wall: 1N,3W,4S,4N,6W,6N ! Wall"
write(666,'(A)') "dn = 2E-2    ! Spacing of first Cell"
close(666)
write(*,*) "done"

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
