module types
   TYPE :: t_start_end
      REAL(KIND=8) ::xstart,xend,ystart,yend
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
INTEGER, PARAMETER :: N1 = 10
INTEGER, PARAMETER :: N2 = 20
!!!!!!!!!!!!!!!! NI VALUES 
INTEGER, PARAMETER :: N3 = 10
INTEGER, PARAMETER :: N4 = 10

REAL(KIND=8), PARAMETER :: xm1 = -1.0D0
REAL(KIND=8), PARAMETER :: x0  =  0.0D0
REAL(KIND=8), PARAMETER :: x1  =  5.0D-1

REAL(KIND=8), PARAMETER :: y0 = 2.0D0
REAL(KIND=8), PARAMETER :: y1 = 2.5D0
REAL(KIND=8), PARAMETER :: y2 = 3.5D0

TYPE(t_start_end) :: block_start_end(NUM_OF_BLOCKS)

INTEGER :: B

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!                                                                                        !!!!!
!!!!!!                    BERECHUNG DER KLEINSTEN ZELLGRÖßE IN RADIALER RICHTUNG              !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!       Annahme: konstante Vergrößerung für bestimmte Anzahl an Zellen, dann konst.      !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*) "================================================" &
          ,"GRIDGEN FOR SIMPLE INJECTOR"  &
          ,"================================================"

call set_dimension(2)
call add_block(N3,N1)
block_start_end(1) % xstart = xm1
block_start_end(1) % xend   = x0
block_start_end(1) % ystart = y0
block_start_end(1) % yend   = y1

call add_block(N4, N1)
block_start_end(2) % xstart = x0
block_start_end(2) % xend   = x1
block_start_end(2) % ystart = y0
block_start_end(2) % yend   = y1

call add_block(N4, N2)
block_start_end(3) % xstart = x0
block_start_end(3) % xend   = x1
block_start_end(3) % ystart = y1
block_start_end(3) % yend   = y2

call allocate_blocks(0)

DO B = 1, NUM_OF_BLOCKS
call make_block_grid(block_start_end(b) ,blocks(b))
END DO

call write_grid()

open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
write(666,'(A)') "wall: 1N 3W ! Wall"
write(666,'(A)') "dn = 5E-2    ! Spacing of first Cell"
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
