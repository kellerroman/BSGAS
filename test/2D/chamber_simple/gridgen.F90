module types
   use mod_gridgen
   TYPE :: t_start_end
      REAL(KIND=8) ::xstart,xend,ystart,yend
!< GITTERPUNKTE -POSITION (I,J,K,COORD)
   END TYPE
end module types

PROGRAM GridGen
use types
IMPLICIT NONE





INTEGER, PARAMETER :: B2S(4) = (/2,3,5,6/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!                                                                                        !!!!!
!!!!!!                              GITTERABMESSUNGEN                                         !!!!!
!!!!!!                                                                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, PARAMETER :: NUM_OF_BLOCKS = 6


!!!!!!!!!!!!!!!! NJ VALUES 
INTEGER, PARAMETER :: N1 = 10
INTEGER, PARAMETER :: N2 = 10
INTEGER, PARAMETER :: N3 = 10
INTEGER, PARAMETER :: N4 = 20
!!!!!!!!!!!!!!!! NI VALUES 
INTEGER, PARAMETER :: N5 = 10
INTEGER, PARAMETER :: N6 = 50
INTEGER, PARAMETER :: N7 = 10
INTEGER, PARAMETER :: N8 = 2
INTEGER, PARAMETER :: N9 = 10

REAL(KIND=8), PARAMETER :: xm1 = -1.0D0
REAL(KIND=8), PARAMETER :: x0  =  0.0D0
!REAL(KIND=8), PARAMETER :: x2  =  5.0D0
!REAL(KIND=8), PARAMETER :: x3  =  5.4D0
!REAL(KIND=8), PARAMETER :: x4  =  5.5D0
REAL(KIND=8), PARAMETER :: x5  =  6.0D0

REAL(KIND=8), PARAMETER :: y0 = 0.0D0
REAL(KIND=8), PARAMETER :: y1 = 1.0D0
REAL(KIND=8), PARAMETER :: y2 = 2.0D0
REAL(KIND=8), PARAMETER :: y3 = 2.5D0
REAL(KIND=8), PARAMETER :: y4 = 5.0D0
REAL(KIND=8), PARAMETER :: y5 = 2.0D0
REAL(KIND=8), PARAMETER :: y6 = 3.0D0

TYPE(t_start_end) :: block_start_end(NUM_OF_BLOCKS)

INTEGER :: B,I,J,N,gj

REAL(KIND=8) :: y_max,dy
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

call add_block( N6+N7+N8+N9, N1)
block_start_end(2) % xstart = x0
block_start_end(2) % xend   = x5
block_start_end(2) % ystart = y0
block_start_end(2) % yend   = y1

call add_block( N6+N7+N8+N9, N2)
block_start_end(3) % xstart = x0
block_start_end(3) % xend   = x5
block_start_end(3) % ystart = y1
block_start_end(3) % yend   = y2

call add_block( N5, N3)
block_start_end(4) % xstart = xm1
block_start_end(4) % xend   = x0
block_start_end(4) % ystart = y2
block_start_end(4) % yend   = y3

call add_block( N6+N7+N8+N9, N3)
block_start_end(5) % xstart = x0
block_start_end(5) % xend   = x5
block_start_end(5) % ystart = y2
block_start_end(5) % yend   = y3

call add_block( N6+N7+N8+N9, N4)
block_start_end(6) % xstart = x0
block_start_end(6) % xend   = x5
block_start_end(6) % ystart = y3
block_start_end(6) % yend   = y4

call allocate_blocks(0)

DO B = 1, NUM_OF_BLOCKS
call make_block_grid(block_start_end(b),blocks(B))
END DO

do i = n6+1, n6+n7+n8+n9+1
   y_max = y4 + (y5-y4) * dble(i-n6-1) / dble(n7)
   if (i < n6+n7+1) then
   else if (i< n6+n7+n8+1) then
   y_max = y5
   else
   write(*,*) i,i-n6-n7-n8-1
   y_max = y5 + (y6-y5) * dble(i-n6-n7-n8-1) / dble(n9)

   end if

   write(*,*) i,y_max
   dy = y_max / (N1+N2+N3+N4)
   gj = 0
   do n = 1,4
      b = B2S(n)
      do j = 1, blocks(b) % nPkts(2)
         gj = gj + 1
         blocks(b) % xyzs(i,j,1,2) = dy * (gj-1)
      end do
      gj = gj -1
   end do
end do


call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
write(666,'(A)') "wall: 1N 3W 4S 4N 6W 6N ! Wall"
write(666,'(A)') "dn = 5E-2    ! Spacing of first Cell"
close(666)
write(*,*) "done"

WRITE(*,*) "======================================================" &
          ,"FINISHED"  &
          ,"======================================================"

END PROGRAM GridGen

subroutine make_block_grid(se,block)
use types
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
