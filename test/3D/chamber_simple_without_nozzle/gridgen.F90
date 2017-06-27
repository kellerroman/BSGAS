program gridgen_o_grid
   use mod_gridgen
implicit none
type :: t_roundblock
   real(kind=8) ::xstart,xend,ystart,yend
   integer :: le_re
end type
INTEGER, PARAMETER :: NUM_OF_BLOCKS = 14
!!!!!!!!!!!!!!!! NJ VALUES 
integer, parameter :: NJ1 = 15
integer, parameter :: NJ2 = 20
integer, parameter :: NJ3 = 15
integer, parameter :: NJ4 = 30
!!!!!!!!!!!!!!!! NI VALUES 
integer, parameter :: NI1 = 12
integer, parameter :: NI2 = 50

integer, parameter :: NK = 16
integer, parameter :: nVar   = 0
real(kind=8), parameter :: PI = 3.1415927/180.D0

real(kind=8), parameter :: XM1 = -1.0D0
real(kind=8), parameter :: X0  =  0.0D0
real(kind=8), parameter :: X1  =  6.0D0

!real(kind=8), parameter :: Y0 = 0.0D0
real(kind=8), parameter :: Y1 = 1.0D0
real(kind=8), parameter :: Y2 = 2.0D0
real(kind=8), parameter :: Y3 = 2.5D0
real(kind=8), parameter :: Y4 = 5.0D0

TYPE(t_roundblock) :: block_start_end(NUM_OF_BLOCKS)
integer :: im, jm, km

integer :: i, j, k, b
real(kind = 8) :: spkty,epkty,spktz,epktz,lx,ly, winkel, length

WRITE(*,*) "GRIDGEN FÜR EIN EINFACHES O_GRID"

call set_dimension(3)

call add_block(NI1, NK / 2,       NK / 2)
call add_block(NI1, NJ1 - NK / 2, NK / 2)
call add_block(NI1, NJ1 - NK / 2, NK / 2)
block_start_end(1:3) % xstart = XM1
block_start_end(1:3) % xend   =  X0

call add_block(NI2, NK / 2,       NK / 2)
call add_block(NI2, NJ1 - NK / 2, NK / 2)
call add_block(NI2, NJ1 - NK / 2, NK / 2)

block_start_end(4:6) % xstart =  X0 
block_start_end(4:6) % xend   =  X1

b = 7
call add_block(NI2, NJ2,       NK / 2)
block_start_end(b) % xstart =  X0 
block_start_end(b) % xend   =  X1
block_start_end(b) % ystart =  Y1 
block_start_end(b) % yend   =  Y2
block_start_end(b) % le_re  = 1

b = b + 1
call add_block(NI2, NJ2,       NK / 2)
block_start_end(b) % xstart =  X0 
block_start_end(b) % xend   =  X1
block_start_end(b) % ystart =  Y1 
block_start_end(b) % yend   =  Y2
block_start_end(b) % le_re  = 2

b = b + 1
call add_block(NI1, NJ3,       NK / 2)
block_start_end(b) % xstart =  XM1 
block_start_end(b) % xend   =  X0
block_start_end(b) % ystart =  Y2 
block_start_end(b) % yend   =  Y3
block_start_end(b) % le_re  = 1

b = b + 1
call add_block(NI1, NJ3,       NK / 2)
block_start_end(b) % xstart =  XM1 
block_start_end(b) % xend   =  X0
block_start_end(b) % ystart =  Y2 
block_start_end(b) % yend   =  Y3
block_start_end(b) % le_re  = 2

b = b + 1
call add_block(NI2, NJ3,       NK / 2)
block_start_end(b) % xstart =  X0  
block_start_end(b) % xend   =  X1
block_start_end(b) % ystart =  Y2 
block_start_end(b) % yend   =  Y3
block_start_end(b) % le_re  = 1

b = b + 1
call add_block(NI2, NJ3,       NK / 2)
block_start_end(b) % xstart =  X0  
block_start_end(b) % xend   =  X1
block_start_end(b) % ystart =  Y2 
block_start_end(b) % yend   =  Y3
block_start_end(b) % le_re  = 2

b = b + 1
call add_block(NI2, NJ4,       NK / 2)
block_start_end(b) % xstart =  X0  
block_start_end(b) % xend   =  X1
block_start_end(b) % ystart =  Y3 
block_start_end(b) % yend   =  Y4
block_start_end(b) % le_re  = 1

b = b + 1
call add_block(NI2, NJ4,       NK / 2)
block_start_end(b) % xstart =  X0  
block_start_end(b) % xend   =  X1
block_start_end(b) % ystart =  Y3 
block_start_end(b) % yend   =  Y4
block_start_end(b) % le_re  = 2
call allocate_blocks(nVar)

do b = 7, NUM_OF_BLOCKS
call make_block_grid(block_start_end(b) , blocks(b) )
end do

do b = 1,4,3
   im = blocks(b) % nPkts(1)
   jm = blocks(b) % nPkts(2)
   km = blocks(b) % nPkts(3)
   length = block_start_end(b) % xend - block_start_end(b) % xstart
   do k = 1, km
      do j = 1, jm
         do i = 1, im
            blocks(b) % xyzs(i,j,k,1) = dble(i-1)*length      /dble(im-1) + &
            block_start_end(b) % xstart
            blocks(b) % xyzs(i,j,k,2) = dble(j-1)*Y1*0.3d0/dble(jm-1)
            blocks(b) % xyzs(i,j,k,3) = dble(k-1)*Y1*0.3d0/dble(km-1)
         end do
      end do
   end do
end do
do b = 2,5,3
   im = blocks(b) % nPkts(1)
   jm = blocks(b) % nPkts(2)
   km = blocks(b) % nPkts(3)
   length = block_start_end(b) % xend - block_start_end(b) % xstart
   do k = 1, km
      do j = 1, jm
         do i = 1,im
            winkel =  (90.0D0 - 45.0D0 * DBLE(k-1) /DBLE(km-1) )* pi
            spkty = Y1*0.3D0
            spktz = DBLE(k-1)*Y1*0.3D0/DBLE(km-1)
            epkty = Y1 * sin(winkel)
            epktz = Y1 * cos(winkel)
            if (k == 1) then
               epkty = Y1
               epktz = 0.0d0
            end if
            lx = epkty-spkty
            ly = epktz-spktz
            blocks(b) % xyzs(i,j,k,1) = dble(i-1)*length      /dble(im-1) + &
            block_start_end(b) % xstart
            blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(j-1)/DBLE(jm-1) * lx 
            blocks(b) % xyzs(i,j,k,3) = spktz + DBLE(j-1)/DBLE(jm-1) * ly
         end do
      end do
   end do
end do

do b = 3,6,3
   im = blocks(b) % nPkts(1)
   jm = blocks(b) % nPkts(2)
   km = blocks(b) % nPkts(3)
   length = block_start_end(b) % xend - block_start_end(b) % xstart
   do k = 1, km
      do j = 1, jm
         do i = 1,im
            winkel =  (45.0d0 - 45.0d0 * dble(k-1) / dble(km-1)) * pi
            spkty = dble(km-k) * Y1 * 0.3d0 / dble(km-1)
            spktz = Y1 * 0.3d0
            epkty = Y1 * sin(winkel)
            epktz = Y1 * cos(winkel)
            if (k == km) then
               epkty = 0.0d0
               epktz = Y1
            end if
            lx = epkty-spkty
            ly = epktz-spktz

            blocks(b) % xyzs(i,j,k,1) = dble(i-1)*length      /dble(im-1) + &
            block_start_end(b) % xstart
            blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(jm-1) * lx 
            blocks(b) % xyzs(i,j,k,3) = spktz + dble(j-1)/dble(jm-1) * ly
         end do
      end do
   end do
end do

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
write(666,'(A)') "wall: 2N, 3N ,7W,8W,9S,9N,10S,10N,13E,13N,14E,14N"
write(666,'(A)') "dn = 2E-2    ! Spacing of first Cell"
close(666)
write(*,*) "fdone"


contains
subroutine make_block_grid(se,block)
implicit none
real(kind=8), parameter :: PI = 3.1415927/180.D0
TYPE(t_roundblock) :: se
TYPE(tblock) :: block


integer :: i,j,k
real(kind=8) :: dx,dy, angle, pos,aa

dx = ( se % xend - se % xstart ) / block % nCells(1)
dy = ( se % yend - se % ystart ) / block % nCells(2)

if (se % le_re == 1) then
   aa = 90.0D0
else
   aa = 45.0D0
end if

do k = 1, block % nPkts(3) 
   angle = (aa - 45.0D0  * dble(k-1)/dble(block % nCells(3))) * PI
   do j = 1, block % nPkts(2)
      pos = se % ystart + dy * (j-1)
      do i = 1, block % nPkts(1)
         block % xyzs(i,j,k,1) = se % xstart + dx * (i-1)
         block % xyzs(i,j,k,2) = pos * sin(angle)
         block % xyzs(i,j,k,3) = pos * cos(angle)
      end do
   end do
end do
end subroutine make_block_grid
end program
