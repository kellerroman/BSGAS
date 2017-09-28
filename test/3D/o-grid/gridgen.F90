program gridgen_o_grid
   use mod_gridgen
implicit none
integer, parameter :: NI = 8
integer, parameter :: NJ = 32
integer, parameter :: NK = 16
!integer, parameter :: NK = 8
!integer, parameter :: NY = 7
integer, parameter :: NB = 3
integer, parameter :: nVar   = 0
real(kind = 8), parameter :: Pi = 3.1415927/180.D0
real(kind = 8), parameter :: RADIUS = 2.0000000D-0
real(kind = 8), parameter :: LENGTH = 1.0000000D-0
integer :: im, jm, km

integer :: i, j, k, b
real(kind = 8) :: spkty,epkty,spktz,epktz,lx,ly, winkel
WRITE(*,*) "GRIDGEN FÜR EIN EINFACHES O_GRID"

call set_dimension(3)

call add_block(NI, NK / 2,      NK / 2)
call add_block(NI, NJ - NK / 2, NK / 2)
call add_block(NI, NJ - NK / 2, NK / 2)

call allocate_blocks(nVar)

write(*,'(2(I4,1X))') ((blocks(i) % nCells(j),J=1,2),I=1,NB)

b = 1
im = blocks(b) % nPkts(1)
jm = blocks(b) % nPkts(2)
km = blocks(b) % nPkts(3)
do k = 1, km
   do j = 1, jm
      do i = 1, im
         blocks(b) % xyzs(i,j,k,1) = dble(i-1)*LENGTH      /dble(im-1)
         blocks(b) % xyzs(i,j,k,2) = dble(j-1)*RADIUS*0.3d0/dble(jm-1)
         blocks(b) % xyzs(i,j,k,3) = dble(k-1)*RADIUS*0.3d0/dble(km-1)
      end do
   end do
end do
b = 2
im = blocks(b) % nPkts(1)
jm = blocks(b) % nPkts(2)
km = blocks(b) % nPkts(3)
do k = 1, km
   do j = 1, jm
      do i = 1,im
         winkel =  (90.0D0 - 45.0D0 * DBLE(k-1) /DBLE(km-1) )* pi
         spkty = RADIUS*0.3D0
         spktz = DBLE(k-1)*RADIUS*0.3D0/DBLE(km-1)
         epkty = RADIUS * sin(winkel)
         epktz = RADIUS * cos(winkel)
         if (k == 1) then
            epkty = RADIUS
            epktz = 0.0d0
         end if
         lx = epkty-spkty
         ly = epktz-spktz
         blocks(b) % xyzs(i,j,k,1) = dble(i-1)*LENGTH      /dble(im-1)
         blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(j-1)/DBLE(jm-1) * lx 
         blocks(b) % xyzs(i,j,k,3) = spktz + DBLE(j-1)/DBLE(jm-1) * ly
      end do
   end do
end do

b = 3
im = blocks(b) % nPkts(1)
jm = blocks(b) % nPkts(2)
km = blocks(b) % nPkts(3)
do k = 1, km
   do j = 1, jm
      do i = 1,im
         winkel =  (45.0d0 - 45.0d0 * dble(k-1) / dble(km-1)) * pi
         spkty = dble(km-k) * radius * 0.3d0 / dble(km-1)
         spktz = RADIUS * 0.3d0
         epkty = RADIUS * sin(winkel)
         epktz = RADIUS * cos(winkel)
         if (k == km) then
            epkty = 0.0d0
            epktz = radius
         end if
         lx = epkty-spkty
         ly = epktz-spktz

         blocks(b) % xyzs(i,j,k,1) = dble(i-1)*LENGTH      /dble(im-1)
         blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(jm-1) * lx 
         blocks(b) % xyzs(i,j,k,3) = spktz + dble(j-1)/dble(jm-1) * ly
      end do
   end do
end do

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
write(666,'(A)') "wall: 2N, 3N ! Wall"
write(666,'(A)') "dn = 1E-2    ! Spacing of first Cell"
write(666,'(A)') "inflow: 1E, 2E, 3E ! INFLOW"
write(666,'(A)') "dn = 5E-2    ! Spacing of first Cell"
close(666)
write(*,*) "done"


end program
