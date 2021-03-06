program gridgen_o_grid
   use mod_gridgen
implicit none
integer, parameter :: NB = 3
integer, parameter :: nVar   = 0
real(kind = 8), parameter :: Pi = 3.1415927/180.D0

integer:: nk = 16
integer:: ny = 32
real(kind = 8) :: radius
integer :: ni, nj

integer :: i, j, k, b
real(kind = 8) :: spktx,spkty,epkty,epktx,lx,ly, winkel
character(len = 1) :: arg
integer :: testcase
WRITE(*,*) "GRIDGEN FÜR EIN EINFACHES O_GRID"

testcase = 1
call get_command_argument(1, arg)
if (len_trim(arg) /= 0) then
   read(arg,*) testcase
end if
if (testcase == 1) then
   nk = 16
   ny = 24
   radius = 2.0000000d-0
else
   nk = 32
   ny = 240
   radius = 19.0500000d-3
end if
call set_dimension(2)

call add_block(nk / 2, nk / 2)
call add_block(nk / 2, ny)
call add_block(nk / 2, ny)

call allocate_blocks(nVar)

write(*,'(2(I4,1X))') ((blocks(i) % nCells(j),J=1,2),I=1,NB)

b = 1
ni = blocks(b) % nPkts(1)
nj = blocks(b) % nPkts(2)
 k = 1
do j = 1, nj
   do i = 1,ni
      blocks(b) % xyzs(i,j,k,1) = dble(i-1)*RADIUS*0.3d0/dble(ni-1)
      blocks(b) % xyzs(i,j,k,2) = dble(j-1)*RADIUS*0.3d0/dble(nj-1)
   end do
end do

b = 2
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (45.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      spkty = dble(ni-i) * radius * 0.3d0 / dble(ni-1)
      spktx = radius * 0.3d0
      epktx = radius * sin(winkel)
      epkty = radius * cos(winkel)
      if (i == ni) then
         epktx = radius
         epkty = 0.0d0
      end if
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx 
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly
   end do
end do

b = 3
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  45.0D0 * DBLE(I-1) /DBLE(NI-1) * pi
      spktx = DBLE(I-1)*RADIUS*0.3D0/DBLE(NI-1)
      spkty = RADIUS*0.3D0
      epktx = RADIUS * SIN(WINKEL)
      epkty = RADIUS * COS(WINKEL)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx 
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly
   END DO
END DO

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
if (testcase == 1) then
!   write(666,'(A)') "wall: 2N, 3N ! Wall"
!   write(666,'(A)') "dn = 1E-2    ! Spacing of first Cell"
else
   write(666,'(A)') "wall: 2N, 3N ! Wall"
   write(666,'(A)') "dn = 1E-6    ! Spacing of first Cell"
end if
close(666)
write(*,*) "done"


end program
