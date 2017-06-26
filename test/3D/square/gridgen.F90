program gridgen

   use mod_gridgen

implicit none

real(kind=8), parameter :: a2d = 180.0D0 / 3.1415927D0
integer :: imax = 11
integer :: jmax = 11
integer :: kmax = 11
integer, parameter :: nVar   = 0

real(kind=8), parameter :: length = 1.0D0
real(kind=8) :: mat(2,2)
real(kind=8) :: temp(2)
real(kind=8) :: winkel = 90.0D0

character(len = 100) :: arg

integer :: i,j,k,b,ir,jr,kr,di,dj,dk

integer :: square_case, nBlock
debug = .true.
write(*,'(A)') "GRID GEN for 3D square test case"
i = 1
square_case = -1
do
   call get_command_argument(i, arg)
   if (len_trim(arg) == 0) exit
   if( i == 1) read(arg,*) square_case
   if( i == 2) read(arg,*) winkel
   i = i+1
end do

call set_dimension(3)
select case(square_case) 
case(1)
   call add_block(imax-1,jmax-1,kmax-1)
   nBlock = 1
case(2)
   call add_block((imax-1)/2,jmax-1,kmax-1)
   call add_block((imax-1)/2,jmax-1,kmax-1)
   nBlock = 2
case(3)
   call add_block((imax-1)/2,(jmax-1)/2,kmax-1)
   call add_block((imax-1)/2,(jmax-1)/2,kmax-1)
   call add_block((imax-1)/2,(jmax-1)/2,kmax-1)
   call add_block((imax-1)/2,(jmax-1)/2,kmax-1)
   nBlock = 4
case(4)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   call add_block((imax-1)/2,(jmax-1)/2,(kmax-1)/2)
   nBlock = 8
case(5)
   call add_block(imax-1,jmax-1,kmax-1)
   nBlock = 1
case default
   write(*,*) "Case unknown"
   stop 1
end select

call allocate_blocks(nVar)

winkel = winkel / a2d
mat(1,1) = + cos(winkel)
mat(2,1) = - sin(winkel)
mat(1,2) = + sin(winkel)
mat(2,2) = + cos(winkel)

di = 0
dj = 0
dk = 0
do b = 1, nBlock
   do k = 1, blocks(b) % nPkts(3)
      do j = 1, blocks(b) % npkts(2)
         do i = 1,blocks(b) % npkts(1)
            ir = i + di
            jr = j + dj
            kr = k + dk
            blocks(b) % xyzs(i,j,k,1) = length/dble(imax-1) * dble(ir-1)!* dble(i-1)
            blocks(b) % xyzs(i,j,k,2) = length/dble(jmax-1) * dble(jr-1)
            blocks(b) % xyzs(i,j,k,3) = length/dble(kmax-1) * dble(kr-1)
            if (square_case == 5) then
               temp = blocks(1) % xyzs(i,j,k,1:2)
               blocks(b) % xyzs(i,j,k,1) = mat(1,1) * temp(1) + mat(2,1) * temp(2)
               blocks(b) % xyzs(i,j,k,2) = mat(1,2) * temp(1) + mat(2,2) * temp(2)
               blocks(b) % xyzs(i,j,k,3) = 0.0d0
            end if
         end do
      end do
   end do
   di = di + blocks(b) % nCells(1)
   if (b == 2 .or. b == 6) then
      di = 0
      dj = dj + blocks(b) % nCells(2)
   end if
   if (b == 4) then
      di = 0
      dj = 0
      dk = dk + blocks(b) % nCells(3)
   end if
end do

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
select case(square_case) 
case(1,5)
write(666,'(A)') "wall: 1S, 1E, 1B! Wall"
case(2)
write(666,'(A)') "wall: 1S, 2S, 2E, 1B,2B ! Wall"
case(3)
write(666,'(A)') "wall: 1S,2S,2E,4E,1F,2F,3F,4F ! Wall"
case(4)
write(666,'(A)') "wall: 1S,2S,5S,6S,2E,4E,6E,8E,1F,2F,3F,4F ! Wall"
case default
   write(*,*) "Case unknown"
   stop 1
end select
write(666,'(A)') "dn = 5E-2    ! Spacing of first Cell"
close(666)

write(*,'(A)') "done"

end program
