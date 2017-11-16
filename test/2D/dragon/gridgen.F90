program gridgen

   use mod_gridgen

implicit none

real(kind=8), parameter :: a2d = 180.0D0 / 3.1415927D0
integer :: imax
integer :: jmax
integer, parameter :: nVar   = 0

real(kind=8), parameter :: length = 1.0D0
real(kind=8) :: mat(2,2)
real(kind=8) :: winkel = 90.0D0

character(len = 100) :: arg

integer :: i,j,k,b,ir,jr,di,dj

integer :: square_case, nBlock
debug = .true.
write(*,'(A)') "GRID GEN for 2D test grid"
i = 1
square_case = -1
do
   call get_command_argument(i, arg)
   if (len_trim(arg) == 0) exit
   if( i == 1) read(arg,*) square_case
   if( i == 2) read(arg,*) winkel
   i = i+1
end do

call set_dimension(2)
select case(square_case) 
case(1)
   imax =  5
   jmax = 11
   call add_block(imax-1,jmax-1)
   nBlock = 1
case(2)
   imax =  3
   jmax = 11
   call add_block(imax-1,jmax-1)
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
 k = 1
do b = 1, nBlock
   do j = 1, blocks(b) % npkts(2)
      do i = 1,blocks(b) % npkts(1)
         ir = i + di
         jr = j + dj
         blocks(b) % xyzs(i,j,k,1) = length / dble(imax-1) * dble(ir-1)
         blocks(b) % xyzs(i,j,k,2) = length / dble(jmax-1) * dble(jr-1) &
                                   - length / 2.0d0 / dble(imax-1) * dble(ir-1) &
                                   * dble(jr-1) / dble(jmax-1)
         blocks(b) % xyzs(i,j,k,3) = 0.0d0
      end do
   end do
   di = di + blocks(b) % nCells(1)
   if (b == 2) then
      di = 0
      dj = dj + blocks(b) % nCells(2)
   end if
end do

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall "
!select case(square_case) 
!case(1,4)
!write(666,'(A)') "wall: 1S, 1E ! Wall"
!write(666,'(A)') "dn = 4E-2    ! Spacing of first Cell"
!case(2)
!write(666,'(A)') "wall: 1S, 2S, 2E ! Wall"
!write(666,'(A)') "dn = 4E-2    ! Spacing of first Cell"
!case(3)
!write(666,'(A)') "wall: 1S, 2S, 2E, 4E ! Wall"
!write(666,'(A)') "dn = 4E-2    ! Spacing of first Cell"
!case(5)
!write(666,'(A)') "wall: 1W ! Wall"
!write(666,'(A)') "dn = 2.0E-1    ! Spacing of first Cell"
!case default
!   write(*,*) "Case unknown"
!   stop 1
!end select
close(666)

write(*,'(A)') "done"

end program
