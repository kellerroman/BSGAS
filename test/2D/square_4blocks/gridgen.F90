program gridgen

   use mod_gridgen

implicit none

integer :: imax = 11
integer :: jmax = 11
integer, parameter :: nVar   = 0

real(kind=8), parameter :: length = 3.0D0


integer :: i,j,k,b,ir,di,jr,dj

debug = .true.
write(*,'(A)') "GRID GEN for 2D test grid"
call set_dimension(2)
call add_block((imax-1)/2,(jmax-1)/2)
call add_block((imax-1)/2,(jmax-1)/2)
call add_block((imax-1)/2,(jmax-1)/2)
call add_block((imax-1)/2,(jmax-1)/2)
call allocate_blocks(nVar)

di = 0
dj = 0
 k = 1
do b = 1, 4
   do j = 1, blocks(b) % npkts(2)
      do i = 1,blocks(b) % npkts(1)
         ir = i + di
         jr = j + dj
         blocks(b) % xyzs(i,j,k,1) = length/dble(imax-1) * dble(ir-1)!* dble(i-1)
         blocks(b) % xyzs(i,j,k,2) = length/dble(jmax-1) * dble(jr-1)
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
write(*,'(A)') "done"

end program
