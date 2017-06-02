program gridgen

   use mod_gridgen

implicit none

integer :: imax = 11
integer :: jmax =  9
integer, parameter :: nVar   = 0

real(kind=8), parameter :: length = 3.0D0


integer :: i,j,k,b,ir,di

debug = .true.
write(*,'(A)') "GRID GEN for 2D test grid"
call set_dimension(2)
call add_block((imax-1)/2,jmax-1)
call add_block((imax-1)/2,jmax-1)
call allocate_blocks(nVar)

di = 0
 k = 1
do b = 1, 2
   do j = 1, blocks(b) % npkts(2)
      do ir = 1,blocks(b) % npkts(1)
         i = ir + di
         blocks(b) % xyzs(ir,j,k,1) = length/dble(imax-1) * dble(i-1)!* dble(i-1)
         blocks(b) % xyzs(ir,j,k,2) = length/dble(jmax-1) * dble(j-1)
         blocks(b) % xyzs(ir,j,k,3) = 0.0d0
      end do
   end do
   di = di + blocks(b) % nCells(1)
end do

call write_grid()
write(*,'(A)') "done"

end program
