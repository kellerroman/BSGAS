program gridgen

   use mod_gridgen

implicit none

integer :: imax = 11
integer :: jmax = 10
integer, parameter :: nVar   = 0

real(kind=8), parameter :: length = 3.0D0


integer :: i,j,k

debug = .true.
write(*,'(A)') "GRID GEN for 2D test grid"
call set_dimension(2)
call add_block(imax-1,jmax-1)
call allocate_blocks(nVar)


 k = 1
do j = 1, jmax
   do i = 1,imax
      blocks(1) % xyzs(i,j,k,1) = length/dble(imax-1) * dble(i-1)!* dble(i-1)
      blocks(1) % xyzs(i,j,k,2) = length/dble(jmax-1) * dble(j-1)
      blocks(1) % xyzs(i,j,k,3) = 0.0d0
   end do
end do

call write_grid()
write(*,'(A)') "done"

end program