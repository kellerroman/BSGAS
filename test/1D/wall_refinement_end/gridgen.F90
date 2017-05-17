program gridgen

   use mod_gridgen

implicit none

integer :: imax = 11
integer, parameter :: nVar   = 0

real(kind=8), parameter :: length = 3.0D0


integer :: i,j,k

debug = .true.
write(*,'(A)') "GRID GEN for 1D test grid"
call set_dimension(1)
call add_block(imax-1)
call allocate_blocks(nVar)


 k = 1
 j = 1
do i = 1,imax
   blocks(1) % xyzs(i,j,k,1) = length/dble(imax-1) * dble(i-1)!* dble(i-1)
   blocks(1) % xyzs(i,j,k,2) = 0.0D0
   blocks(1) % xyzs(i,j,k,3) = 0.0d0
end do

call write_grid()
write(*,'(A)') "done"

end program
