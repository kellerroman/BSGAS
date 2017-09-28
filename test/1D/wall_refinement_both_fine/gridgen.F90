program gridgen

   use mod_gridgen

implicit none

integer :: imax = 200
integer, parameter :: nVar   = 0

real(kind=8), parameter :: length = 1.016D-1


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
open(666,file="bc.cfg")
write(666,'(A)') "wall: 1N,1S!,2S,2W,3N ! Wall"
write(666,'(A)') "dn = 1E-6    ! Spacing of first Cell"
close(666)

WRITE(*,*) "======================================================" &
          ,"FINISHED"  &
          ,"======================================================"

end program
