program gridgen

   use mod_gridgen

implicit none

integer, parameter :: nVar   = 0

real(kind=8) :: length

!!!!!!!!!!!!!!!! nj values 
integer                 :: nj1
integer                 :: nj2

real(kind=8)            :: y1
real(kind=8)            :: y3

integer :: bfs_case
integer :: i,j,k

character(len = 1) :: arg

debug = .true.
WRITE(*,*) "======================================================" &
          ,"==       GRIDGEN for Backward Facing Step           =="  &
          ,"======================================================"
bfs_case = 1
call get_command_argument(1, arg)
if (len_trim(arg) /= 0) then
   read(arg,*) bfs_case
end if

if (bfs_case == 1) then
   nj1 = 100
   nj2 = 200

   y1 = -1.270d-02
   y3 =  1.016d-01
else
   nj1 =  50
   nj2 = 100

   y1 = -6.350d-03
   y3 =  5.080d-02
end if
call set_dimension(1)
call add_block(nj1 + nj2)
call allocate_blocks(nVar)


 k = 1
 j = 1
 length  = abs(y1)
do i = 1,nj2+1
   blocks(1) % xyzs(i,j,k,2) = length/dble(nj2) * dble(i-1) + y1
   blocks(1) % xyzs(i,j,k,1) = 0.0D0
   blocks(1) % xyzs(i,j,k,3) = 0.0d0
end do
 length  = abs(y3)
do i = 1,nj1+1
   blocks(1) % xyzs(i+nj2,j,k,2) = length/dble(nj1) * dble(i-1)
   blocks(1) % xyzs(i,j,k,1)     = 0.0D0
   blocks(1) % xyzs(i,j,k,3)     = 0.0d0
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
