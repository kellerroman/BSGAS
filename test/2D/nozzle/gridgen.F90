program gridgen_nozzle_simple
   use mod_gridgen
implicit none

real(kind = 8), parameter :: POS_THROAT = 0.3

integer, parameter :: NI   = 20
integer, parameter :: NI_THROAT = int(NI * POS_THROAT)

integer, parameter :: NJ = 10

real(kind = 8), parameter :: Y_IN     = 2.0D0
real(kind = 8), parameter :: Y_THROAT = 1.0D0
real(kind = 8), parameter :: Y_END    = 3.0D0
real(kind = 8), parameter :: X_END    = 7.0D0
real(kind = 8), parameter :: X_THROAT = X_END * POS_THROAT

integer :: i ,j, k,b
real(kind = 8) ::  y_max
write(*,*) "gridgen for a simple nozzle 29.01.2015"

WRITE(*,'(A,3(1X,A,":",1X,I0))') "GITTERDIMENSIONEN:", "NI",NI,"NI_THROAT",NI_THROAT,"NJ",NJ

WRITE(*,*) "INLET : Y:", Y_IN
WRITE(*,*) "THROAT: X:", X_THROAT,"Y:",Y_THROAT
WRITE(*,*) "EXIT  : X:", X_END   ,"Y:",Y_END

call set_dimension(2)
call add_block(NI,NJ)
call allocate_blocks(0)
k = 1
b = 1
do j = 0, nj
   do i = 0,ni
      if (i < NI_THROAT) then
         y_max = Y_IN - (Y_IN - Y_THROAT) * dble(i) / dble(NI_THROAT)

      else
         y_max = Y_THROAT + (Y_END-Y_THROAT) * dble(i-NI_THROAT) / dble(ni-NI_THROAT)
      end if
      blocks(b) % xyzs(i+1,j+1,k,1) = X_END * dble(i)/dble(ni)
      blocks(b) % xyzs(i+1,j+1,k,2) = Y_MAX * dble(j)/dble(nj)
   end do
end do

call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
write(666,'(A)') "wall: 1N ! Wall"
write(666,'(A)') "dn = 5E-2    ! Spacing of first Cell"
close(666)
write(*,*) "done"


end program gridgen_nozzle_simple
