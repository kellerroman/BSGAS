module spring
use const
use types
use unstr, only: git
use help_routines, only: alloc
implicit none
integer(INT_KIND), parameter :: N_SPRINGS = 2
real(REAL_KIND), parameter :: dw_soll = 1.0E-1_REAL_KIND
real(REAL_KIND), allocatable :: springs(:,:)
contains

subroutine init_springs
implicit none

call alloc(springs,N_SPRINGS,git % nedge)
springs = 1.0E+0_REAL_KIND / dble(N_SPRINGS)

end subroutine init_springs

subroutine calc_edge_springs
implicit none
integer :: e

call wall_refinment
call edge_streching

do e = 1, git % nedge
   git % edge_springs(e) = sum(springs(:,e))
end do
end subroutine calc_edge_springs

subroutine wall_refinment
implicit none
integer :: e
e = git % nedge
springs(1,e) = springs(1,e) & 
      * exp(1E-2 * (git % edge_lengths(e) - dw_soll))
springs(1,e) = max(0.0E0_REAL_KIND,springs(1,e))

!write(*,*)  springs(1,e), git % edge_lengths(e),exp(1E-2 * (git % edge_lengths(e) - dw_soll))
end subroutine wall_refinment


subroutine edge_streching
implicit none
real(REAL_KIND), parameter :: dl_max = 1.2E+0_REAL_KIND
integer :: e
integer :: ne  !neighbor edge

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
do e = 1, git % nedge
   el = git % edge_lengths(e)
   nel = 1E10_REAL_KIND
   if (e > 1) then
      ne = e - 1
      nel = min(nel,git % edge_lengths(ne))
   end if
   if (e < git % nedge) then
      ne = e + 1
      nel = min(nel,git % edge_lengths(ne))
   end if
   fkt = el / nel
   fkt = fkt - dl_max
   springs(2,e) = max(springs(2,e) * exp(1E-2 * fkt), 0.0E+0_REAL_KIND)
end do
end subroutine edge_streching
end module spring
