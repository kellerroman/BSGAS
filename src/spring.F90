module spring
use const
use types
use unstr, only: git
use help_routines, only: alloc
implicit none
integer(INT_KIND), parameter :: N_SPRINGS = 3
real(REAL_KIND),   parameter :: SPRING_MIN = 1.0E-10_REAL_KIND
real(REAL_KIND), allocatable :: springs(:,:)

real(REAL_KIND) :: cell_inc
real(REAL_KIND) :: cell_parallel_inc
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
call edge_parallel_streching

do e = 1, git % nedge
   git % edge_springs(e) = sum(springs(:,e))
end do
end subroutine calc_edge_springs

subroutine wall_refinment
implicit none
integer :: e,i
do i = 1, git % nWallEdge
   e = git % wall_edges(i)
   springs(1,e) = springs(1,e) & 
         * exp(1E-2 * (git % edge_lengths(e) - git % wall_edge_dns(i)))
   springs(1,e) = max(0.0E0_REAL_KIND,springs(1,e))

!   write(*,*)  e,springs(1,e), git % edge_lengths(e),exp(1E-2 * (git % edge_lengths(e) - dw_soll))
end do
end subroutine wall_refinment


subroutine edge_streching
implicit none
integer :: e
integer :: ne  !neighbor edge
integer :: n

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
do e = 1, git % nedge
   el = git % edge_lengths(e)
   nel = 1E10_REAL_KIND
   do n = 1, git % edge_nneighbor(e)
      ne = git % edge_neighbor(n,e)
      nel = min(nel,git % edge_lengths(ne))
   end do
   fkt = el / nel
   fkt = fkt - cell_inc
   springs(2,e) = max(springs(2,e) * exp(1E-2 * fkt),SPRING_MIN)
end do
end subroutine edge_streching
subroutine edge_parallel_streching
implicit none
integer :: e
integer :: ne  !neighbor edge
integer :: n

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
do e = 1, git % nedge
   el = git % edge_lengths(e)
   nel = 1E10_REAL_KIND
   do n = 1, git % edge_nparallel(e)
      ne = git % edge_parallel(n,e)
      nel = min(nel,git % edge_lengths(ne))
   end do
   fkt = el / nel
   fkt = fkt - cell_parallel_inc
   !springs(3,e) = max(springs(3,e) * exp(1E-2 * fkt), SPRING_MIN)
end do
end subroutine edge_parallel_streching
end module spring
