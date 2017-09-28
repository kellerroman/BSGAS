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
integer :: e

call alloc(springs,N_SPRINGS,git % nedge)
do e = 1, git % nedge
   springs(:,e) = 1.0E+0_REAL_KIND / dble(N_SPRINGS) / git % edge_lengths(e)
end do

!write(*,*) maxval(git % wall_edge_dns), minval(git % wall_edge_dns)
end subroutine init_springs

subroutine calc_edge_springs(max_spring,min_spring)
implicit none
real(REAL_KIND), intent(out) :: max_spring
real(REAL_KIND), intent(out) :: min_spring
integer :: e
!integer :: en

call wall_refinment
call edge_streching
call edge_parallel_streching

do e = 1, git % nedge
   git % edge_springs(e) = sum(springs(:,e))
end do
max_spring = maxval( git % edge_springs)
min_spring = minval( git % edge_springs)
!e = 2148
!write(*,*) e,git % edge_lengths(e), springs(:,e), git % point_coords(:,git % edge_points(1,e))&
!         , git % point_coords(:,git % edge_points(2,e))
!en = git % edge_neighbor(1,e)
!write(*,*) en, git % edge_lengths(en), springs(:,en), git % point_coords(:,git % edge_points(1,en))&
!         , git % point_coords(:,git % edge_points(2,en))
!en = git % edge_neighbor(2,e)
!write(*,*) en, git % edge_lengths(en), springs(:,en), git % point_coords(:,git % edge_points(1,en))&
!         , git % point_coords(:,git % edge_points(2,en))
!en = git % edge_neighbor(2,en)
!write(*,*) en, git % edge_lengths(en), springs(:,en), git % point_coords(:,git % edge_points(1,en))&
!         , git % point_coords(:,git % edge_points(2,en))
!en = git % edge_neighbor(2,e)
!write(*,*) git % point_edges(:,git % edge_points(2,en))
end subroutine calc_edge_springs

subroutine wall_refinment
implicit none
integer :: e,i
do i = 1, git % nWallEdge
   e = git % wall_edges(i)
   springs(1,e) = springs(1,e) & 
         * exp(1E-0 * (git % edge_lengths(e) - git % wall_edge_dns(i)))
   springs(1,e) = max(0.0E0_REAL_KIND,springs(1,e))

!   if (git % edge_lengths(e) < 1E-2) then !git % wall_edge_dns(i)) then
!      write(*,'(I5,1X,I7,1X,Es10.3,4(1X,F5.1),3(1X,ES10.3))') &
!                 i,e, git % edge_lengths(e)& 
!                ,0.5D0 * (git % point_refs(:,git % edge_points(1,e)) &
!                         +git % point_refs(:,git % edge_points(2,e)) )&
!                , springs(:,e)
!   end if

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
   springs(3,e) = max(springs(3,e) * exp(1E-2 * fkt), SPRING_MIN)
end do
end subroutine edge_parallel_streching
end module spring
