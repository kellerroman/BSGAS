module spring
use const
use types
use unstr, only: git
use help_routines, only: alloc
implicit none
integer(INT_KIND), parameter :: N_SPRINGS      = 3
real(REAL_KIND)              :: faktor_wall
real(REAL_KIND)              :: faktor_strech
real(REAL_KIND)              :: faktor_para
real(REAL_KIND)              :: spring_max
real(REAL_KIND), parameter   :: SPRING_MIN     = 1.00E-10_REAL_KIND
real(REAL_KIND), parameter   :: SPRING_INC     = 1.05E-00_REAL_KIND
real(REAL_KIND), parameter   :: INV_SPRING_INC = 1.0E0_REAL_KIND / SPRING_INC
real(REAL_KIND), allocatable :: springs(:,:)
real(REAL_KIND), allocatable :: edge_values(:,:)

real(REAL_KIND) :: cell_inc
real(REAL_KIND) :: cell_parallel_inc
contains

subroutine init_springs
implicit none
integer :: e

call alloc(springs,N_SPRINGS,git % nedge)
call alloc(edge_values,N_SPRINGS,git % nedge)
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

call wall_refinement
call edge_streching
call edge_parallel_streching

!$OMP PARALLEL DO
do e = 1, git % nedge
   springs(:,e)          = min(spring_max,springs(:,e))
   git % edge_springs(e) = sum(           springs(:,e))
end do
!$OMP END PARALLEL DO
max_spring = maxval( git % edge_springs)
min_spring = minval( git % edge_springs)
end subroutine calc_edge_springs

subroutine wall_refinement
implicit none
integer :: e,i
!$OMP PARALLEL DO PRIVATE(e)
do i = 1, git % nWallEdge
   e = git % wall_edges(i)
   springs(1,e) = springs(1,e) & 
         * exp(faktor_wall * (git % edge_lengths(e) - git % wall_edge_dns(i)))
   edge_values(1,e) = git % edge_lengths(e) - git % wall_edge_dns(i)
   springs(1,e) = max(0.0E0_REAL_KIND ,springs(1,e))
end do
!$OMP END PARALLEL DO
end subroutine wall_refinement

subroutine edge_streching
implicit none
integer :: e
integer :: ne  !neighbor edge
integer :: n

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
real(REAL_KIND) :: delta
!$OMP PARALLEL DO PRIVATE(el,nel,fkt,delta)
do e = 1, git % nedge
   el = git % edge_lengths(e)
   nel = 1E10_REAL_KIND
   do n = 1, git % edge_nneighbor(e)
      ne = git % edge_neighbor(n,e)
      nel = min(nel, git % edge_lengths(ne))
   end do
   fkt = el / nel - cell_inc
   edge_values(2,e) = fkt
   delta = exp(faktor_strech * fkt)
   !delta = max (delta,INV_SPRING_INC )
   !delta = min (delta,    SPRING_INC )
   springs(2,e) = max(springs(2,e) * delta, SPRING_MIN)
end do
!$OMP END PARALLEL DO
end subroutine edge_streching

subroutine edge_parallel_streching
implicit none
integer :: e
integer :: ne  !neighbor edge
integer :: n

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
!$OMP PARALLEL DO PRIVATE(el,nel,fkt)
do e = 1, git % nedge
   el = git % edge_lengths(e)
   nel = 1E10_REAL_KIND
   do n = 1, git % edge_nparallel(e)
      ne = git % edge_parallel(n,e)
      nel = min(nel,git % edge_lengths(ne))
   end do
   fkt = el / nel
   fkt = fkt - cell_parallel_inc
   edge_values(3,e) = fkt
   springs(3,e) = max(springs(3,e) * exp(faktor_para * fkt), SPRING_MIN)
end do
!$OMP END PARALLEL DO
end subroutine edge_parallel_streching

end module spring
