module spring
use const
use types
use unstr, only: git
use help_routines, only: alloc
implicit none

enum, bind(C)
   enumerator :: CONTROLTYPE_EXP= 1, CONTROLTYPE_DENOM, CONTROLTYPE_PPID
   ! Spring Control types:
   ! EPX : Expotential, deviation is put into a exponential funktion to change the spring
   ! DENOM : relative deviation is used to correct the spring, underrelaxation is also availble
   ! PPID : Pseudo PID, old value is used to control the spring value
end enum

integer(INT_KIND), parameter :: N_SPRINGS      = 4

real(REAL_KIND), parameter   :: SPRING_MIN     = 1.00E-2_REAL_KIND
real(REAL_KIND), parameter   :: SPRING_SMOOTH_MAX = SPRING_MIN * 1.0E+12_REAL_KIND

integer(INT_KIND)            :: spring_control_type_wall   = CONTROLTYPE_DENOM
integer(INT_KIND)            :: spring_control_type_strech = CONTROLTYPE_DENOM
integer(INT_KIND)            :: spring_control_type_para   = CONTROLTYPE_DENOM
integer(INT_KIND)            :: spring_control_type_smooth = CONTROLTYPE_DENOM
integer(INT_KIND)            :: smooth_springs = 0

real(REAL_KIND)              :: cell_inc
real(REAL_KIND)              :: cell_parallel_inc
real(REAL_KIND)              :: faktor_wall
real(REAL_KIND)              :: faktor_strech
real(REAL_KIND)              :: faktor_para
real(REAL_KIND)              :: faktor_smooth
real(REAL_KIND)              :: spring_max
real(REAL_KIND)              :: relax_fkt
real(REAL_KIND)              :: relax_inv

integer(INT_KIND), allocatable :: driving_forces(:)

real(REAL_KIND), allocatable :: springs(:,:)
real(REAL_KIND), allocatable :: edge_values(:,:)
real(REAL_KIND), allocatable :: spring_old(:)
real(REAL_KIND), allocatable :: spring_res(:)

contains

subroutine init_springs
implicit none
integer :: e,i

call alloc(springs,N_SPRINGS,git % nedge)
call alloc(edge_values,N_SPRINGS,git % nedge)
call alloc(driving_forces,git % nedge)
call alloc(spring_old,git % nedge)
call alloc(spring_res,git % nedge)
do e = 1, git % nedge
   springs(2:,e) = SPRING_MIN / git % edge_lengths(e)  ! Length is taken to put grid at equilibrium at start up
   springs(1 ,e) = SPRING_MIN ! Wall-Edges do noz have automatic decrease so start with very small values
end do

! Set WallEdges to the VAlue of the other springs to see an effect in the first iteration
do i = 1, git % nWallEdge
   e = git % wall_edges(i)
   springs(1,e) = springs(2,e)
end do
end subroutine init_springs

subroutine calc_edge_springs(max_spring,min_spring,max_spring_res,sum_spring_res)
!use control, only: iter
implicit none
real(REAL_KIND), intent(out) :: max_spring
real(REAL_KIND), intent(out) :: min_spring
real(REAL_KIND), intent(out) :: max_spring_res
real(REAL_KIND), intent(out) :: sum_spring_res
integer :: e,n,ne

!if (mod(iter,spring_intervall) /= 0) return
relax_inv = 1.00E+00_REAL_KIND - relax_fkt

call wall_refinement
call edge_streching
call edge_parallel_streching
call edge_smoothing

!$OMP PARALLEL DO
do e = 1, git % nedge
   spring_old(e) = git % edge_springs(e)
   !springs(:,e)          = min(spring_max,springs(:,e))
   git % edge_springs(e) = maxval(           springs(:,e))
   driving_forces(e)     = maxloc(           springs(:,e),1)
end do
!$OMP END PARALLEL DO
! smoothing
if (smooth_springs == 1) then
!$OMP PARALLEL DO PRIVATE(ne,n)
   do e = 1, git % nedge
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         springs(2,ne) = max(git % edge_springs(e) / 1.3D0, springs(2,ne))
         springs(2,e) = max(git % edge_springs(ne) / 1.3D0, springs(2,e))
      end do
      do n = 1, git % edge_nparallel(e)
         ne = git % edge_parallel(n,e)
         springs(3,ne) = max(git % edge_springs(e) / 1.3D0, springs(3,ne))
         springs(3,e) = max(git % edge_springs(ne) / 1.3D0, springs(3,e))
      end do
   end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(ne,n)
   do e =  git % nedge,1,-1
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         springs(2,ne) = max(git % edge_springs(e) / 1.3D0, springs(2,ne))
         springs(2,e) = max(git % edge_springs(ne) / 1.3D0, springs(2,e))
         git % edge_springs(ne) = maxval(           springs(:,ne))
         driving_forces(ne)     = maxloc(           springs(:,ne),1)
         git % edge_springs(e)  = maxval(           springs(:,e))
         driving_forces(e)      = maxloc(           springs(:,e),1)
      end do
      do n = 1, git % edge_nparallel(e)
         ne = git % edge_parallel(n,e)
         springs(3,ne) = max(git % edge_springs(e) / 1.3D0, springs(3,ne))
         springs(3,e) = max(git % edge_springs(ne) / 1.3D0, springs(3,e))
         git % edge_springs(ne) = maxval(           springs(:,ne))
         driving_forces(ne)     = maxloc(           springs(:,ne),1)
         git % edge_springs(e)  = maxval(           springs(:,e))
         driving_forces(e)      = maxloc(           springs(:,e),1)
      end do
   end do
!$OMP END PARALLEL DO
else if (smooth_springs == 2) then
!$OMP PARALLEL DO PRIVATE(ne,n)
   do e = 1, git % nedge
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         git % edge_springs(ne) = max(git % edge_springs(e) / 1.3D0, git % edge_springs(ne))
         git % edge_springs(e) = max(git % edge_springs(ne) / 1.3D0, git % edge_springs(e))
      end do
      do n = 1, git % edge_nparallel(e)
         ne = git % edge_parallel(n,e)
         git % edge_springs(ne) = max(git % edge_springs(e) / 1.3D0, git % edge_springs(ne))
         git % edge_springs(e) = max(git % edge_springs(ne) / 1.3D0, git % edge_springs(e))
      end do
   end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(ne,n)
   do e =  git % nedge,1,-1
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         git % edge_springs(ne) = max(git % edge_springs(e) / 1.3D0, git % edge_springs(ne))
         git % edge_springs(e) = max(git % edge_springs(ne) / 1.3D0, git % edge_springs(e))
      end do
      do n = 1, git % edge_nparallel(e)
         ne = git % edge_parallel(n,e)
         git % edge_springs(ne) = max(git % edge_springs(e) / 1.3D0, git % edge_springs(ne))
         git % edge_springs(e) = max(git % edge_springs(ne) / 1.3D0, git % edge_springs(e))
      end do
   end do
!$OMP END PARALLEL DO
end if !Smoothing
!$OMP PARALLEL DO
do e = 1, git % nedge
   spring_res(e) = abs(spring_old(e) - git % edge_springs(e))
end do
!$OMP END PARALLEL DO
max_spring = maxval( git % edge_springs)
min_spring = minval( git % edge_springs)
max_spring_res = maxval(spring_res)
sum_spring_res = sum(spring_res)
end subroutine calc_edge_springs

subroutine wall_refinement
implicit none
integer :: e,i
real(REAL_KIND) :: fkt
if (spring_control_type_wall == CONTROLTYPE_EXP) then
!$OMP PARALLEL DO PRIVATE(e)
   do i = 1, git % nWallEdge
      e = git % wall_edges(i)
      fkt = exp(faktor_wall * (git % edge_lengths(e) - git % wall_edge_dns(i)))
      edge_values(1,e) = git % edge_lengths(e) / git % wall_edge_dns(i) - 1.0D0
      springs(1,e) = max(SPRING_MIN, springs(1,e) * fkt)
   end do
!$OMP END PARALLEL DO
else
!$OMP PARALLEL DO PRIVATE(e)
   do i = 1, git % nWallEdge
      e = git % wall_edges(i)
      fkt = git % edge_lengths(e) / git % wall_edge_dns(i)
      edge_values(1,e) = fkt
      fkt = fkt * relax_fkt + relax_inv
      springs(1,e) = max(SPRING_MIN, springs(1,e) * fkt)
   end do
!$OMP END PARALLEL DO
end if
end subroutine wall_refinement

subroutine edge_streching
implicit none
integer :: e
integer :: ne  !neighbor edge
integer :: n

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
real(REAL_KIND) :: delta

if (spring_control_type_strech == CONTROLTYPE_EXP) then
!$OMP PARALLEL DO PRIVATE(el,nel,fkt,ne,n,delta)
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
      springs(2,e) = max(springs(2,e) * delta, SPRING_MIN)
   end do
!$OMP END PARALLEL DO
else
!$OMP PARALLEL DO PRIVATE(el,nel,fkt,ne,n,delta)
   do e = 1, git % nedge
      el = git % edge_lengths(e)
      nel = 1E10_REAL_KIND
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         nel = min(nel, git % edge_lengths(ne))
      end do
      fkt = el / nel / cell_inc
      edge_values(2,e) = fkt
      fkt = fkt * relax_fkt + relax_inv
      springs(2,e) = max(springs(2,e) * fkt, SPRING_MIN)
   end do
!$OMP END PARALLEL DO
end if
end subroutine edge_streching

subroutine edge_parallel_streching
implicit none
integer :: e
integer :: ne  !neighbor edge
integer :: n

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
if (spring_control_type_para == CONTROLTYPE_EXP) then
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
      !springs(3,e) = max(springs(3,e) * fkt, SPRING_MIN)
   end do
!$OMP END PARALLEL DO
else
!$OMP PARALLEL DO PRIVATE(el,nel,fkt)
   do e = 1, git % nedge
      el = git % edge_lengths(e)
      nel = 1E10_REAL_KIND
      do n = 1, git % edge_nparallel(e)
         ne = git % edge_parallel(n,e)
         nel = min(nel,git % edge_lengths(ne))
      end do
      fkt = el / nel / cell_parallel_inc
      edge_values(3,e) = fkt
      fkt = fkt * relax_fkt + relax_inv
      springs(3,e) = max(springs(3,e) * fkt, SPRING_MIN)
   end do
!$OMP END PARALLEL DO
end if
end subroutine edge_parallel_streching

subroutine edge_smoothing
implicit none
integer :: e
integer :: ne  !neighbor edge
integer :: n

real(REAL_KIND) :: el,nel
real(REAL_KIND) :: fkt
real(REAL_KIND) :: delta

if (spring_control_type_smooth == CONTROLTYPE_EXP) then
!$OMP PARALLEL DO PRIVATE(el,nel,fkt,ne,n,delta)
   do e = 1, git % nedge
      write(*,*) "not implemented yet", __FILE__,__LINE__
      stop 1
      el = git % edge_lengths(e)
      nel = 1E10_REAL_KIND
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         nel = min(nel, git % edge_lengths(ne))
      end do
      fkt = el / nel - 1.0E0_REAL_KIND
      edge_values(4,e) = fkt
      delta = exp(faktor_smooth * fkt)
      springs(4,e) = max(springs(4,e) * delta, SPRING_MIN)
   end do
!$OMP END PARALLEL DO
else if (spring_control_type_smooth == CONTROLTYPE_PPID) then
!$OMP PARALLEL DO PRIVATE(el,nel,fkt,ne,n,delta)
   do e = 1, git % nedge
      el = git % edge_lengths(e)
      nel = 1E10_REAL_KIND
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         nel = min(nel, git % edge_lengths(ne))
      end do
      fkt =  - nel + 2.0d0 * el - edge_values(4,e)
      edge_values(4,e) = el
      delta = faktor_smooth * fkt
      springs(4,e) = min(max(springs(4,e) * delta, SPRING_MIN),SPRING_SMOOTH_MAX*el)
   end do
!$OMP END PARALLEL DO
else
!$OMP PARALLEL DO PRIVATE(el,nel,fkt,ne,n,delta)
   do e = 1, git % nedge
      el = git % edge_lengths(e)
      nel = 1E10_REAL_KIND
      do n = 1, git % edge_nneighbor(e)
         ne = git % edge_neighbor(n,e)
         nel = min(nel, git % edge_lengths(ne))
      end do
      fkt = el / nel
      edge_values(4,e) = fkt
      fkt = fkt * relax_fkt + relax_inv
      springs(4,e) = min(max(springs(4,e) * fkt, SPRING_MIN),SPRING_SMOOTH_MAX*el)
      !springs(4,e) = SPRING_MIN
   end do
!$OMP END PARALLEL DO
end if
end subroutine edge_smoothing
end module spring
