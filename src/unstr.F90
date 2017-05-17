module unstr
use const
use types
use help_routines, only: alloc
implicit none
type(t_unstr) :: git
contains
subroutine strukt2unstr(blocks)
implicit none
type(t_block) :: blocks(:)

integer :: i,j,k,np,ne,b

integer :: p, p1, pe

write(*,*) "Block structured number of Blocks:",ubound(blocks)

b = 1

if (blocks(b) % nPoints(2) == 1) then
   git % dimension = 1
else if  ( blocks(b) % nPoints(3) == 1) then
   git % dimension = 2
else
   git % dimension = 3
end if

write(*,*) "Grid is of Dimension:", git % dimension

git % npoint = product(blocks(b) % nPoints)

if (git % dimension == 1) then
   git % nedge = git % npoint - 1
end if

write(*,*) "Number of Points & Edges in unstr Grid:",git % npoint, git % nedge

call alloc(git % point_coords    , 3, git % npoint)
call alloc(git % point_refs      , 4, git % npoint)
call alloc(git % point_nedges       , git % npoint)
call alloc(git % point_edges     , 6, git % npoint)
call alloc(git % point_edge_signs, 6, git % npoint)
call alloc(git % point_forces    , 3, git % npoint)
call alloc(git % edge_lengths       , git % nedge)
call alloc(git % edge_points     , 2, git % nedge)
call alloc(git % edge_springs       , git % nedge)
call alloc(git % edge_vectors    , 3, git % nedge)
call alloc(git % edge_forces     , 3, git % nedge)

git % point_forces = 0.0e0_REAL_KIND
git % point_nedges = 0
git % point_forces = 0.0e0_REAL_KIND

np = 0
ne = 0
do k = 1, blocks(b) % nPoints(3)
   do j = 1, blocks(b) % nPoints(2)
      do i = 1, blocks(b) % nPoints(1)
         np = np + 1
         git % point_coords(:,np) = blocks(b) % coords(i,j,k,:)
         git % point_refs(1,np) = b
         git % point_refs(2,np) = i
         git % point_refs(3,np) = j
         git % point_refs(4,np) = k
         if (i > 1) then
            ne = ne + 1
            ! Add points to edge and edge to points
            do p = 1, 2
               p1 = np - 2 + p ! id of point, since it is i-direction it is np-1 & np
               ! Add point to edge
               git % edge_points(p,ne) = p1
               ! Add edge to the first point
               !Increase points edge count
               pe  = git % point_nedges(p1) + 1
               git % point_nedges(p1) = pe
               ! add edge to point at current edge count
               git % point_edges(pe,p1) = ne
               ! add sign of the resulting edge force, for p2 -> -1 p1 -> 1
               if (p == 1) then
                  git % point_edge_signs(pe,p1) =  1.0E0_REAL_KIND
               else
                  git % point_edge_signs(pe,p1) = -1.0E0_REAL_KIND
               end if

            end do
         end if
      end do
   end do
end do

do np = 1, git % nPoint
   pe = git % point_nedges(np)
   write(*,*) np, git % point_coords(:,np),pe , git % point_edges(1:pe,np)
end do

do ne = 1, git % nedge
   write(*,*) ne, git % edge_points(:,ne)
end do

git % edge_springs = 1.0e0_REAL_KIND

git % edge_lengths = -1
git % edge_vectors = -1
git % edge_forces = -1

end subroutine strukt2unstr

subroutine calc_edge_length()
implicit none
integer :: ne,p1,p2
real(REAL_KIND) :: x1,x2,y1,y2,z1,z2
real(REAL_KIND) :: dx,dy,dz
do ne = 1, git % nedge
   p1 = git % edge_points(1,ne)
   p2 = git % edge_points(2,ne)
   x1 = git % point_coords(1,p1)
   x2 = git % point_coords(1,p2)
   y1 = git % point_coords(2,p1)
   y2 = git % point_coords(2,p2)
   z1 = git % point_coords(3,p1)
   z2 = git % point_coords(3,p2)
   dx = x2-x1
   dy = y2-y1
   dz = z2-z1
   git % edge_vectors(1,ne) = dx
   git % edge_vectors(2,ne) = dy
   git % edge_vectors(3,ne) = dz
   git % edge_lengths(ne) = sqrt( dx * dx &
                                + dy * dy &
                                + dz * dz )
   !write(*,*) ne, git % edge_lengths(ne), git % edge_vectors(:,ne)
end do
end subroutine calc_edge_length
subroutine calc_edge_forces(max_f)
implicit none
real(REAL_KIND), intent(out) :: max_f

integer :: ne

max_f = 0.0E0_REAL_KIND
do ne = 1, git % nedge
   git % edge_forces(:,ne) = git % edge_springs(ne) * git % edge_vectors(:,ne)
   max_f = max(max_f,sqrt( git % edge_forces(1,ne) * git % edge_forces(1,ne) &
                         + git % edge_forces(2,ne) * git % edge_forces(2,ne) &
                         + git % edge_forces(3,ne) * git % edge_forces(3,ne) ))
end do
end subroutine calc_edge_forces

subroutine calc_point_forces(max_f)
implicit none
real(REAL_KIND), intent(out) :: max_f

integer :: np, ne, edge
real(REAL_KIND) :: tmp(3)
max_f = 0.0E0_REAL_KIND
do np = 2, git % nPoint-1
    tmp = 0.0E0_REAL_KIND
   do ne = 1, git % point_nedges(np)
      edge = git % point_edges(ne,np)
      tmp = tmp + git % edge_forces(:,edge) * git % point_edge_signs(ne,np)
   end do
   git % point_forces(:,np) = tmp
   max_f = max(max_f,sqrt(tmp(1)*tmp(1)+tmp(2)*tmp(2)+tmp(3)*tmp(3)))
   !write(*,*) np,git % point_nedges(np), git % point_forces(:,np)
end do
end subroutine calc_point_forces
subroutine move_points(max_f)
implicit none
real(REAL_KIND), parameter :: POINT_WEIGTH =  1.5E0_REAL_KIND
real(REAL_KIND), intent(in) :: max_f
real(REAL_KIND) :: min_l, fk

integer :: np

min_l = minval(git % edge_lengths)
fk = min(1.0E0_REAL_KIND,min_l / max_f)
!write(*,*) "FK",min_l,max_f,fk
do np = 2, git % npoint-1
   git % point_coords(:,np) = git % point_coords(:,np) + git % point_forces(:,np) / POINT_WEIGTH * fk

end do
end subroutine move_points
end module unstr
