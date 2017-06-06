module unstr
use const
use types
use help_routines, only: alloc, vec_common
implicit none
type(t_unstr) :: git
contains
subroutine strukt2unstr(blocks)
implicit none
type(t_block) :: blocks(:)

integer :: nBlock
integer :: i,j,k,np,e,b,nb,nbe,nbn,nben,nbne


integer :: p, p1, pe, eop, e2, p2, ps

integer :: nne, ne

logical :: edge_exists,do_connect

nBlock = ubound(blocks,1)

write(*,*) "Block structured number of Blocks:",nBlock

b = 1
git % nPoint = 0
git % nEdge  = 0

if (blocks(b) % nPoints(2) == 1) then
   git % dimension = 1
else if  ( blocks(b) % nPoints(3) == 1) then
   git % dimension = 2
else
   git % dimension = 3
end if

write(*,*) "Grid is of Dimension:", git % dimension

do b = 1, nBlock
   i = blocks(b) % nPoints(1)
   j = blocks(b) % nPoints(2)
   k = blocks(b) % nPoints(3)

   if (blocks(b) % boundary_cond(WEST) % bc_type > 0) then
      i = i - 1
   end if
   if (blocks(b) % boundary_cond(SOUTH) % bc_type > 0) then
      j = j - 1
   end if
   if (git % dimension ==3) then
      if (blocks(b) % boundary_cond(FRONT) % bc_type > 0) then
         k = k - 1
      end if
   end if

   git % nPoint = git % nPoint + i * j * k
   write(*,*) i,j,k

   if (git % dimension == 1) then
      git % nEdge = git % npoint - 1
   else if (git % dimension == 2) then
      git % nEdge = git % nEdge &
                  + blocks(b) % nCells(1) * j &
                  + blocks(b) % nCells(2) * i
   end if


end do

write(*,*) "Number of Points & Edges in unstr Grid:",git % npoint, git % nedge

call alloc(git % point_coords          , 3, git % npoint)
call alloc(git % point_refs            , 4, git % npoint)
call alloc(git % point_nedges             , git % npoint)
call alloc(git % point_edges           , 6, git % npoint)
call alloc(git % point_edge_signs      , 6, git % npoint)
call alloc(git % point_forces          , 3, git % npoint)
call alloc(git % point_move_rest          , git % npoint)
call alloc(git % point_move_rest_vector, 3, git % npoint)
call alloc(git % edge_lengths             , git % nedge)
call alloc(git % edge_points           , 2, git % nedge)
call alloc(git % edge_springs             , git % nedge)
call alloc(git % edge_vectors          , 3, git % nedge)
call alloc(git % edge_forces           , 3, git % nedge)
call alloc(git % edge_nneighbor           , git % nedge)
call alloc(git % edge_neighbor         , 2, git % nedge)


git % point_forces = 0.0e0_REAL_KIND
git % point_nedges = 0
git % point_forces = 0.0e0_REAL_KIND

git % edge_nneighbor = 0
git % edge_neighbor  = -1
np = 0
e = 0
do b = 1, nBlock
   do k = 1, blocks(b) % nPoints(3)
      do j = 1, blocks(b) % nPoints(2)
         do i = 1, blocks(b) % nPoints(1)
            if (blocks(b) % refs(i,j,k) == -1) then
   !====================================================================================================
   !==========================      CREATE POINT    ====================================================
   !====================================================================================================
               ! creation of Point
               np = np + 1
               git % point_coords(:,np) = blocks(b) % coords(i,j,k,:)
               git % point_refs(1,np) = b
               git % point_refs(2,np) = i
               git % point_refs(3,np) = j
               git % point_refs(4,np) = k
               blocks(b) % refs(i,j,k) = np

   !====================================================================================================
   !==========================   TRANSFER POINT TO OTHER BLOCKS ========================================
   !====================================================================================================
               ! Check if this point is also a point in another grid
               !!! CORNER POINT:
               if (i == blocks(b) % nPoints(1) .and. j == blocks(b) % nPoints(2)) then
                  nbe = blocks(b) % boundary_cond(EAST) % bc_type   ! Neighbor EAST
                  nbn = blocks(b) % boundary_cond(NORTH) % bc_type   ! Neighbor NORTH
                  ! Both sides have a connecting Block, Ckeck if the both have a
                  ! diagonal neighbor
                  ! If this is the same point, also transfer point there
                  if (nbe > 0 .and. nbn > 0) then
                     if   (      blocks(b) % boundary_cond(EAST) % permutation == 1  & 
                           .and. blocks(b) % boundary_cond(NORTH) % permutation == 1) then
                        nben = blocks(nbe) % boundary_cond(NORTH) % bc_type
                        nbne = blocks(nbn) % boundary_cond(EAST) % bc_type
                        if (nben > 0 .and. nbne > 0 .and. nben == nbne) then
                           if   (      blocks(nbn) % boundary_cond(EAST) % permutation == 1  & 
                                 .and. blocks(nbe) % boundary_cond(NORTH) % permutation == 1) then
                              write(*,*) "Transfering Point to diagonal Block:",b,i,j,k, " to ", nben,1,1,k
                              blocks(nben) % refs(1,1,k) = np

                           else
                              write(*,*) "Error in UNSTR: Permutation diagonal noch nicht implementiert" &
                                        ,__LINE__,__FILE__
                              write(*,*) "Point:",np,"REF:",b,i,j,k,"to",nb,"PERM:", blocks(nbn) % boundary_cond(EAST) % permutation
                              stop 1
                           end if
                        end if
                     else if   ( blocks(b) % boundary_cond(EAST) % permutation == 2  & 
                           .and. blocks(b) % boundary_cond(NORTH) % permutation == 1) then
                        nben = blocks(nbe) % boundary_cond(WEST) % bc_type
                        nbne = blocks(nbn) % boundary_cond(EAST) % bc_type
                        if (nbe == nbne .and. nbn == nben) then
                           write(*,*) "O-Grid Detected"
                        else if (nben > 0 .and. nbne > 0 .and. nben == nbne) then
                           if   (      blocks(nbn) % boundary_cond(EAST) % permutation == 1  & 
                                 .and. blocks(nbe) % boundary_cond(NORTH) % permutation == 1) then
                              write(*,*) "Transfering Point to diagonal Block:",b,i,j,k, " to ", nben,1,1,k
                              blocks(nben) % refs(1,1,k) = np

                           else
                              write(*,*) "Error in UNSTR: Permutation diagonal noch nicht implementiert" &
                                        ,__LINE__,__FILE__
                              write(*,*) "Point:",np,"REF:",b,i,j,k,"to",nb,"PERM:", blocks(nbn) % boundary_cond(EAST) % permutation
                              stop 1
                           end if
                        end if
                     else
                        write(*,*) "Error in UNSTR: Permutation diagonal noch nicht implementiert" &
                                  ,__LINE__,__FILE__
                               write(*,*) "Point:",np,"REF:",b,i,j,k,"to",nb,"PERM:" &
                                        , blocks(b) % boundary_cond(EAST) % permutation &
                                        , blocks(b) % boundary_cond(NORTH) % permutation 
                        stop 1
                     end if
                  end if
               end if
               if (i == blocks(b) % nPoints(1)) then
                  nb = blocks(b) % boundary_cond(EAST) % bc_type
                  if (nb > 0) then
                     select case (blocks(b) % boundary_cond(EAST) % permutation) 
                     case (1)
                        write(*,*) "Transfering Point to connected Block:",b,i,j,k, " to ", nb,1,j,k
                        blocks(nb) % refs(1,j,k) = np
                     case (2)
                        write(*,*) "Transfering Point to connected Block:",b,i,j,k, " to ", nb,1+blocks(nb) % nPoints(1) - j,1,k
                        blocks(nb) % refs(1 + blocks(b) % nPoints(1) - j,1,k) = np
                     case default
                        write(*,*) "Error in UNSTR: Permutation in EAST-Richtung noch nicht implementiert" &
                                  ,__LINE__,__FILE__
                        write(*,*) "Point:",np,"REF:",b,i,j,k,"to",nb,"PERM:", blocks(b) % boundary_cond(EAST) % permutation
                        stop 1
                     end select
                  end if
               end if
               if (j == blocks(b) % nPoints(2)) then
                  nb = blocks(b) % boundary_cond(NORTH) % bc_type
                  if (nb > 0) then
                     if (blocks(b) % boundary_cond(NORTH) % permutation == 1) then
                        write(*,*) "Transfering Point to connected Block:",b,i,j,k, " to ", nb,i,1,k
                        blocks(nb) % refs(i,1,k) = np
                     else
                        write(*,*) "Error in UNSTR: Permutation in NORTH-Richtung noch nicht implementiert" &
                                  ,__LINE__,__FILE__
                        write(*,*) "Point:",np,"REF:",b,i,j,k,"to",nb,"PERM:", blocks(b) % boundary_cond(EAST) % permutation
                        stop 1
                     end if
                  end if
               end if
            end if

   !====================================================================================================
   !==========================      CREATE EDGES    ====================================================
   !====================================================================================================
            if (i > 1) then
               edge_exists = .false.
               ! check if the edge already exists
               ! This can be the case if we are at a block boundary and there is an other block connected
               if (    (j == 1                      .and. blocks(b) % boundary_cond(SOUTH) % bc_type > 0)  &       ! at south boundary and a connected block
                  .or. (j == blocks(b) % nPoints(2) .and. blocks(b) % boundary_cond(NORTH) % bc_type > 0)) then    ! at north boundary and a connected block

               write(*,*) "Checking connection from ",i,j,k," to ",i-1,j,k
                  ! A connection from p (i,j,k) to ps (i-1,j,k) shall be created
               
                  p = blocks(b) % refs(i,j,k)
                  ! check all existing edges from point p and if one edge already connects to ps skip creation of this edge
                  ne = git % point_nedges(p)
                  if (ne >0) then
                     ps = blocks(b) % refs(i-1,j,k)
                     do pe = 1,ne
                        e2 = git % point_edges(pe,p)
                        do p2 = 1,2
                           p1 = git % edge_points(p2,e2)
                           if (ps == p1) then
                              edge_exists = .true.
                              write(*,*) "Connection exists: ",e2
                              exit
                           end if
                        end do
                        if (edge_exists) exit
                     end do
                  end if 
               end if
               if (.not. edge_exists) then
                  e = e + 1
                  ! Add points to edge and edge to points
                  do p = 1, 2
                     ! first point i-1, second point i
                     p1 = blocks(b) % refs(i-(2-p),j,k)
                     ! Add point to edge
                     git % edge_points(p,e) = p1
                     ! Add edge to the first point
                     !Increase points edge count
                     pe  = git % point_nedges(p1) + 1
                     git % point_nedges(p1) = pe
                     ! add edge to point at current edge count
                     git % point_edges(pe,p1) = e
                     ! add sign of the resulting edge force, for p2 -> -1 p1 -> 1
                     git % point_edge_signs(pe,p1) =  dble (1-(p-1)*2)
                  end do
   !====================================================================================================
   !==========================   CONNECT NEIGHBOR EDGES ================================================
   !====================================================================================================
                  if (i > 1) then
                     ne = -1
                     p1 = git % edge_points(1,e) ! Left point of Edge
                     do_connect = .true.
                     if (i > 2) then
                        ps = blocks(b) % refs(i - 2,j,k) ! Point_Soll Point we are looking for
                     else if (blocks(b) % boundary_cond(WEST) % bc_type > 0) then
                        ! eventually there is a edge on another blockto connect to
                        nb = blocks(b) % boundary_cond(WEST) % bc_type
                        if (blocks(b) % boundary_cond(WEST) % permutation == 1) then
                           ps = blocks(nb) % refs(blocks(nb) % nCells(1),j,k)
                        else
                           do_connect = .false.
                           write(*,*) "Error in UNSTR: Permutation in i-Richtung noch nicht implementiert" &
                                     ,__LINE__,__FILE__
                           stop 1
                        end if
                     else
                        do_connect = .false.
                     end if
                     if (do_connect) then
                        do eop = 1, git % point_nedges(p1)
                           e2 = git % point_edges(eop,p1) ! eop'th Edge of Point p1
                           p2 = git % edge_points(1,e2)
                           if (p2 == ps) then
                              ne = e2
                              exit
                           end if
                        end do
                        if (ne == -1) then
                           write(*,*) "Error Neighbor edge not found"
                           write(*,*) b,i,j,k
                           write(*,*) p1, git % point_edges(:,p1)
                           stop 1
                        end if

                        write(*,*) b,i,j,k, np,e ,p1,ps
                        write(*,*) "Connecting: ",b,i,j,k," with ",git % point_refs(:,p2)
                        nne = git % edge_nneighbor(e) + 1
                        git % edge_nneighbor(e) = nne
                        git % edge_neighbor(nne,e) = ne
                        nne = git % edge_nneighbor(ne) + 1
                        git % edge_nneighbor(ne) = nne
                        git % edge_neighbor(nne,ne) = e
                     end if
                  end if
               end if
            end if
   !====================================================================================================
   !==========================      CREATE EDGES    ====================================================
   !====================================================================================================
            if (j > 1) then
               edge_exists = .false.
               ! check if the edge already exists
               ! This can be the case if we are at a block boundary and there is an other block connected
               if (    (i == 1                      .and. blocks(b) % boundary_cond(WEST) % bc_type > 0)  &       ! at West boundary and a connected block
                  .or. (i == blocks(b) % nPoints(1) .and. blocks(b) % boundary_cond(EAST) % bc_type > 0)) then    ! at east boundary and a connected block
                  ! A connection from p (i,j,k) to ps (i,j-1,k) shall be created
                  p = blocks(b) % refs(i,j,k)
                  ! check all existing edges from point p and if one edge already connects to ps skip creation of this edge
                  ne = git % point_nedges(p)
                  if (ne > 0) then
                     ps = blocks(b) % refs(i,j-1,k)
                     ! loop over all edges of p
                     do pe = 1,ne
                        e2 = git % point_edges(pe,p)
                        do p2 = 1,2
                           p1 = git % edge_points(p2,e2)
                           if (ps == p1) then
                              edge_exists = .true.
                              exit
                           end if
                        end do
                        if (edge_exists) exit
                     end do
                  end if 
               end if
               if (.not. edge_exists) then
                  e = e + 1
                  ! Add points to edge and edge to points
                  do p = 1, 2
                     ! first point j-1, second point j
                     p1 = blocks(b) % refs(i,j-(2-p),k)
                     ! Add point to edge
                     git % edge_points(p,e) = p1
                     ! Add edge to the first point
                     !Increase points edge count
                     pe  = git % point_nedges(p1) + 1
                     git % point_nedges(p1) = pe
                     ! add edge to point at current edge count
                     git % point_edges(pe,p1) = e
                     ! add sign of the resulting edge force, for p2 -> -1 p1 -> 1
                     git % point_edge_signs(pe,p1) =  dble (1-(p-1)*2)
                  end do
   !====================================================================================================
   !==========================   CONNECT NEIGHBOR EDGES ================================================
   !====================================================================================================
                  if (j > 1) then
                     ne = -1
                     p1 = git % edge_points(1,e) ! Lower point of Edge
                     do_connect = .true.
                     ! if j > 2 there exists another edge in the negative j direction
                     ! Unfortunatelly there is no direct way to get this edge, thus
                     ! we compare all edges of the lower point and see if there first point's
                     ! reference is j-2
                     if (j > 2) then
                        ps = blocks(b) % refs(i,j - 2,k) ! Point_Soll Point we are looking for
                     else if (blocks(b) % boundary_cond(SOUTH) % bc_type > 0) then
                        ! eventually there is a edge on another blockto connect to
                        nb = blocks(b) % boundary_cond(SOUTH) % bc_type
                        if (blocks(b) % boundary_cond(SOUTH) % permutation == 1) then
                           ps = blocks(nb) % refs(i,blocks(nb) % nCells(2),k)
                        else if (blocks(b) % boundary_cond(SOUTH) % permutation == 2) then
                           ps = blocks(nb) % refs(blocks(nb) % nCells(1), 1 + blocks(nb) % nPoints(2) - i,k)
                        else
                           do_connect = .false.
                           write(*,*) "Error in UNSTR: Permutation in i-Richtung noch nicht implementiert" &
                                     ,__LINE__,__FILE__
                           stop 1
                        end if
                     else
                        do_connect = .false.
                     end if
                     if (do_connect) then
                        do eop = 1, git % point_nedges(p1)
                           e2 = git % point_edges(eop,p1) ! eop'th Edge of Point p1
                           p2 = git % edge_points(1,e2)
                           if (p2 == ps) then
                              ne = e2
                              exit
                           end if
                        end do
                        if (ne == -1) then
                           write(*,*) "Error Neighbor edge not found",__LINE__,__FILE__
                           write(*,*) b,i,j,k
                           write(*,*) ps, git % point_edges(:,p1)
                           stop 1
                        end if
                        
                        nne = git % edge_nneighbor(e) + 1 ! Increasing neighbor edge count
                        git % edge_nneighbor(e) = nne     ! Increasing neighbor edge count
                        git % edge_neighbor(nne,e) = ne   ! Referncing new neighbor edge
                        nne = git % edge_nneighbor(ne) + 1! Increasing neighbor edge count of neighbor edge
                        git % edge_nneighbor(ne) = nne
                        git % edge_neighbor(nne,ne) = e   ! REferenceing current edge in neighbor edge's neighbor edge array
                     end if
                  end if
              end if
           end if
         end do
      end do
   end do
end do
!!!!!!!! TEST ARRAY ASSUMPTIONS
if (np /= git % nPoint) then
   write(*,*) "Number of Points wrongly approximated"
   write(*,*) np,git % nPoint
   stop 1
end if
if (e /= git % nEdge) then
   write(*,*) "Number of Edges wrongly approximated"
   write(*,*) e, git % nEdge
   stop 1
end if
write(*,*) "Points"
do np = 1, git % nPoint
   pe = git % point_nedges(np)
   write(*,'("#P ",I4," Ref ",4I4," NEdge ",I2," Edges: ",6I4)' ) &
            np, git % point_refs(:,np),pe, git % point_edges(1:pe,np)
end do
write(*,*) "Edges"
p1 = 0
p2 = 0
do e = 1, git % nedge
   !write(*,'("#E ",I4," Points: ",2I4," Neighbors: ",2I4)') e, git % edge_points(:,e), git % edge_neighbor(:,e)
   if (git % edge_nneighbor(e) == 2) then
      p1 = p1 + 1
   write(*,'("#E ",I4," Points: ",2I4," Neighbors: ",2I4, " Pointrange: ",4I4)') &
         e, git % edge_points(:,e), git % edge_neighbor(:,e) &
         , git % edge_points(1,git % edge_neighbor(1,e)) &
         , git % edge_points(:,e) &
         , git % edge_points(2,git % edge_neighbor(2,e)) 
   else
      p2 = p2 + 1
   write(*,'("#E ",I4," Points: ",2I4," Neighbors: ",I4,4X, " Pointrange: ",4I4)') &
         e, git % edge_points(:,e), git % edge_neighbor(1,e) &
         , git % edge_points(:,e) &
         , git % edge_points(:,git % edge_neighbor(1,e)) 
   end if
end do
write(*,*) p1,p2
git % edge_springs = 1.0e0_REAL_KIND

git % edge_lengths = -1
git % edge_vectors = -1
git % edge_forces = -1

end subroutine strukt2unstr

subroutine calc_edge_length()
implicit none
integer :: e,p1,p2
real(REAL_KIND) :: x1,x2,y1,y2,z1,z2
real(REAL_KIND) :: dx,dy,dz
do e = 1, git % nedge
   p1 = git % edge_points(1,e)
   p2 = git % edge_points(2,e)
   x1 = git % point_coords(1,p1)
   x2 = git % point_coords(1,p2)
   y1 = git % point_coords(2,p1)
   y2 = git % point_coords(2,p2)
   z1 = git % point_coords(3,p1)
   z2 = git % point_coords(3,p2)
   dx = x2-x1
   dy = y2-y1
   dz = z2-z1
   git % edge_vectors(1,e) = dx
   git % edge_vectors(2,e) = dy
   git % edge_vectors(3,e) = dz
   git % edge_lengths(e) = sqrt( dx * dx &
                               + dy * dy &
                               + dz * dz )
   !write(*,*) e, git % edge_lengths(e), git % edge_vectors(:,e)
end do
end subroutine calc_edge_length

subroutine calc_edge_forces(max_f)
implicit none
real(REAL_KIND), intent(out) :: max_f

integer :: e

max_f = 0.0E0_REAL_KIND
do e = 1, git % nedge
   git % edge_forces(:,e) = git % edge_springs(e) * git % edge_vectors(:,e)
   max_f = max(max_f,sqrt( git % edge_forces(1,e) * git % edge_forces(1,e) &
                         + git % edge_forces(2,e) * git % edge_forces(2,e) &
                         + git % edge_forces(3,e) * git % edge_forces(3,e) ))
end do
end subroutine calc_edge_forces

subroutine calc_point_forces(max_f)
implicit none
real(REAL_KIND), intent(out) :: max_f

integer :: np, e, edge
real(REAL_KIND) :: tmp(3), sp
max_f = 0.0E0_REAL_KIND
do np = 1, git % nPoint
    tmp = 0.0E0_REAL_KIND
    !write(*,*) np
   do e = 1, git % point_nedges(np)
      edge = git % point_edges(e,np)
      tmp = tmp + git % edge_forces(:,edge) * git % point_edge_signs(e,np)
      !write(*,*) "->   ",e,edge,git % edge_forces(:,edge) , git % point_edge_signs(e,np)
   end do
   !write(*,*) tmp, git % point_move_rest_vector(:,np)
   if (git % point_move_rest(np)) then

      sp = tmp(1)*git % point_move_rest_vector(1,np) &
         + tmp(2)*git % point_move_rest_vector(2,np) &
         + tmp(3)*git % point_move_rest_vector(3,np)

      tmp = git % point_move_rest_vector(:,np) * sp
!      write(*,*) "restricting:",np,git % point_refs(2:4,np),sp,tmp
   end if
      
   git % point_forces(:,np) = tmp
   max_f = max(max_f,sqrt(tmp(1)*tmp(1)+tmp(2)*tmp(2)+tmp(3)*tmp(3)))
   !write(*,*) np,git % point_nedges(np), git % point_forces(:,np)
end do
end subroutine calc_point_forces
subroutine move_points(max_f)
implicit none
real(REAL_KIND), parameter :: POINT_WEIGTH =  2.5E0_REAL_KIND
real(REAL_KIND), intent(in) :: max_f
real(REAL_KIND) :: min_l, fk

integer :: np

min_l = minval(git % edge_lengths)
fk = min(1.0E0_REAL_KIND,min_l / max_f)
!write(*,*) "FK",min_l,max_f,fk
do np = 1, git % npoint
   git % point_coords(:,np) = git % point_coords(:,np) + git % point_forces(:,np) / POINT_WEIGTH * fk
end do
end subroutine move_points
end module unstr
