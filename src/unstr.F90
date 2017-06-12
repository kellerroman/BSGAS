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
integer :: i,j,k,np,e,b,nb

integer :: is,ie
integer :: js,je
integer :: ks,ke
integer :: di,dj,dk

integer :: ci,dii,dij,dik
integer :: cj,dji,djj,djk
integer :: ck,dki,dkj,dkk

integer :: p, p1, pe, eop, e2, p2, ps, f

integer :: nne, ne

logical :: edge_exists,do_connect
                  
real(REAL_KIND) :: temp

type(t_same) :: sp

nBlock = ubound(blocks,1)

!write(*,*) "Block structured number of Blocks:",nBlock

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
   !write(*,*) blocks(b) % nPoints,i,j,k

   if (git % dimension == 1) then
      git % nEdge = git % npoint - 1
   else if (git % dimension == 2) then
      git % nEdge = git % nEdge &
                  + blocks(b) % nCells(1) * j &
                  + blocks(b) % nCells(2) * i
   end if


end do

p = 0
temp = 0.0d0
do b = 1, nBlock
   p = p + product(blocks(b) % nPoints)
   do k = 1, blocks(b) % nPoints(3)
      do j = 1, blocks(b) % nPoints(2)
         do i = 1, blocks(b) % nPoints(1)
            nb = blocks(b) % nSamePoints(i,j,k)
            temp = temp + dble(nb) / dble(nb+1)
         end do
      end do
   end do
end do
git % npoint = p-int(temp)
write(*,*) "Points in structured Grid:             ",p
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
call alloc(git % edge_nparallel           , git % nedge)
call alloc(git % edge_parallel         , 2, git % nedge)


git % point_forces = 0.0e0_REAL_KIND
git % point_nedges = 0
git % point_forces = 0.0e0_REAL_KIND

git % edge_nneighbor = 0
git % edge_neighbor  = -1
git % edge_nparallel = 0
git % edge_parallel  = -1
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
               np = np + 1
               git % point_coords(:,np) = blocks(b) % coords(i,j,k,:)
               git % point_refs(1,np) = b
               git % point_refs(2,np) = i
               git % point_refs(3,np) = j
               git % point_refs(4,np) = k
               blocks(b) % refs(i,j,k) = np

   !====================================================================================================
   !==========================   TRANSFER POINT TO OTHER SAME POINTS ===================================
   !====================================================================================================
               do p = 1, blocks(b) % nSamePoints(i,j,k)
                  sp = blocks(b) % SamePoints(p,i,j,k)
                  blocks(sp % b) % refs(sp % i,sp % j, sp % k) = np
               end do
            end if

   !====================================================================================================
   !====================================================================================================
   !==========================      CREATE EDGES    ====================================================
   !==========================      I-DIRECTION     ====================================================
   !====================================================================================================
   !====================================================================================================
            if (i > 1) then
               edge_exists = .false.
               ! check if the edge already exists
               ! This can be the case if we are at a block boundary and there is an other block connected
               if (    (j == 1                      .and. blocks(b) % boundary_cond(SOUTH) % bc_type > 0)  &       ! at south boundary and a connected block
                  .or. (j == blocks(b) % nPoints(2) .and. blocks(b) % boundary_cond(NORTH) % bc_type > 0)) then    ! at north boundary and a connected block
               !write(*,*) "Checking connection from ",i,j,k," to ",i-1,j,k
                  ! A connection from p (i,j,k) to ps (i-1,j,k) shall be created
                  p = blocks(b) % refs(i,j,k)
                  ! check all existing edges from point p and if one edge already connects to ps skip creation of this edge
                  ne = git % point_nedges(p)
                  if (ne > 0) then
                     ps = blocks(b) % refs(i-1,j,k)
                     do pe = 1,ne
                        e2 = git % point_edges(pe,p)
                        do p2 = 1,2
                           p1 = git % edge_points(p2,e2)
                           if (ps == p1) then
                              edge_exists = .true.
                              !write(*,*) "Connection exists: ",e2
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
                  ne = -1
                  p1 = git % edge_points(1,e) ! Left point of Edge
                  do_connect = .true.
                  if (i > 2) then
                     ps = blocks(b) % refs(i - 2,j,k) ! Point_Soll Point we are looking for
                  else if (blocks(b) % boundary_cond(WEST) % bc_type > 0) then
                  associate (bc => blocks(b) % boundary_cond(WEST))
                     ! eventually there is a edge on another blockto connect to
                     nb = bc% bc_type
                     ps = blocks(nb) % refs( &
                           bc % ci + bc % dii * (i-2) + bc % dij * j + bc % dik * k &
                          ,bc % cj + bc % dji * (i-2) + bc % djj * j + bc % djk * k &
                          ,bc % ck + bc % dki * (i-2) + bc % dkj * j + bc % dkk * k )
                  end associate
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

                     !write(*,*) b,i,j,k, np,e ,p1,ps
                     !write(*,*) "Connecting: ",b,i,j,k," with ",git % point_refs(:,p2)
                     nne = git % edge_nneighbor(e) + 1
                     git % edge_nneighbor(e) = nne
                     git % edge_neighbor(nne,e) = ne
                     nne = git % edge_nneighbor(ne) + 1
                     git % edge_nneighbor(ne) = nne
                     git % edge_neighbor(nne,ne) = e
                  end if
   !====================================================================================================
   !==========================   CONNECT PARALLEL EDGES ================================================
   !====================================================================================================
                  ne = -1
                  do_connect = .true.
                  if (j > 1) then
                     p1 = blocks(b) % refs(i-1,j-1,k)
                     ps = blocks(b) % refs(i,j-1,k) ! Point_Soll Point we are looking for
                  else if (blocks(b) % boundary_cond(SOUTH) % bc_type > 0) then
                     ! eventually there is a edge on another blockto connect to
                  associate (bc => blocks(b) % boundary_cond(SOUTH))
                     ! eventually there is a edge on another blockto connect to
                     nb = bc % bc_type
                     p1 = blocks(nb) % refs(i-1,blocks(nb) % nCells(2),k)
                     p2 = blocks(nb) % refs( &
                           bc % ci + bc % dii * (i-1) + bc % dij * (j-1) + bc % dik * k &
                          ,bc % cj + bc % dji * (i-1) + bc % djj * (j-1) + bc % djk * k &
                          ,bc % ck + bc % dki * (i-1) + bc % dkj * (j-1) + bc % dkk * k )
                     write(*,*) p1,p2
                     ps = blocks(nb) % refs(i  ,blocks(nb) % nCells(2),k)
                     p2 = blocks(nb) % refs( &
                           bc % ci + bc % dii * (i  ) + bc % dij * (j-1) + bc % dik * k &
                          ,bc % cj + bc % dji * (i  ) + bc % djj * (j-1) + bc % djk * k &
                          ,bc % ck + bc % dki * (i  ) + bc % dkj * (j-1) + bc % dkk * k )
                     write(*,*) ps,p2
                     stop
                  end associate
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
                        p2 = git % edge_points(2,e2)
                        if (p2 == ps) then
                           ne = e2
                           exit
                        end if
                     end do
                     if (ne == -1) then
                        write(*,*) "Error Parallel Neighbor edge not found"
                        write(*,*) b,i,j,k
                        write(*,*) git % point_refs(:,p1)
                        write(*,*) git % point_refs(:,ps)
                        write(*,*) p1, git % point_edges(:,p1)
                        stop 1
                     end if

                     !write(*,*) b,i,j,k, np,e ,p1,ps
                     !write(*,*) "Connecting: ",b,i,j,k," with ",git % point_refs(:,p2)
                     nne = git % edge_nparallel(e) + 1
                     git % edge_nparallel(e) = nne
                     git % edge_parallel(nne,e) = ne

                     nne = git % edge_nparallel(ne) + 1
                     git % edge_nparallel(ne) = nne
                     git % edge_parallel(nne,ne) = e
                  end if
               end if
            end if
   !====================================================================================================
   !====================================================================================================
   !==========================      CREATE EDGES    ====================================================
   !==========================      J-DIRECTION     ====================================================
   !====================================================================================================
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
   !====================================================================================================
   !==========================   CONNECT PARALLEL EDGES ================================================
   !====================================================================================================
                  ne = -1
                  do_connect = .true.
                  if (i > 1) then
                     p1 = blocks(b) % refs(i-1,j-1,k)
                     ps = blocks(b) % refs(i-1,j,k) ! Point_Soll Point we are looking for
                  else if (blocks(b) % boundary_cond(WEST) % bc_type > 0) then
                     ! eventually there is a edge on another blockto connect to
                  associate (bc => blocks(b) % boundary_cond(WEST))
                     nb = bc % bc_type
                     p1 = blocks(nb) % refs(i-1,blocks(nb) % nCells(2),k)
                     p2 = blocks(nb) % refs( &
                           bc % ci + bc % dii * (i-1) + bc % dij * (j-1) + bc % dik * k &
                          ,bc % cj + bc % dji * (i-1) + bc % djj * (j-1) + bc % djk * k &
                          ,bc % ck + bc % dki * (i-1) + bc % dkj * (j-1) + bc % dkk * k )
                     write(*,*) p1,p2
                     ps = blocks(nb) % refs(i  ,blocks(nb) % nCells(2),k)
                     p2 = blocks(nb) % refs( &
                           bc % ci + bc % dii * (i-1) + bc % dij * (j  ) + bc % dik * k &
                          ,bc % cj + bc % dji * (i-1) + bc % djj * (j  ) + bc % djk * k &
                          ,bc % ck + bc % dki * (i-1) + bc % dkj * (j  ) + bc % dkk * k )
                     write(*,*) ps,p2
                     stop
                  end associate
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
                        p2 = git % edge_points(2,e2)
                        if (p2 == ps) then
                           ne = e2
                           exit
                        end if
                     end do
                     if (ne == -1) then
                        write(*,*) "Error Parallel Neighbor edge not found"
                        write(*,*) b,i,j,k
                        write(*,*) git % point_refs(:,p1)
                        write(*,*) git % point_refs(:,ps)
                        write(*,*) p1, git % point_edges(:,p1)
                        stop 1
                     end if

                     !write(*,*) b,i,j,k, np,e ,p1,ps
                     !write(*,*) "Connecting: ",b,i,j,k," with ",git % point_refs(:,p2)
                     nne = git % edge_nparallel(e) + 1
                     git % edge_nparallel(e) = nne
                     git % edge_parallel(nne,e) = ne

                     nne = git % edge_nparallel(ne) + 1
                     git % edge_nparallel(ne) = nne
                     git % edge_parallel(nne,ne) = e
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


! special cases of Edge connection
do b = 1, nBlock
   do f = 1, 4
   associate(bc => blocks(b) % boundary_cond(f))
      if (bc % bc_type > 0) then
         nb = bc % bc_type
         is = bc % is
         ie = bc % ie
         js = bc % js
         je = bc % je
         ks = bc % ks
         ke = bc % ke
         di = bc % id
         dj = bc % jd
         dk = bc % kd
         do k = ks,ke
            do j = js,je
               do i = is,ie
                  p1 = blocks(b) % refs(i,j,k)
                  p2 = blocks(b) % refs(i-di,j-dj,k-dk)
                  e = -1
                  do ne = 1, git % point_nedges(p1)
                     e2 = git % point_edges(ne,p1)
                     do np = 1,2
                        if (p2 == git % edge_points(np,e2)) then
                           e = e2
                           exit
                        end if
                     end do
                     if (e /= -1) exit
                  end do
                  if (e == -1) then
                     write(*,*) "Edge at boundary not found",b,f,i,j,k,di,dj,dk
                     stop 1
                  end if
                  if (git % edge_nneighbor(e) /= 2) then
                     ci  = bc % ci; dii = bc % dii; dij = bc % dij; dik = bc % dik
                     cj  = bc % cj; dji = bc % dji; djj = bc % djj; djk = bc % djk
                     ck  = bc % ck; dki = bc % dki; dkj = bc % dkj; dkk = bc % dkk
                     
                     p2 = blocks(nb) % refs( &
                                 ci + dii * (i+di) + dij * (j+dj) + dik * (k+dk) &
                                ,cj + dji * (i+di) + djj * (j+dj) + djk * (k+dk) &
                                ,ck + dki * (i+di) + dkj * (j+dj) + dkk * (k+dk) )
                     e2 = -1
                     do ne = 1, git % point_nedges(p1)
                        nne = git % point_edges(ne,p1)
                        do np = 1,2
                           if (p2 == git % edge_points(np,nne)) then
                              e2 = nne
                              exit
                           end if
                        end do
                        if (e2 /= -1) exit
                     end do
                     if (e2 == -1) then
                        write(*,*) "Edge at boundary not found",b,f,i,j,k,p1,p2
                        stop 1
                     end if
                     !write(*,'("Edge ",I0,"@ ",5(I0,1X),"should have two neighbors")') e,b,f,i,j,k,e2
                     nne = git % edge_nneighbor(e) + 1
                     git % edge_neighbor(nne,e) = e2
                     git % edge_nneighbor(e) = nne
                     nne = git % edge_nneighbor(e2) + 1
                     git % edge_neighbor(nne,e2) = e
                     git % edge_nneighbor(e2) = nne
                  end if
               end do
            end do
         end do

      end if
   end associate
   end do
end do


!write(*,*) "Points"
!do np = 1, git % nPoint
!   pe = git % point_nedges(np)
!   write(*,'("#P ",I8," Ref ",4I8," NEdge ",I8," Edges: ",6I8)' ) &
!            np, git % point_refs(:,np),pe, git % point_edges(1:pe,np)
!end do
!write(*,*) "Edges"
!do e = 1, git % nedge
!   !write(*,'("#E ",I4," Points: ",2I4," Neighbors: ",2I4)') e, git % edge_points(:,e), git % edge_neighbor(:,e)
!   if (git % edge_nneighbor(e) == 2) then
!   write(*,'("#E ",I8," Points: ",2I8," Neighbors: ",2I8, " Parallel: ",4I8)') &
!         e, git % edge_points(:,e), git % edge_neighbor(:,e), git % edge_parallel(:,e)
!   else
!   write(*,'("#E ",I8," Points: ",2I8," Neighbors: ",I8,8X, " Parallel: ",4I8)') &
!         e, git % edge_points(:,e), git % edge_neighbor(1,e), git % edge_parallel(:,e)
!   end if
!end do

git % edge_springs = 1.0e0_REAL_KIND

git % edge_lengths = -1
git % edge_vectors = -1
git % edge_forces  = -1

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
