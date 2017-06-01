module boundary
use help_routines
use types
implicit none

contains
subroutine init_boundary(git, blocks)
type(t_block) :: blocks(:)
type(t_unstr) :: git

integer :: nBlock
integer :: b,i,j,k

integer :: p,dn


real(REAL_KIND) :: v1(3), v2(3)

nBlock = ubound(blocks,1)
if (blocks(1) % nPoints(3) > 1) then
   write(*,*) "Boundary: 3D not supported yet",__FILE__,__LINE__
   stop 1
end if
!====================================================================================================
!==========================      CREATE MOVEMENT RESTRICTION INFORMATION  ===========================
!====================================================================================================
git % point_move_rest = .false.
k = 1
do b = 1, nBlock
   !!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!     WEST SIDE !!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!
   if (blocks(b) % boundary_cond(1) % bc_type <= 0) then !!! NO BLOCK CONNECTION WEST SIDE
      i = 1
      do j = 2, blocks(b) % nCells(2)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         v1 = blocks(b) % coords(i,j+1,k,:) - git % point_coords(:,p)
         v2 = git % point_coords(:,p) - blocks(b) % coords(i,j-1,k,:)
         call vec_common(v1,v2)
         git % point_move_rest_vector(:,p) = v1
      end do
      ! CORNER POINTS
      do j = 1, blocks(b) % nPoints(2), blocks(b) % nCells(2)
         if (j == 1) then
            dn = +1
         else 
            dn = -1
         end if
         p = blocks(b) % refs(i,j,k)
         v1 = blocks(b) % coords(i,j+dn,k,:) - git % point_coords(:,p)
         ! If restriction is already active, then see if same dir as current
         if (git % point_move_rest(p)) then
            v2 = git % point_move_rest_vector(:,p)
            call vec_common(v1,v2)
            git % point_move_rest_vector(:,p) = v1
         else
            git % point_move_rest(p) = .true.
            call vec_normalize(v1)
            git % point_move_rest_vector(:,p) = v1
         end if
      end do
   end if
   if (blocks(b) % boundary_cond(2) % bc_type <= 0) then !!! NO BLOCK CONNECTION EAST SIDE
      i = blocks(b) % nPoints(1)
      do j = 2, blocks(b) % nCells(2)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         v1 = blocks(b) % coords(i,j+1,k,:) - git % point_coords(:,p)
         v2 = git % point_coords(:,p) - blocks(b) % coords(i,j-1,k,:)
         call vec_common(v1,v2)
         git % point_move_rest_vector(:,p) = v1
      end do
      ! CORNER POINTS
      do j = 1, blocks(b) % nPoints(2), blocks(b) % nCells(2)
         if (j == 1) then
            dn = +1
         else 
            dn = -1
         end if
         p = blocks(b) % refs(i,j,k)
         v1 = blocks(b) % coords(i,j+dn,k,:) - git % point_coords(:,p)
         ! If restriction is already active, then see if same dir as current
         if (git % point_move_rest(p)) then
            v2 = git % point_move_rest_vector(:,p)
            call vec_common(v1,v2)
            git % point_move_rest_vector(:,p) = v1
         else
            git % point_move_rest(p) = .true.
            call vec_normalize(v1)
            git % point_move_rest_vector(:,p) = v1
         end if
      end do
   end if
   if (blocks(b) % boundary_cond(3) % bc_type <= 0) then !!! NO BLOCK CONNECTION SOUTH SIDE
      j= 1
      do i = 2, blocks(b) % nCells(1)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         v1 = blocks(b) % coords(i+1,j,k,:) - git % point_coords(:,p)
         v2 = git % point_coords(:,p) - blocks(b) % coords(i-1,j,k,:)
         call vec_common(v1,v2)
         git % point_move_rest_vector(:,p) = v1
      end do
      ! CORNER POINTS
      do i = 1, blocks(b) % nPoints(1), blocks(b) % nCells(1)
         if (i == 1) then
            dn = +1
         else 
            dn = -1
         end if
         p = blocks(b) % refs(i,j,k)
         v1 = blocks(b) % coords(i+dn,j,k,:) - git % point_coords(:,p)
         ! If restriction is already active, then see if same dir as current
         if (git % point_move_rest(p)) then
            v2 = git % point_move_rest_vector(:,p)
            call vec_common(v1,v2)
            git % point_move_rest_vector(:,p) = v1
         else
            git % point_move_rest(p) = .true.
            call vec_normalize(v1)
            git % point_move_rest_vector(:,p) = v1
         end if
      end do
   end if
   if (blocks(b) % boundary_cond(4) % bc_type <= 0) then !!! NO BLOCK CONNECTION NORTH SIDE
      j= blocks(b) % nPoints(2)
      do i = 2, blocks(b) % nCells(1)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         v1 = blocks(b) % coords(i+1,j,k,:) - git % point_coords(:,p)
         v2 = git % point_coords(:,p) - blocks(b) % coords(i-1,j,k,:)
         call vec_common(v1,v2)
         git % point_move_rest_vector(:,p) = v1
      end do
      ! CORNER POINTS
      do i = 1, blocks(b) % nPoints(1), blocks(b) % nCells(1)
         if (i == 1) then
            dn = +1
         else 
            dn = -1
         end if
         p = blocks(b) % refs(i,j,k)
         v1 = blocks(b) % coords(i+dn,j,k,:) - git % point_coords(:,p)
         ! If restriction is already active, then see if same dir as current
         if (git % point_move_rest(p)) then
            v2 = git % point_move_rest_vector(:,p)
            call vec_common(v1,v2)
            git % point_move_rest_vector(:,p) = v1
         else
            git % point_move_rest(p) = .true.
            call vec_normalize(v1)
            git % point_move_rest_vector(:,p) = v1
         end if
      end do
   end if
end do
end subroutine init_boundary
subroutine init_walledges(git, blocks)
type(t_block) :: blocks(:)
type(t_unstr) :: git
integer :: nBlock
integer :: b,i,j,k

integer :: bc                 ! loop var over block boundaries
integer :: is,ie
integer :: js,je
integer :: ks,ke
integer :: di,dj,dk
integer :: p1,p2,pt,ep
integer :: e, eop

integer :: ne                 ! Number of edges
integer :: nwe ! NUMBER WALL EDGES

real(REAL_KIND) :: ref(4)
nBlock = ubound(blocks,1)
!====================================================================================================
!====================    CREATE LIST OF WALL EDGES WITH SPECIFIC REQUIRED LENGTH   ==================
!====================================================================================================
git % nWallEdge = 0
do b = 1, nBlock
   do bc = 1,4
      if (blocks(b) % boundary_cond(bc) % bc_type <= 0) then
         i = blocks(b) % nPoints(1)
         j = blocks(b) % nPoints(2)
         k = blocks(b) % nPoints(3)
         if (bc <= 2) then
            i = 1
         else if (bc <= 4) then
            j = 1
         else
            k = 1
         end if
         write(*,*) b,bc,git % nWallEdge, i,j,k
         git % nWallEdge = git % nWallEdge + i * j * k
      end if
   end do
end do

call alloc(git % wall_edges         , git % nWallEdge)
nwe = 0
do b = 1, nBlock
   do bc = 1,4
      if (blocks(b) % boundary_cond(bc) % bc_type <= 0) then
         is = 1
         ie = blocks(b) % nPoints(1)
         js = 1
         je = blocks(b) % nPoints(2)
         ks = 1
         ke = blocks(b) % nPoints(3)
         di = 0
         dj = 0
         dk = 0
         if (bc == 1) then
            ie = is
            di = 1
         else if (bc == 2) then
            is = ie
            di = -1
         else if (bc == 3) then
            je = js
            dj = 1
         else if (bc == 4) then
            js = je
            dj = -1
         else if (bc == 5) then
            ke = ks
            dk = 1
         else
            ks = ke
            dk = -1
         end if
            
         do k = ks,ke
            do j = js,je
               do i = is,ie
                  p1 = blocks(b) % refs(i,j,k)
                  p2 = blocks(b) % refs(i+di,j+dj,k+dk)
                  ! WALL EDGE
                  ne = git % point_nedges(p1)
                  do eop = 1, ne
                     e = git % point_edges(eop,p1)
                     do ep = 1, 2
                        pt = git % edge_points(ep,e)

                        if (p2 == pt) then
                           nwe = nwe + 1
                           git % wall_edges(nwe)  = e
                        end if
                     end do

                  end do
               end do
            end do
         end do
      end if
   end do
end do
if (nwe /= git % nWallEdge) then
   write(*,*) "Number of WallEdges wrongly approximated"
   write(*,*) nwe, git % nWallEdge
   !stop 1
end if
write(*,*) "Walledges"
do nwe = 1, git % nWallEdge
   e = git%wall_edges(nwe)
   p1 = git % edge_points(1,e)
   p2 = git % edge_points(2,e)
   ref = dble(git % point_refs(:,p1)+git % point_refs(:,p2)) / 2.0E0_REAL_KIND
   write(*,'("# ",I4," #E: ",I4," p: ",2I4," Ref:",4(1X,F5.1))') nwe,e,git % edge_points(:,e),ref
end do
end subroutine init_walledges
end module boundary
