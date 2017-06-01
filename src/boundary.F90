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

integer :: p, dn

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
end module boundary
