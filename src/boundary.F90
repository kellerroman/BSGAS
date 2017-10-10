module boundary
use help_routines
use types
implicit none
   character(len=*), parameter :: BC_SHORT(6) = ["w","e","s","n","f","b"]
contains

subroutine read_boundary(blocks)
implicit none
!================================================================================================================================
!===      reads boundary conditions file (bc.cfg)                                                                             ===
!===                                                                                                                          ===
!===      AUTHORS:      Ronan Keller                                                                                          ===
!===      START  :      04.06.2017                                                                                            ===
!===                                                                                                                          ===
!================================================================================================================================
character(len=*), parameter :: filename_bc = "bc.cfg"

type(t_block) :: blocks(:)
integer :: nBlock
integer :: fu, io_stat, pos
logical :: fexists
character(len=100) :: line
character(len=100) :: block_list
character(len=100) :: temp
real(REAL_KIND) :: dn_value

integer :: stat
! stat = 0: Expecting new boundary condition
! stat = 1: Expecting a variable definiton (dn)
integer :: b,f

nBlock = ubound(blocks,1)
inquire(file=trim(filename_bc),exist=fexists)
if (.not. fexists) then
   write(*,*) "Boundary file: '"//trim(filename_bc)//"' not found!",__FILE__,__LINE__
   stop 1
end if
stat = 0
open(newunit = fu, file= trim(filename_bc))
do 
   read(fu,'(A)',iostat=io_stat) line

   if (io_stat < 0) then
      !write(*,*) "End of File"
      exit
   end if
   line = trim(adjustl(line))
   pos = index(line,"!")
   if (pos > 0 ) then
      line = line(1:pos-1)
   end if
   if (len_trim(line) == 0) cycle
   line = lower_case(line)
   if (stat == 0) then
      pos = index(line,":")
      if (pos == 0) then
         write(*,*) "Expecting new Boundary Type"
         stop 1
      end if

      block_list = trim(adjustl(line(pos+1:)))
      stat = 1
   else if (stat == 1) then
      pos = index(line,"=")
      if (pos == 0) then
         write(*,*) "Expecting Value Definition"
         stop 1
      end if
      temp =trim(adjustl(line(1:pos-1)))
      if (temp =="dn") then
         temp = trim(adjustl(line(pos+1:)))
         read(temp,*) dn_value
      else
         write(*,*) "Variable "//temp//" unknown"
         stop 1
      end if
      stat = 0
      fexists = .false.
      !separate the individuel blocks that are separated by a comma
      do 
         pos = index(block_list,",")
         if (pos > 0) then
            temp = trim(adjustl(block_list(1:pos-1)))
            block_list = trim(adjustl(block_list(pos+1:)))
         else
            temp = trim(adjustl(block_list))
            fexists = .true.
         end if
         pos = 1
         do  
            if (ichar(temp(pos:pos)) > 57) then
               read(temp(1:pos-1),*) b
               if (b > nBlock) then
                  write(*,*) "Block from Boundary file does not exists"
                  stop 1
               end if
               f = 0
               do 
                  f = f + 1
                  if (BC_SHORT(f) == trim(temp(pos:pos))) then
                     exit
                  else if (f > 5) then
                     write(*,*) "Face not Recognized:", trim(temp(pos:pos))
                     stop 1
                  end if
               end do
               if (blocks(b) % boundary_cond(f) % bc_type > 0) then
                  write(*,*) "Block",b," Face ",f," is connected to another block and cannot be set as a special boundary"
                  stop 1
               else
                  blocks(b) % boundary_cond(f) % bc_type = -1
                  blocks(b) % boundary_cond(f) % dn = dn_value
               end if
               exit
            else
               pos = pos + 1
            end if
         end do

         if (fexists) exit

      end do
   else
      write(*,*) "Error with Status"
   end if
end do
close(fu)
end subroutine read_boundary

subroutine init_boundary(git, blocks)
use structured_grid, only: number_of_face
implicit none
integer, parameter :: NV_dir(2,4) = reshape([-1,-1,0,-1,0,0,-1,0],[2,4])
!< Discribes the relative position of the 4 neighbor NormalVectors with respect
!to two Vectors. Example n-t NORMAL VECTOR Position:
! fi = i + NV_dir(1,n) * d1(1) + NV_dir(2,n) * d2(1)
! fj = j + NV_dir(1,n) * d1(2) + NV_dir(2,n) * d2(2)
! fk = k + NV_dir(1,n) * d1(3) + NV_dir(2,n) * d2(3)
type(t_block) :: blocks(:)
type(t_unstr) :: git

type :: t_norm
   real(REAL_KIND),allocatable :: nv(:,:,:,:)
end type t_norm

type(t_norm), allocatable :: norms(:)
integer :: nBlock
integer :: b,i,j,k,f

integer :: p,dn

integer :: is,ie,id
integer :: js,je,jd
integer :: ks,ke,kd

integer :: fi,fj,fk,n

integer :: d,vid(3),vec_count,vv2(3),d1(3),d2(3)

real(REAL_KIND) :: v1(3), v2(3),vec(3,2)

logical :: is_3D

nBlock = ubound(blocks,1)
allocate(norms(nBlock))
if (blocks(1) % nPoints(2) == 1) then
   do b = 1, nBlock
      do i = 1,blocks(b) % nPoints(1), blocks(b) % nCells(1)
         p = blocks(b) % refs(i,1,1)
         git % point_move_rest(p) = .true.
         git % point_move_rest_type(p) = 1
         git % point_move_rest_vector(:,p) = 0.0d0
         write(*,*) "Restiricting Movement for",i,p
      end do
   end do
   return
else if (blocks(1) % nPoints(3) > 1) then
   is_3D = .true.
else
   is_3D = .false.
end if
!====================================================================================================
!==========================      CREATE MOVEMENT RESTRICTION INFORMATION  ===========================
!====================================================================================================
git % point_move_rest = .false.
git % point_move_rest_type = 0
k = 1
do b = 1, nBlock
   if (is_3D) then
   
   allocate(norms(b) % nv(3,blocks(b) % nPoints(1), blocks(b) % nPoints(2), blocks(b) % nPoints(3)))
   do f = 1, number_of_face
   associate(bc => blocks(b) % boundary_cond(f))
      if (bc % bc_type <= 0) then
         is = bc % is
         ie = bc % ie
         js = bc % js
         je = bc % je
         ks = bc % ks
         ke = bc % ke
         id = bc % id
         jd = bc % jd
         kd = bc % kd
         ! dimension not normal to surface are 
         vid(1) = 1-abs(id)
         vid(2) = 1-abs(jd)
         vid(3) = 1-abs(kd)
         do k = ks,ke - vid(3)
            do j = js,je - vid(2)
               do i = is,ie - vid(1)
                  vec_count = 0
                  do d = 1,3
                     if (vid(d) == 0) cycle
                     vv2 = 0
                     vv2(d) = 1
                     ! not last cell
                     if (vv2(1) * i + vv2(2) * j + vv2(3) * k < vv2(1) * ie + vv2(2) * je + vv2(3) * ke) then
                        vec_count = vec_count + 1
                        vec(:,vec_count) = blocks(b) % coords(i+vv2(1),j+vv2(2),k+vv2(3),:) &
                                         - blocks(b) % coords(i       ,j       ,k       ,:) 
                     end if
                  end do
                  call cross_product(vec(:,1),vec(:,2),norms(B) % nv(:,i,j,k))
               end do
            end do
         end do
         d1 = vid
         d2 = vid
         ! vid has to dimensions with 1 , d1 and d2 should only face in one
         ! direction each -> d1, last "1" is deleted, d2 first "1" is deleted
         ! example: vid =[1,0,1] -> d1 = [1,0,0]; d2 = [0,0,1]
         do d = 3,1,-1
            if (d1(d) == 1) then
               d1 (d) = 0
               exit
            end if
         end do
         do d = 1,3
            if (d2(d) == 1) then
               d2 (d) = 0
               exit
            end if
         end do
         do k = ks,ke
            do j = js,je
               do i = is,ie
                  p = blocks(b) % refs(i,j,k)
                  ! loop over the 4 neighbor Faces
                  do n = 1,4
                     fi = i + NV_dir(1,n) * d1(1) + NV_dir(2,n) * d2(1)
                     fj = j + NV_dir(1,n) * d1(2) + NV_dir(2,n) * d2(2)
                     fk = k + NV_dir(1,n) * d1(3) + NV_dir(2,n) * d2(3)
                     if (  fi <= ie - vid(1) .and. fi >= is &
                     .and. fj <= je - vid(2) .and. fj >= js &
                     .and. fk <= ke - vid(3) .and. fk >= ks ) then
                        v1 = norms(b) % nv(:,fi,fj,fk)
                        if (git % point_move_rest(p)) then
                           v2 = git % point_move_rest_vector(:,p)
                           select case (git % point_move_rest_type(p) ) 
                           case (3) 
                              if (.NOT. vec_same(v1,v2)) then   
                                 git % point_move_rest_type(p) = 2
                                 call cross_product(v1,v2,git % point_move_rest_vector(:,p))
                              end if
                           case (2)
                              if (abs(scalar_product(v1,v2)) >= 10E-8) then
                                 git % point_move_rest_type(p) = 1
                                 git % point_move_rest_vector(:,p) = [0,0,0]
                              end if
                           case (1)
                           end select
                        else
                           git % point_move_rest(p) = .true.
                           git % point_move_rest_type(p) = 3
                           git % point_move_rest_vector(:,p) = v1
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end if
   end associate
   end do 
   else !is_3D

   if (blocks(b) % boundary_cond(1) % bc_type <= 0) then !!! NO BLOCK CONNECTION WEST SIDE
      i = 1
      do j = 2, blocks(b) % nCells(2)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         git % point_move_rest_type(p) = 2
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
         git % point_move_rest_type(p) = 2
      end do
   end if
   if (blocks(b) % boundary_cond(2) % bc_type <= 0) then !!! NO BLOCK CONNECTION EAST SIDE
      i = blocks(b) % nPoints(1)
      do j = 2, blocks(b) % nCells(2)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         git % point_move_rest_type(p) = 2
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
         git % point_move_rest_type(p) = 2
      end do
   end if
   if (blocks(b) % boundary_cond(3) % bc_type <= 0) then !!! NO BLOCK CONNECTION SOUTH SIDE
      j= 1
      do i = 2, blocks(b) % nCells(1)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         git % point_move_rest_type(p) = 2
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
         git % point_move_rest_type(p) = 2
      end do
   end if
   if (blocks(b) % boundary_cond(4) % bc_type <= 0) then !!! NO BLOCK CONNECTION NORTH SIDE
      j= blocks(b) % nPoints(2)
      do i = 2, blocks(b) % nCells(1)
         p = blocks(b) % refs(i,j,k)
         git % point_move_rest(p) = .true.
         git % point_move_rest_type(p) = 2
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
         git % point_move_rest_type(p) = 2
      end do
   end if
   end if !is_3D
end do
end subroutine init_boundary

subroutine init_walledges(git, blocks)
use structured_grid, only: number_of_face
implicit none
type(t_block) :: blocks(:)
type(t_unstr) :: git
integer :: nBlock
integer :: b,i,j,k

integer :: f                 ! loop var over block boundaries
integer :: is,ie
integer :: js,je
integer :: ks,ke
integer :: di,dj,dk
integer :: p1,p2,pt,ep
integer :: e, eop

integer :: ne                ! Number of edges
integer :: nwe ! NUMBER WALL EDGES

!real(REAL_KIND) :: ref(4)
nBlock = ubound(blocks,1)
!====================================================================================================
!====================    CREATE LIST OF WALL EDGES WITH SPECIFIC REQUIRED LENGTH   ==================
!====================================================================================================
git % nWallEdge = 0
!write(*,*) "Walledges"
do b = 1, nBlock
   do f = 1,number_of_face
      if (blocks(b) % boundary_cond(f) % bc_type < 0) then
         i = blocks(b) % nPoints(1)
         j = blocks(b) % nPoints(2)
         k = blocks(b) % nPoints(3)
         if (f <= EAST) then
            i = 1
         else if (f <= NORTH) then
            j = 1
         else if (f <= BACK) then
            k = 1
         end if
         do ne = 1,number_of_face,2
            if (blocks(b) % boundary_cond(ne) % bc_type > 0) then
               if      (ne <= EAST ) then
                  i = max(i - 1,1)
               else if (ne <= NORTH) then
                  j = max(j - 1,1)
               else if (ne <= BACK ) then
                  k = max(k - 1,1)
               end if
            end if
         end do
         git % nWallEdge = git % nWallEdge + i * j * k
      end if
   end do
end do

call alloc(git % wall_edges         , git % nWallEdge)
call alloc(git % wall_edge_dns      , git % nWallEdge)
nwe = 0
do b = 1, nBlock
   do f = 1,number_of_face
   associate(bc => blocks(b) % boundary_cond(f))
      if (bc % bc_type < 0) then
         is = bc % is
         ie = bc % ie
         js = bc % js
         je = bc % je
         ks = bc % ks
         ke = bc % ke
         di = -bc % id
         dj = -bc % jd
         dk = -bc % kd
         ! ignore edges a boundary if there is a connected block   
         do ne = 1, number_of_face,2
            if (blocks(b) % boundary_cond(ne) % bc_type > 0) then
               if (ne <= EAST) then
                  if (is /= ie) &
                  is = is + 1
               else if (ne <= NORTH) then
                  if (js /= je) &
                  js = js + 1
               else if (ne <= BACK) then
                  if (ks /= ke) &
                  ks = ks + 1
               end if
            end if
         end do
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
                           git % wall_edge_dns(nwe) = blocks(b) % boundary_cond(f) % dn
                        end if
                     end do

                  end do
               end do
            end do
         end do
      end if
   end associate
   end do
end do
if (nwe /= git % nWallEdge) then
   write(*,*) "Number of WallEdges wrongly approximated"
   write(*,*) nwe, git % nWallEdge
   stop 1
end if
!do nwe = 1, git % nWallEdge
!   e = git%wall_edges(nwe)
!   p1 = git % edge_points(1,e)
!   p2 = git % edge_points(2,e)
!   write(*,'("# ",I4," #E: ",I4," p: ",2I4," Ref:",4(1X,F5.1))') nwe,e,git % edge_points(:,e) &
!      ,dble(git % point_refs(:,p1)+git % point_refs(:,p2)) / 2.0E0_REAL_KIND
!end do
end subroutine init_walledges
end module boundary
