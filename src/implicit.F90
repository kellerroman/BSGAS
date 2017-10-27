module implicit
! **************************************************************************************************
! ***                          Module Implicit                                                   ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   20.10.2017
! Last changes: 20.10.2017
! Version:      V0.1.0
! --------------------------------------------------------------------------------------------------
! Description:
!   This module contains routines and variables associated with the implicit
!   solver, initialization, matrix values, and call of the external solver
!   (SuperLU)
!   
!   
! --------------------------------------------------------------------------------------------------
! Comments and Notes:
!   
! --------------------------------------------------------------------------------------------------
! References:
!
! --------------------------------------------------------------------------------------------------
! Author and Change History:
!   - 2017-10-20,RK : Creation of Implicit Module
!
! **************************************************************************************************
use types
use help_routines
use screen_io
implicit none

integer, parameter :: NRHS = 1
character(len=1), parameter :: NUM2DIM(3) = ["X","Y","Z"]

integer(INT_KIND) :: nPoint

integer(INT_KIND), allocatable :: Point_ids(:)
integer(INT_KIND), allocatable :: Point_dims(:)

integer(INT_KIND), allocatable :: GridPkt2MatPkt(:,:)
!< Id of the Point in matrix id

integer(INT_KIND), allocatable :: edge_nmat_entry(:)
integer(INT_KIND), allocatable :: edge_pos_mat_entry(:,:)

integer(INT_KIND), allocatable :: edge_nrhs_entry(:)
integer(INT_KIND), allocatable :: edge_pos_rhs_entry(:,:)
real(REAL_KIND)  , allocatable :: edge_val_rhs_entry(:,:)

integer(INT_KIND) :: nNonZero

real(REAL_KIND)  , allocatable :: values(:)
real(REAL_KIND)  , allocatable :: rhs   (:)

real(REAL_KIND)  , allocatable :: res   (:)
real(REAL_KIND)  , allocatable :: sol   (:)

integer(INT_KIND), allocatable :: colPos(:)
integer(INT_KIND), allocatable :: rowPos(:)

integer(INT_KIND), allocatable :: mat_dia(:)

integer(KIND = 8)             ::  factors
contains

subroutine init_matrix(git)
implicit none
type(t_unstr), intent(in) :: git

integer(INT_KIND) :: p,d,n
integer(INT_KIND) :: eq, pos
integer(INT_KIND) :: gpkt1, gpkt2, d1 , d2
integer(INT_KIND) :: found, mat_pos
!integer :: iopt, info


nNonZero = 0
nPoint   = 0

call alloc( GridPkt2MatPkt,git % dimension,git % nPoint)
call alloc( edge_nmat_entry, git % nEdge)
call alloc( edge_pos_mat_entry,git % dimension * 4, git % nEdge)

call alloc( edge_nrhs_entry, git % nEdge)
call alloc( edge_pos_rhs_entry,git % dimension, git % nEdge)
call alloc( edge_val_rhs_entry,git % dimension, git % nEdge)

GridPkt2MatPkt = -1

do p = 1, git % nPoint
   do d = 1, git % dimension
      if (git % point_move_dim_rest(d,p)) then
      else
         nPoint = nPoint + 1
      end if
   end do
end do

call sw_info_1_int("Points in Matrix:", nPoint)

call alloc(Point_ids    ,nPoint)
call alloc(Point_dims   ,nPoint)
n = 0
do p = 1, git % nPoint
   do d = 1, git % dimension
      if (git % point_move_dim_rest(d,p)) then
      else
         n = n + 1
         Point_ids (n) = p
         Point_dims(n) = d
         GridPkt2MatPkt(d,p) = n
      end if
   end do
end do

call sw_info_1_int("Maximum id-Dist:",git % max_id_delta)

!do p = 1, nPoint
!   write(*,'(I3,A1,1X)', ADVANCE="NO")  Point_ids(p), NUM2DIM(Point_dims(p))
!end do
!write(*,*)
!do eq = 1, nPoint
!   gpkt1 = Point_ids (eq)
!   d1    = Point_dims(eq)
!   do pos = 1, nPoint
!      d2    = Point_dims(pos)
!      found = 0
!      if (d1 == d2) then
!         if (pos == eq) then
!            do n = 1, git % point_nedges(gpkt1)
!               p = git % point_edges(n,gpkt1)
!               found = found + p
!            end do
!         else
!            gpkt2 = Point_ids (pos)
!            do n = 1, git % point_nedges(gpkt1)
!               if (gpkt2 == git % point_neighbors(n,gpkt1)) then
!                  found = git % point_edges(n,gpkt1)
!                  exit
!               end if
!            end do
!         end if
!      end if
!      if (found /= 0) then
!         write(*,'(" ",I3.3,1X)',ADVANCE="NO") found
!      else  
!         write(*,'("    ",1X)',ADVANCE="NO")
!      end if
!   end do
!   do n = 1, git % point_nedges(gpkt1)
!      if (GridPkt2MatPkt(d1,git % point_neighbors(n,gpkt1)) == -1 ) then
!         write(*,'(" ",I3.3,1X)',ADVANCE="NO") git % point_edges(n,gpkt1)
!      end if
!   end do
!   write(*,'("    ",1X)')
!end do

do eq = 1, nPoint
   gpkt1 = Point_ids (eq)
   d1    = Point_dims(eq)
   do pos = max(1,eq-git % max_id_delta), min(nPoint,eq+git % max_id_delta)
      d2    = Point_dims(pos)
      found = 0
      if (d1 == d2) then
         if (pos == eq) then
            do n = 1, git % point_nedges(gpkt1)
               p = git % point_edges(n,gpkt1)
               found = found + p
            end do
         else
            gpkt2 = Point_ids (pos)
            do n = 1, git % point_nedges(gpkt1)
               if (gpkt2 == git % point_neighbors(n,gpkt1)) then
                  found = git % point_edges(n,gpkt1)
                  exit
               end if
            end do
         end if
      end if
      if (found /= 0) then
   !      write(*,'(" ",I3.3,1X)',ADVANCE="NO") found
         nNonZero = nNonZero + 1
      else  
   !      write(*,'("    ",1X)',ADVANCE="NO")
      end if
   end do
   do n = 1, git % point_nedges(gpkt1)
      if (GridPkt2MatPkt(d1,git % point_neighbors(n,gpkt1)) == -1 ) then
   !      write(*,'(" ",I3.3,1X)',ADVANCE="NO") git % point_edges(n,gpkt1)
      end if
   end do
   !write(*,'("    ",1X)')
end do

call sw_info_1_int("Number of Non-Zero Elements:",nNonZero)

call alloc(colPos       ,nPoint+1)
call alloc(rowPos       ,nNonZero)

call alloc(values       ,nNonZero)
call alloc(rhs          ,nPoint  )
call alloc(mat_dia      ,nPoint  )
call alloc(res          ,nPoint  )
call alloc(sol          ,nPoint  )

mat_pos = 0
do pos = 1, nPoint
   colPos(pos) = mat_pos + 1
   gpkt1 = Point_ids (pos)
   d1    = Point_dims(pos)
   do eq = max(1,pos-git % max_id_delta), min(nPoint,pos+git % max_id_delta)
      d2    = Point_dims(eq)
      if (d1 == d2) then
         if (pos == eq) then
            mat_pos = mat_pos + 1
            rowPos(mat_pos) = eq
            mat_dia(eq) = mat_pos
            do n = 1, git % point_nedges(gpkt1)
               p = git % point_edges(n,gpkt1)
               d = edge_nmat_entry(p) + 1
               edge_nmat_entry(p) = d
               edge_pos_mat_entry(d,p) = mat_pos
            end do
         else
            gpkt2 = Point_ids (eq)
            do n = 1, git % point_nedges(gpkt1)
               if (gpkt2 == git % point_neighbors(n,gpkt1)) then
                  mat_pos = mat_pos + 1
                  rowPos(mat_pos) = eq
                  p = git % point_edges(n,gpkt1)
                  d = edge_nmat_entry(p) + 1
                  edge_nmat_entry(p) = d
                  edge_pos_mat_entry(d,p) = mat_pos
                  exit
               end if
            end do
         end if
      end if
   end do
end do

if (mat_pos /= nNonZero) then
   write(*,*) "Error in NonZero Matrix Value calculation",__FILE__,__LINE__
   stop 1
end if
do eq = 1, nPoint
   d1    = Point_dims(eq)
   gpkt1 = Point_ids (eq)
   do n = 1, git % point_nedges(gpkt1)
      if (GridPkt2MatPkt(d1,git % point_neighbors(n,gpkt1)) == -1 ) then
         p = git % point_edges(n,gpkt1)
         d = edge_nrhs_entry(p) + 1
         gpkt2 = git % point_neighbors(n,gpkt1)
         edge_nrhs_entry(p) = d
         edge_pos_rhs_entry(d,p) = eq
         edge_val_rhs_entry(d,p) = - git % point_coords(d1,gpkt2)
      end if
   end do
end do
colPos(nPoint + 1) = mat_pos + 1

do p = 1, nPoint
    sol(p) = git % point_coords(Point_dims(p), Point_ids(p))
end do
!iopt = 1
!call c_fortran_dgssv( iopt, nPoint, nNonZero, NRHS, values, rowPos, colPos, rhs, nPoint, factors, info )
!!
!if (info .eq. 0) then
!!         write (*,*) 'Factorization succeeded'
!else
!   write(*,*) 'INFO from factorization = ', info
!   stop 1
!endif
end subroutine init_matrix


subroutine calc_matrix(git)
implicit none
type(t_unstr), intent(in) :: git

integer(INT_KIND) :: e,n,pos
values = 0
rhs    = 0
do e = 1, git % nEdge
   do n = 1, edge_nmat_entry(e)
      pos = edge_pos_mat_entry(n,e)
      values(pos) = values(pos) + git % edge_springs(e)
   end do
   do n = 1, edge_nrhs_entry(e)
      pos = edge_pos_rhs_entry(n,e)
      rhs(pos) = rhs(pos) + git % edge_springs(e) * edge_val_rhs_entry(n,e)
   end do
end do
do e = 1, nPoint
   pos = mat_dia(e)
   values(pos) = -values(pos)
end do

end subroutine calc_matrix

subroutine solve_system(max_res,sum_res,point_coords)
implicit none

real(REAL_KIND), intent(out) :: max_res, sum_res
real(REAL_KIND), intent(inout) :: point_coords(:,:)
integer :: iopt, info, i
! First, factorize the matrix. The factors are stored in *factors* handle.
info = 0
iopt = 1
max_res = 1E-10
sum_res = 0
call c_fortran_dgssv( iopt, nPoint, nNonZero, NRHS, values, rowPos, colPos, rhs, nPoint, factors, info )
!
if (info .eq. 0) then
!         write (*,*) 'Factorization succeeded'
else
   write(*,*) 'INFO from factorization = ', info
   stop 1
endif
!
! Second, solve the system using the existing factors.
iopt = 2
call c_fortran_dgssv( iopt, nPoint, nNonZero, NRHS, values, rowPos, colPos, rhs, nPoint, factors, info )
!
if (info .eq. 0) then
!         write (*,*) 'Solve succeeded'
else
   write(*,*) 'INFO from triangular solve = ', info
   stop 1
endif
res = abs(sol - rhs)
sol = rhs
max_res = maxval(res)
sum_res = sum(res)

do i = 1, nPoint
   point_coords(Point_dims(i), Point_ids(i)) = sol(i)
end do

! Last, free the storage allocated inside SuperLU
iopt = 3
call c_fortran_dgssv( iopt, nPoint, nNonZero, NRHS, values, rowPos, colPos, rhs, nPoint, factors, info )
!
end subroutine solve_system

end module implicit
