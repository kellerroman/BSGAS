module structured_grid
use hdf5
use const
use types
use help_routines, only: alloc
implicit none

real(REAL_KIND), parameter  :: EPS = 1.0E-6_REAL_KIND

enum, bind(C)
   enumerator :: FILETYPE_UNKNOWM, FILETYPE_HDF5, FILETYPE_UFO
end enum

character(len=100)          :: filename_grid_in
character(len=100)          :: filename_grid_out
character(len=*), parameter :: GROUP_GRID            = "grid"
character(len=*), parameter :: GROUP_BLOCK           = "block"
character(len=*), parameter :: COORD_NAME(3)         = ["CoordinateX","CoordinateY","CoordinateZ"]
integer(INT_KIND)           :: filetype_grid_in      = FILETYPE_UNKNOWM
integer(INT_KIND)           :: filetype_grid_out     = FILETYPE_UNKNOWM

integer(INT_KIND)           :: nBlock
integer(INT_KIND)           :: nCell
integer(INT_KIND)           :: dimen
integer(INT_KIND)           :: solution_output

integer(INT_KIND)           :: number_of_corner_point
integer(INT_KIND)           :: number_of_face
integer(INT_KIND)           :: number_of_points_per_face
integer(INT_KIND)           :: number_of_permutation

type(t_block), allocatable  :: blocks(:)
contains

subroutine read_grid()
implicit none
integer :: pos
logical :: fexists

inquire(file=trim(filename_grid_in),exist=fexists)

if (.not. fexists) then
   write(*,*) "Grid file: '"//trim(filename_grid_in)//"' not found!",__FILE__,__LINE__
   stop 1
end if

if (filetype_grid_in == FILETYPE_UNKNOWM) then
   pos = index(filename_grid_in,".",.TRUE.)
   if (pos > 0) then
      if (trim(filename_grid_in(pos+1:)) == "h5") then
         filetype_grid_in = FILETYPE_HDF5
      else if (trim(filename_grid_in(pos+1:)) == "bin") then
         filetype_grid_in = FILETYPE_UFO
      else
         write(*,*) "Filetype unkown"
         stop 1
      end if
   else
      write(*,*) "Filetype unkown"
      stop 1
   end if
end if

if (filetype_grid_in == FILETYPE_HDF5) then
   call read_grid_hdf5()
else if (filetype_grid_in == FILETYPE_UFO) then
   call read_grid_ufo()
end if
call connect_blocks()
end subroutine read_grid

subroutine write_grid(iter)
implicit none
integer, intent(in), optional :: iter
integer :: pos

if (present(iter)) then
   if (mod(iter,solution_output) /= 0) return
end if
if (filetype_grid_out == FILETYPE_UNKNOWM) then
   filetype_grid_out = FILETYPE_HDF5 ! Default is always HDF5
   pos = index(filename_grid_out,".",.TRUE.)
   if (pos > 0) then
      if (trim(filename_grid_out(pos+1:)) == "h5") then
         filetype_grid_out = FILETYPE_HDF5
      else if (trim(filename_grid_out(pos+1:)) == "bin") then
         filetype_grid_out = FILETYPE_UFO
      end if
   end if
end if


if (filetype_grid_out == FILETYPE_HDF5) then
   call write_grid_hdf5()
else if (filetype_grid_out == FILETYPE_UFO) then
   call write_grid_ufo()
end if
end subroutine write_grid

subroutine read_grid_hdf5()
implicit none

character(len=len(GROUP_BLOCK)+2) :: block_group
real(REAL_KIND), allocatable :: data_in(:,:,:)
integer     ::   error ! Error flag
integer(hid_t) :: file_id       ! file identifier
integer(hid_t) :: group_id_grid      ! dataset identifier
integer(hid_t) :: group_id_block      ! dataset identifier
integer(hid_t) :: dset_id       ! dataset identifier
integer(hid_t) :: dspace_id     ! dataspace identifier
integer(HSIZE_T) :: dims(3)
integer(HSIZE_T) :: maxdims(3)
integer :: b,d
! Initialize FORTRAN interface.
nCell = 0
call h5open_f(error)

! Open an existing file.
call h5fopen_f (filename_grid_in, h5f_acc_rdwr_f, file_id, error)

call h5gopen_f(file_id,GROUP_GRID,group_id_grid,error)

call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)
allocate(blocks(nBlock))
do b = 1, nBlock
   write(block_group,'(A,I0)') GROUP_BLOCK, b
   call h5gopen_f(group_id_grid,block_group,group_id_block,error)
   
   ! Open an existing dataset.
   call h5dopen_f(group_id_block, COORD_NAME(1), dset_id, error)
   call h5dget_space_f(dset_id,dspace_id,error)
   
   call h5sget_simple_extent_ndims_f(dspace_id,dimen,error)
   call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)

   call h5dclose_f(dset_id, error)
   blocks(b) % nPoints = INT(dims,INT_KIND)

   blocks(b) % nCells = max(1,blocks(b) % nPoints - 1)
   nCell = nCell + product(blocks(b) % nCells)
   
   !!!! MUSS VERSCHOBEN WERDEN FUER MULTIBLOCK
!   call allocate_vars(b)

   call alloc (data_in, blocks(b)%nPoints)
   allocate ( blocks(b) % coords (blocks(b)%nPoints(1),blocks(b)%nPoints(2),blocks(b)%nPoints(3),3))
   call alloc (blocks(b) % refs,blocks(b) % nPoints)
   call alloc (blocks(b) % nSamePoints,blocks(b) % nPoints)
   allocate ( blocks(b) % samePoints(8,blocks(b) % nPoints(1) &
                                      ,blocks(b) % nPoints(2) &
                                      ,blocks(b) % nPoints(3) ))
   blocks(b) % refs = -1
   blocks(b) % nSamePoints = 0
   blocks(b) % samePoints(:,:,:,:) % b = -1
   blocks(b) % samePoints(:,:,:,:) % i = -1
   blocks(b) % samePoints(:,:,:,:) % j = -1
   blocks(b) % samePoints(:,:,:,:) % k = -1

   do d = 1,dimen
      call h5dopen_f(group_id_block, COORD_NAME(d), dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_in, dims, error)
      blocks(b) % coords (1:blocks(b)%nPoints(1) &
                         ,1:blocks(b)%nPoints(2) &
                         ,1:blocks(b)%nPoints(3),d) = data_in
      call h5dclose_f(dset_id, error)
   end do

   deallocate(data_in)
   call h5gclose_f(group_id_block, error)
end do
call h5gclose_f(group_id_grid, error) ! CLOSE GRID GROUP

call h5fclose_f(file_id, error)
end subroutine read_grid_hdf5

subroutine read_grid_ufo()
implicit none
integer :: fu
integer :: b,i,j,k
nCell = 0

open(newunit=fu,file=trim(filename_grid_in),form="unformatted",access="stream")

read(fu) dimen, nBlock
if (dimen == 2) then
   dimen = 3
else
   dimen = 2
end if
allocate(blocks(nBlock))
do b = 1, nBlock
   blocks(b) % nPoints = 1
   read(fu) blocks(b) % nPoints(1:dimen)

   blocks(b) % nCells = max(1,blocks(b) % nPoints - 1)
   nCell = nCell + product(blocks(b) % nCells)
end do
do b = 1, nBlock
   allocate ( blocks(b) % coords (blocks(b)%nPoints(1),blocks(b)%nPoints(2),blocks(b)%nPoints(3),3))
   call alloc (blocks(b) % refs,blocks(b) % nPoints)
   call alloc (blocks(b) % nSamePoints,blocks(b) % nPoints)
   allocate ( blocks(b) % samePoints(8,blocks(b) % nPoints(1) &
                                      ,blocks(b) % nPoints(2) &
                                      ,blocks(b) % nPoints(3) ))
   blocks(b) % refs = -1
   blocks(b) % nSamePoints = 0
   blocks(b) % samePoints(:,:,:,:) % b = -1
   blocks(b) % samePoints(:,:,:,:) % i = -1
   blocks(b) % samePoints(:,:,:,:) % j = -1
   blocks(b) % samePoints(:,:,:,:) % k = -1

   read(fu) (((blocks(b) % coords (i,j,k,1:dimen),i = 1, blocks(b) % nPoints(1)) &
                                                 ,j = 1, blocks(b) % nPoints(2)) & 
                                                 ,k = 1, blocks(b) % nPoints(3))  

end do
close(fu)
end subroutine read_grid_ufo

subroutine connect_blocks()
implicit none
integer :: b,f,nb,nf,p,per,d
integer :: ni,nj,nk

logical :: is3D 

integer :: is,ie,id
integer :: js,je,jd
integer :: ks,ke,kd

integer :: ci,dii,dij,dik
integer :: cj,dji,djj,djk
integer :: ck,dki,dkj,dkk

integer :: i,j,k

logical :: found


real(REAL_KIND), allocatable :: corner_points(:,:,:)
!< Array of the Cornerpoints for each block (dim,point_id,block)
integer, allocatable :: point_on_face(:,:)
!< Array with the points for each face (point_id,face)
!< face: (W,E,S,N,F,B)
integer, allocatable :: permutation(:,:)
!< Different Point arrangment due to misaligned blocks (point_id,permutation)

is3D = .false.

b = 1
if (blocks(b) % nPoints(2) == 1) then
   if (nBlock > 1) then
      write(*,*) "For 1D case only single block supported"
      stop 1
   end if
   number_of_corner_point = 2
   number_of_face = 2
   number_of_points_per_face = 2
   number_of_permutation = 1
   allocate( blocks(b) % boundary_cond(6))
   blocks(b) % boundary_cond (:) % bc_type         = 0
   blocks(b) % boundary_cond (:) % face            = 0
   blocks(b) % boundary_cond (:) % permutation     = 0

   blocks(b) % boundary_cond(1) % is   = 1
   blocks(b) % boundary_cond(1) % ie   = 1
   blocks(b) % boundary_cond(1) % id   = -1 

   blocks(b) % boundary_cond(2) % is   = blocks(b) % nPoints(1)
   blocks(b) % boundary_cond(2) % ie   = blocks(b) % nPoints(1)
   blocks(b) % boundary_cond(2) % id   = 1 

   blocks(b) % boundary_cond(:) % js   = 1
   blocks(b) % boundary_cond(:) % je   = 1
   blocks(b) % boundary_cond(:) % jd   = 0

   blocks(b) % boundary_cond(:) % ks   = 1
   blocks(b) % boundary_cond(:) % ke   = 1
   blocks(b) % boundary_cond(:) % kd   = 0
   return
else if  ( blocks(b) % nPoints(3) == 1) then
   number_of_corner_point = 4
   number_of_face = 4
   number_of_points_per_face = 2
   number_of_permutation = 2
else
   is3D = .true.
   number_of_corner_point = 8
   number_of_face = 6
   number_of_points_per_face = 4
   number_of_permutation = 8
end if

call alloc (corner_points,3,number_of_corner_point,nBlock)
call alloc (point_on_face,number_of_points_per_face,number_of_face)
call alloc (permutation,number_of_points_per_face,number_of_permutation)

do b = 1, nBlock
   ni = blocks(b) % nPoints(1)
   nj = blocks(b) % nPoints(2)
   nk = blocks(b) % nPoints(3)

   allocate( blocks(b) % boundary_cond(6))

   blocks(b) % boundary_cond (:) % bc_type         = 0
   blocks(b) % boundary_cond (:) % face   = 0
   blocks(b) % boundary_cond (:) % permutation     = 0

   corner_points (:,1,b) = blocks(b) % coords ( 1, 1, 1,:)
   corner_points (:,2,b) = blocks(b) % coords (ni, 1, 1,:)
   corner_points (:,3,b) = blocks(b) % coords ( 1,nj, 1,:)
   corner_points (:,4,b) = blocks(b) % coords (ni,nj, 1,:)
   if (is3D) then
      corner_points (:,5,b) = blocks(b) % coords ( 1, 1,nk,:)
      corner_points (:,6,b) = blocks(b) % coords (ni, 1,nk,:)
      corner_points (:,7,b) = blocks(b) % coords ( 1,nj,nk,:)
      corner_points (:,8,b) = blocks(b) % coords (ni,nj,nk,:)
   end if
end do


if (is3D) then
   point_on_face(:,1) = [1,3,5,7] !WEST
   point_on_face(:,2) = [2,4,6,8] !EAST
   point_on_face(:,3) = [1,2,5,6] !SOUTH
   point_on_face(:,4) = [3,4,7,8] !NORTH
   point_on_face(:,5) = [1,2,3,4] !FRONT
   point_on_face(:,6) = [5,6,7,8] !BACK

   permutation (:,1)  = [1,2,3,4] ! No twist
   permutation (:,2)  = [2,3,4,1] ! twist by 90°
   permutation (:,3)  = [3,4,1,2]
   permutation (:,4)  = [4,1,2,3]
   permutation (:,5)  = [4,3,2,1]
   permutation (:,6)  = [3,2,1,4]
   permutation (:,7)  = [2,1,4,3]
   permutation (:,8)  = [1,4,3,2]
else
   point_on_face(:,1) = [1,3] !WEST
   point_on_face(:,2) = [2,4] !EAST
   point_on_face(:,3) = [1,2] !SOUTH
   point_on_face(:,4) = [3,4] !NORTH

   permutation (:,1)  = [1,2] ! no twist
   permutation (:,2)  = [2,1] ! Upside down 
end if

do b = 1, nBlock
   do f = 1, number_of_face
   select case(f)
      case(WEST)
         is = 1; ie = is                    ; id = -1
         js = 1; je = blocks(b) % nPoints(2); jd = 0
         ks = 1; ke = blocks(b) % nPoints(3); kd = 0
      case (EAST)
         is = blocks(b) % nPoints(1); ie = is                    ; id = 1
         js = 1                     ; je = blocks(b) % nPoints(2); jd = 0
         ks = 1                     ; ke = blocks(b) % nPoints(3); kd = 0
      case(SOUTH)
         is = 1; ie = blocks(b) % nPoints(1); id = 0
         js = 1; je = js                    ; jd = -1
         ks = 1; ke = blocks(b) % nPoints(3); kd = 0
      case(NORTH)
         is = 1                     ; ie = blocks(b) % nPoints(1); id = 0
         js = blocks(b) % nPoints(2); je = js                    ; jd = 1
         ks = 1                     ; ke = blocks(b) % nPoints(3); kd = 0
      case(FRONT)
         is = 1                     ; ie = blocks(b) % nPoints(1); id = 0
         js = 1                     ; je = blocks(b) % nPoints(2); jd = 0
         ks = 1                     ; ke = ks                    ; kd = -1
      case(BACK)
         is = 1                     ; ie = blocks(b) % nPoints(1); id = 0
         js = 1                     ; je = blocks(b) % nPoints(2); jd = 0
         ks = blocks(b) % nPoints(3); ke = ks                    ; kd = 1
      case default
         write(*,*) "ERROR face not supported yet"
         stop 1
      end select
      blocks(b) % boundary_cond(f) % is   = is
      blocks(b) % boundary_cond(f) % ie   = ie
      blocks(b) % boundary_cond(f) % id   = id

      blocks(b) % boundary_cond(f) % js   = js
      blocks(b) % boundary_cond(f) % je   = je
      blocks(b) % boundary_cond(f) % jd   = jd

      blocks(b) % boundary_cond(f) % ks   = ks
      blocks(b) % boundary_cond(f) % ke   = ke
      blocks(b) % boundary_cond(f) % kd   = kd
      dii = 0; dij = 0; dik = 0
      dji = 0; djj = 0; djk = 0
      dki = 0; dkj = 0; dkk = 0 
      do nb = 1, nBlock
         found = .false.
         if (b == nb) cycle
         do nf = 1, number_of_face
            do per = 1, number_of_permutation
               found = .true.
               do p = 1, number_of_points_per_face
                  do d = 1,3
                     if (abs(corner_points(d,point_on_face(p                 , f), b) & 
                           - corner_points(d,point_on_face(permutation(p,per),nf),nb) ) > EPS) then
                        found = .false.
                        exit
                     end if
                  end do
                  if (.not. found) then
                     exit
                  end if
               end do
               if (found) then
                  blocks(b) % boundary_cond(f) % bc_type     = nb
                  blocks(b) % boundary_cond(f) % face        = nf
                  blocks(b) % boundary_cond(f) % permutation = per
                  if (f == EAST .and. nf == WEST) then
                     select case (per) 
                     case(1)
                        ci = 1 - blocks(b) % nPoints(1); dii = 1
                        cj = 0                         ; djj = 1
                        ck = 0                         ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == EAST .and. nf == SOUTH) then
                     select case (per) 
                     case(2)
                        ci = 1 + blocks(b) % nPoints(2) ; dij = -1
                        cj = 1 - blocks(b) % nPoints(1) ; dji = 1
                        ck = 0                          ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == WEST .and. nf == EAST) then
                     select case (per) 
                     case(1)
                        ci = blocks(nb) % nPoints(1) - 1; dii = 1
                        cj = 0                          ; djj = 1
                        ck = 0                          ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == WEST .and. nf == SOUTH) then
                     select case (per) 
                     case(1)
                        ci = 0                          ; dij = 1
                        cj = 2                          ; dji = -1
                        ck = 0                          ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == NORTH .and. nf == SOUTH) then
                     select case (per) 
                     case(1)
                        ci = 0                         ; dii = 1
                        cj = 1 - blocks(b) % nPoints(2); djj = 1
                        ck = 0                         ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == NORTH .and. nf == NORTH) then
                     select case (per) 
                     case(2)
                        ci = blocks(nb) % nPoints(1)+1; dii = -1
                        cj = blocks(nb) % nPoints(2)              &
                           + blocks(b) % nPoints(2)
                                                        djj = -1
                        ck = 0                        ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == SOUTH .and. nf == NORTH) then
                     select case (per) 
                     case(1)
                        ci = 0                        ; dii = 1
                        cj = blocks(nb) % nPoints(2)-1; djj = 1
                        ck = 0                        ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == SOUTH .and. nf == SOUTH) then
                     select case (per) 
                     case(2)
                        ci = blocks(nb) % nPoints(1)+1; dii = -1
                        cj = 2                        ; djj = -1
                        ck = 0                        ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == SOUTH .and. nf == EAST) then
                     select case (per) 
                     case(2)
                        ci = blocks(nb) % nPoints(1) - 1; dij = 1
                        cj = blocks(nb) % nPoints(2) + 1; dji = -1
                        ck = 0                          ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == SOUTH .and. nf == WEST) then
                     select case (per) 
                     case(1)
                        ci = 2                          ; dij = -1
                        cj = 0                          ; dji = 1
                        ck = 0                          ; dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == SOUTH .and. nf == BACK) then
                     select case (per) 
                     case(3)
                        ci = 0                          ; dii = 1
                        cj = blocks(nb) % nPoints(2) + 1; djk = -1
                        ck = blocks(nb) % nPoints(3) - 1; dkj = 1
                     case default 
                        goto 666
                     end select
                  else if (f == BACK .and. nf == FRONT) then
                     select case (per) 
                     case(1)
                        ci = 0                         ; dii = 1
                        cj = 0                         ; djj = 1
                        ck = 1 - blocks(b) % nPoints(3); dkk = 1
                     case default 
                        goto 666
                     end select
                  else if (f == BACK .and. nf == SOUTH) then
                     select case (per) 
                     case(3)
                        ci = 0                          ; dii = 1
                        cj = 1 - blocks(b) % nPoints(3) ; djk = 1
                        ck = 1 + blocks(b) % nPoints(2) ; dkj = -1
                     case default 
                        goto 666
                     end select
                  else if (f == FRONT .and. nf == BACK) then
                     select case (per) 
                     case(1)
                        ci = 0                          ; dii = 1
                        cj = 0                          ; djj = 1
                        ck = blocks(nb) % nPoints(3) - 1; dkk = 1
                     case default 
                        goto 666
                     end select
                  else
                     write(*,*) "Case ",FACE_NAMES(f)  &
                               ," and ",FACE_NAMES(nf) &
                               , per &
                               ," not implemented yet: ",__FILE__,__LINE__
                     stop 1
                  end if

                  blocks(b) % boundary_cond(f) % ci   = ci
                  blocks(b) % boundary_cond(f) % dii  = dii
                  blocks(b) % boundary_cond(f) % dij  = dij
                  blocks(b) % boundary_cond(f) % dik  = dik

                  blocks(b) % boundary_cond(f) % cj   = cj
                  blocks(b) % boundary_cond(f) % dji  = dji
                  blocks(b) % boundary_cond(f) % djj  = djj
                  blocks(b) % boundary_cond(f) % djk  = djk

                  blocks(b) % boundary_cond(f) % ck   = ck
                  blocks(b) % boundary_cond(f) % dki  = dki
                  blocks(b) % boundary_cond(f) % dkj  = dkj
                  blocks(b) % boundary_cond(f) % dkk  = dkk

                  do k = ks,ke
                     do j = js,je
                        do i = is,ie
                           ni = ci + dii * i + dij * j + dik * k
                           nj = cj + dji * i + djj * j + djk * k
                           nk = ck + dki * i + dkj * j + dkk * k
                           call addSamePoint(b,i,j,k,nb,ni,nj,nk)
                        end do
                     end do
                  end do
                  exit
               end if
            end do
            if (found) then
               exit
            end if
         end do
         if (found) then
            exit
         end if
      end do
   end do
end do
return

666   continue
      write(*,*) "Per: ",per, "for ",FACE_NAMES(f)," and ",FACE_NAMES(nf)," not implemented yet: ",__FILE__,__LINE__
      stop 1
   
end subroutine connect_blocks

subroutine addSamePoint(b,i,j,k,nb,ni,nj,nk)
implicit none
integer, intent(in) :: b,i,j,k,nb,ni,nj,nk

integer :: nsp_b,nsp_nb,n,nsp_t
type(t_same), allocatable :: sp_b(:)
type(t_same)              :: temp   

nsp_b = blocks(b) % nSamePoints(i,j,k)
nsp_nb = blocks(nb) % nSamePoints(ni,nj,nk)
! check if already in list
do n = 1, nsp_b
   if (nb == blocks(b) % SamePoints(n,i,j,k) % b .and. & 
       ni == blocks(b) % SamePoints(n,i,j,k) % i .and. & 
       nj == blocks(b) % SamePoints(n,i,j,k) % j .and. & 
       nk == blocks(b) % SamePoints(n,i,j,k) % k ) then
      !write(*,*) b,i,j,k," is already connected to ",nb,ni,nj,nk,n
      return
   end if
end do
!write(*,*) "Working @",b,i,j,k,nsp_b,nsp_nb

if (nsp_b > 0) then
   allocate(sp_b(nsp_b))
   do n =  1, nsp_b 
      sp_b(n) = blocks(b) % SamePoints(n,i,j,k)
   end do
end if
   
if (nsp_nb > 0) then
   do n =  1, nsp_nb 
      temp = blocks(nb) % SamePoints(n,ni,nj,nk)
      nsp_b = nsp_b + 1
      blocks(b) % SamePoints(nsp_b,i,j,k) = temp
      nsp_t = blocks(temp % b) % nSamePoints(temp % i,temp % j,temp % k)
      nsp_t = nsp_t + 1
      blocks(temp % b) % SamePoints(nsp_t,temp % i,temp % j, temp % k) = t_same(b,i,j,k)
      blocks(temp % b) % nSamePoints(temp % i,temp % j,temp % k) = nsp_t
   end do
end if

if (allocated(sp_b)) then
   do n =  1, ubound(sp_b,1) 
      nsp_nb = nsp_nb + 1 
      blocks(nb) % SamePoints(nsp_nb,ni,nj,nk) = sp_b(n)
      nsp_t = blocks(sp_b(n) % b) % nSamePoints(sp_b(n) % i,sp_b(n) % j,sp_b(n) % k)
      nsp_t = nsp_t + 1
      blocks(sp_b(n) % b) % SamePoints(nsp_t,sp_b(n) % i,sp_b(n) % j, sp_b(n) % k) = t_same(nb,ni,nj,nk)
      blocks(sp_b(n) % b) % nSamePoints(sp_b(n) % i,sp_b(n) % j,sp_b(n) % k) = nsp_t
   end do
   deallocate(sp_b)
end if

nsp_b = nsp_b + 1
!write(*,*) "+Adding:",nb,ni,nj,nk," to ",b,i,j,k,"@",nsp_b
blocks(b) % SamePoints(nsp_b,i,j,k) = t_same(nb,ni,nj,nk)

nsp_nb = nsp_nb + 1
!write(*,*) "-Adding:",b,i,j,k," to ",nb,ni,nj,nk,"@",nsp_nb
blocks(nb) % SamePoints(nsp_nb,ni,nj,nk) = t_same(b,i,j,k)

! Updating the Number of Same points
blocks(b) % nSamePoints(i,j,k) = nsp_b
blocks(nb) % nSamePoints(ni,nj,nk) = nsp_nb
end subroutine addSamePoint

subroutine write_grid_hdf5()
implicit none

integer, parameter                :: RANK = 3
character(len=len(GROUP_BLOCK)+2) :: block_group

integer     ::   error ! Error flag
integer(hid_t) :: file_id       ! file identifier
integer(hid_t) :: group_id_grid      ! dataset identifier
integer(hid_t) :: group_id_block      ! dataset identifier
integer(hid_t) :: dset_id       ! dataset identifier
integer(hid_t) :: dspace_id     ! dataspace identifier
integer(HSIZE_T) :: dims(3)
integer :: b,d

! Initialize FORTRAN interface.
nCell = 0
call h5open_f(error)

! create a new file using default properties.
call h5fcreate_f(filename_grid_out, h5f_acc_trunc_f, file_id, error)

! Create a group grid in the file.
call h5gcreate_f(file_id,  GROUP_GRID, group_id_grid,  error)

!   Write Coordinates into grid/blockN/
do b = 1, nBlock
   write(block_group,'(A,I0)') GROUP_BLOCK,b
   dims  = blocks(b) % nPoints 

   ! create the dataspace.
   call h5screate_simple_f(RANK, dims, dspace_id, error)

   ! Create a group named for block1 in the file.
   call h5gcreate_f(group_id_grid, block_group, group_id_block, error)
   
   do d = 1, RANK
      call h5dcreate_f(group_id_block, COORD_NAME(d), h5t_native_double, dspace_id, &
                      dset_id, error)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, blocks(b) % coords(:,:,:,d), dims, error)
      call h5dclose_f(dset_id, error)
   end do

   ! terminate access to the data space of the grid and create a new one for the data
   call h5sclose_f(dspace_id, error)
   ! Close the group.
   call h5gclose_f(group_id_block, error)
end do
call h5gclose_f(group_id_grid, error) ! CLOSE GRID GROUP

call h5fclose_f(file_id, error)
end subroutine write_grid_hdf5

subroutine write_grid_ufo()
implicit none
integer :: fu
integer :: b,i,j,k
integer :: dimen_temp

open(newunit=fu,file=trim(filename_grid_out),form="unformatted",access="stream")

if (dimen == 3) then
   dimen_temp = 2
else
   dimen_temp = 0
end if
write(fu) dimen_temp, nBlock
do b = 1, nBlock
  write(fu) blocks(b) % nPoints(1:dimen)
end do
do b = 1, nBlock
   write(fu) (((blocks(b) % coords (i,j,k,1:dimen),i = 1, blocks(b) % nPoints(1)) &
                                                  ,j = 1, blocks(b) % nPoints(2)) & 
                                                  ,k = 1, blocks(b) % nPoints(3))  

end do
close(fu)
end subroutine write_grid_ufo
end module structured_grid  
