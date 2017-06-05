module structured_grid
use hdf5
use const
use types
use help_routines, only: alloc
implicit none

character(len=100)          :: filename_grid_in
character(len=*), parameter :: GROUP_GRID            = "grid"
character(len=*), parameter :: GROUP_BLOCK           = "block"
character(len=*), parameter :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]

integer(INT_KIND)           :: nBlock
type(t_block), allocatable  :: blocks(:)
integer(INT_KIND)           :: dimen

contains

subroutine read_grid()
implicit none
call read_grid_hdf5()
call connect_blocks()
end subroutine read_grid

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
logical :: fexists

inquire(file=trim(filename_grid_in),exist=fexists)
if (.not. fexists) then
   write(*,*) "Grid file: '"//trim(filename_grid_in)//"' not found!",__FILE__,__LINE__
   stop 1
end if
! Initialize FORTRAN interface.
call h5open_f(error)

! Open an existing file.
call h5fopen_f (filename_grid_in, h5f_acc_rdwr_f, file_id, error)

call h5gopen_f(file_id,GROUP_GRID,group_id_grid,error)

call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)
allocate(blocks(nBlock))
write(*,*) "Number of Blocks:", nBlock
write(*,'(A3,3(1X,A3))') "B#", "NI","NJ","NK"
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
   !nCell = nCell + product(blocks(b) % nCells)
   
   !!!! MUSS VERSCHOBEN WERDEN FUER MULTIBLOCK
!   call allocate_vars(b)
   write(*,'(I3,3(1X,I3))') b, blocks(b) % nCells

   call alloc (data_in, blocks(b)%nPoints)
   allocate ( blocks(b) % coords (blocks(b)%nPoints(1),blocks(b)%nPoints(2),blocks(b)%nPoints(3),3))
   call alloc (blocks(b) % refs,blocks(b) % nPoints)
   blocks(b) % refs = -1

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

subroutine connect_blocks()
implicit none
real(REAL_KIND),parameter :: EPS = 1.0e-8_REAL_KIND
integer :: b,f,nb,nf,p,per,d
integer :: ni,nj,nk
logical :: is3D 
integer :: number_of_corner_point
integer :: number_of_face
integer :: number_of_points_per_face
integer :: number_of_permutation

logical :: found


real(REAL_KIND), allocatable :: corner_points(:,:,:)
!< Array of the Cornerpoints for each block (dim,point_id,block)
integer, allocatable :: point_on_face(:,:)
!< Array with the points for each face (point_id,face)
!< face: (W,E,S,N,F,B)
integer, allocatable :: permutation(:,:)
!< Different Point arrangment due to misaligned blocks (point_id,permutation)
b = 1
if (blocks(b) % nPoints(2) == 1) then
   return 
else if  ( blocks(b) % nPoints(3) == 1) then
   is3D = .false.
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

   allocate( blocks(b) % boundary_cond(number_of_face))

   blocks(b) % boundary_cond (:) % bc_type         = 0
   blocks(b) % boundary_cond (:) % neighbor_face   = 0
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
                  !write(*,*) "found:",b,f,nb,nf,per
                  blocks(b) % boundary_cond(f) % bc_type       = nb
                  blocks(b) % boundary_cond(f) % neighbor_face = nf
                  blocks(b) % boundary_cond(f) % permutation   = per
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

end subroutine connect_blocks

end module structured_grid  
