module structured_grid
use hdf5
use const
use types
implicit none

character(len=*), parameter :: filename_grid_in = "grid.h5"
character(len=*), parameter :: GROUP_GRID            = "grid"
character(len=*), parameter :: GROUP_BLOCK           = "block"
character(len=*), parameter :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]

integer(INT_KIND) :: nBlock
type(t_block), allocatable :: blocks(:)
integer(INT_KIND) :: dimen

contains

subroutine read_grid()
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

   allocate (data_in(blocks(b)%nPoints(1),blocks(b)%nPoints(2),blocks(b)%nPoints(3)))
   allocate ( blocks(b) % coords (blocks(b)%nPoints(1),blocks(b)%nPoints(2),blocks(b)%nPoints(3),3))
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
end subroutine read_grid
end module structured_grid  
