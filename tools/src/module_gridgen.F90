module mod_gridgen
use hdf5
implicit none
private
integer         , parameter :: MAX_BLOCK             = 10 !Maximale Anzahl Bloecke (Dimension des temp block arrays)
integer         , parameter :: VARNAME_LENGTH        = 20
character(len=*), parameter :: FILENAME              = "grid.h5"  ! file name HDF5
character(len=*), parameter :: GROUP_GRID            = "grid"
character(len=*), parameter :: GROUP_DATA            = "data"
character(len=*), parameter :: GROUP_BLOCK           = "block"

character(len=*) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]

integer :: dimension
integer :: nblock = 0
integer :: nVar = 4
integer, allocatable :: ncells_temp(:,:)
logical :: blocks_allocated = .false.
logical :: debug = .false.
type :: tblock
   integer                       :: ncells(3)
   integer                       :: npkts(3)
   real(kind = 8), allocatable   :: xyzs(:,:,:,:)
end type
type(tblock), allocatable :: blocks(:)
interface add_block
   module procedure add_block_scalar1
   module procedure add_block_scalar2
   module procedure add_block_scalar3
   module procedure add_block_array3
end interface add_block

public:: blocks,  add_block, allocate_blocks, debug, write_grid, set_dimension, tblock

contains
   subroutine set_dimension(given_dimension)
      implicit none
      integer, intent(in) :: given_dimension

      dimension = given_dimension

      allocate(ncells_temp(MAX_BLOCK,dimension) )
   end subroutine set_dimension

   subroutine add_block_scalar1(ni)
      implicit none
      integer, intent(in) :: ni
      if (dimension /= 1) then
         write(*,*) "Cannot add 1D Block"
         write(*,*) "Dimension was set to:",dimension
         stop 1
      end if
      if (blocks_allocated) then
         write(*,*) "add_block: Cannot add Block"
         write(*,*) "Blocks are allready allocated"
         stop 1
      end if
      nblock = nblock + 1
      ncells_temp(nblock,1)  = ni
      if (debug) &
         write(*,*) "Adding Block: NB:",nblock,"(ni)",ni
   end subroutine

   subroutine add_block_scalar2(ni,nj)
      implicit none
      integer, intent(in) :: ni,nj
      if (dimension /= 2) then
         write(*,*) "Cannot add 2D Block"
         write(*,*) "Dimension was set to:",dimension
         stop 1
      end if
      if (blocks_allocated) then
         write(*,*) "add_block: Cannot add Block"
         write(*,*) "Blocks are allready allocated"
         stop 1
      end if
      nblock = nblock + 1
      ncells_temp(nblock,1)  = ni
      ncells_temp(nblock,2)  = nj
      if (debug) &
         write(*,*) "Adding Block: NB:",nblock,"(ni,nj)",ni,nj
   end subroutine

   subroutine add_block_array3(ncells)
      implicit none
      integer, intent(in) :: ncells(3)
      call add_block(ncells(1),ncells(2),ncells(3))
   end subroutine

   subroutine add_block_scalar3(ni,nj,nk)
      implicit none
      integer, intent(in) :: ni,nj,nk
      if (dimension /= 3) then
         write(*,*) "Cannot add 3D Block"
         write(*,*) "Dimension was set to:",dimension
         stop 1
      end if
      if (blocks_allocated) then
         write(*,*) "add_block: Cannot add Block"
         write(*,*) "Blocks are allready allocated"
         stop 1
      end if
      nblock = nblock + 1
      ncells_temp(nblock,1)  = ni
      ncells_temp(nblock,2)  = nj
      ncells_temp(nblock,3)  = nk
      if (debug) &
         write(*,*) "Adding Block: NB:",nblock,"(ni,nj,nk)",ni,nj,nk
   end subroutine

   subroutine allocate_blocks(anzahl_var)
      implicit none
      integer, intent(in) :: anzahl_var
      integer :: i
      if (blocks_allocated) then
         write(*,*) "allocate_blocks: Cannot allocate Blocks"
         write(*,*) "Blocks are allready allocated"
         stop 1
      end if
      if (nblock > 0) then
         blocks_allocated = .true.
         !if (anzahl_var < 5) then
         !   write(*,*) "allocate_blocks: Cannot allocate Blocks"
         !   write(*,*) "Passed nVar is to small",anzahl_var
         !else
            nvar = anzahl_var
         !end if

         if (debug) &
            write(*,*) "Allocate Blocks",nblock,"nVar:",nVar

         allocate( blocks(nblock))
         do i = 1,nblock
            blocks(i) % nCells = 1
            blocks(i) % nPkts  = 1

            blocks(i) % nCells(1:dimension) = ncells_temp(i,:)
            blocks(i) % nPkts (1:dimension) = ncells_temp(i,:) + 1
            if (debug) then
               write(*,*) i, blocks(i)%ncells, blocks(i) % npkts
            end if
            allocate(blocks(i) % xyzs(blocks(i) %  npkts(1),blocks(i) %  npkts(2),blocks(i) %  npkts(3),3))
         end do
      else
         write(*,*) "allocate_blocks: Cannot allocate Blocks"
         write(*,*) "no Blocks defined"
         stop 1
      end if
   end subroutine
   subroutine write_grid()
   implicit none

   integer, parameter                :: RANK = 3
   integer(hsize_t), dimension(RANK) :: DIMS
   integer     ::   error ! error flag

   integer(hid_t) :: file_id       ! file identifier
   integer(hid_t) :: dset_id       ! dataset identifier
   integer(hid_t) :: group_id      ! dataset identifier
   integer(hid_t) :: group_id2     ! dataset identifier
   integer(hid_t) :: dspace_id     ! dataspace identifier
      
   character(len=7) :: block_group

   integer :: nb 
   integer :: nd
   
   if (debug) &
      write(*,*) "write_grid: Writing HDF5 output File"

   ! initialize fortran interface.
   call h5open_f(error)
   
   ! create a new file using default properties.
   call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)

   ! Create a group grid in the file.
   call h5gcreate_f(file_id,  GROUP_GRID, group_id,  error)

   !   Write Coordinates into grid/blockN/
   do nb = 1, nBlock
      write(block_group,'(A,I0)') GROUP_BLOCK,nb
      dims  = blocks(nb) % npkts 

      ! create the dataspace.
      call h5screate_simple_f(rank, dims, dspace_id, error)
   
      ! Create a group named for block1 in the file.
      call h5gcreate_f(group_id, block_group, group_id2, error)
      
      do nd = 1, RANK
         call h5dcreate_f(group_id2, COORD_NAME(nd), h5t_native_double, dspace_id, &
                         dset_id, error)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, blocks(nb) % xyzs(:,:,:,nd), dims, error)
         call h5dclose_f(dset_id, error)
      end do

      ! terminate access to the data space of the grid and create a new one for the data
      call h5sclose_f(dspace_id, error)
      ! Close the group.
      call h5gclose_f(group_id2, error)
   end do

   ! Close the group.
   call h5gclose_f(group_id, error)

   ! close the file.
   call h5fclose_f(file_id, error)
   
   ! close fortran interface.
   call h5close_f(error)

   end subroutine

end module
