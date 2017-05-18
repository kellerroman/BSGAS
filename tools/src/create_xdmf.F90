   program create_xmdf
      use hdf5
   implicit none
   character(len=100)              :: filename_sol = "grid.h5"
   character(len=100)              :: filename_xdmf = "grid.xdmf"

   character(len=*), parameter :: GROUP_GRID            = "grid"
   character(len=*), parameter :: GROUP_BLOCK           = "block"
   
   character(len=*) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]


   integer     ::   error ! Error flag
   integer(hid_t) :: file_id       ! file identifier
   integer(hid_t) :: dset_id       ! dataset identifier
   integer(hid_t) :: group_id      ! dataset identifier
   integer(hid_t) :: group_id2     ! dataset identifier
   integer(hid_t) :: dspace_id     ! dataspace identifier
   integer(HSIZE_T),allocatable :: dims(:,:)
   integer(HSIZE_T),allocatable :: dims_data(:,:)
   integer(HSIZE_T) :: maxdims(3)
   character(len=len(GROUP_BLOCK)+2) :: block_group
   integer :: fu
   integer :: nBlock

   integer :: nb, i
   integer :: nd,d
   integer:: rank
   logical                         :: fexists

   character(len = 100) :: arg
   
   DO
      CALL get_command_argument(i, arg)
      IF (LEN_TRIM(arg) == 0) EXIT
      if( i == 1) then
         read(arg,*) filename_sol
      else if (i == 2) then
         read(arg,*) filename_xdmf
      end if   
      i = i+1
   END DO
   write(*,*) "create_xdmf: Writing XDMF output File "//trim(filename_xdmf)//" for "//trim(filename_sol)
   inquire(file=trim(filename_sol),exist=fexists)
   
   if(.not. fexists) then
     write(*,*) "Solution File konnte nicht gefunden werden: "//TRIM(filename_sol)
     stop 1
   end if
      
   open(newunit=fu,file=filename_xdmf)
   ! Initialize FORTRAN interface.
   call h5open_f(error)

   ! Open an existing file.
   call h5fopen_f (filename_sol, h5f_acc_rdwr_f, file_id, error)

   call h5gopen_f(file_id,GROUP_GRID,group_id,error)

   call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)

   allocate(dims(3,nBlock))
   allocate(dims_data(3,nBlock))

   write(*,*) "Number of Blocks:", nBlock
   do nb = 1, nBlock

      write(block_group,'(A,I0)') GROUP_BLOCK,nb

      call h5gopen_f(group_id,block_group,group_id2,error)
      
      ! Open an existing dataset.
      call h5dopen_f(group_id2, COORD_NAME(1), dset_id, error)
      call h5dget_space_f(dset_id,dspace_id,error)
      
      call h5sget_simple_extent_ndims_f(dspace_id,RANK,error)
      call h5sget_simple_extent_dims_f(dspace_id,dims(:,nb),maxdims,error)

      call h5dclose_f(dset_id, error)
      call h5gclose_f(group_id2, error)
      write(*,*) nb, dims(:,nb)

   end do
   call h5gclose_f(group_id, error)

!   call h5gopen_f(file_id,GROUP_DATA,group_id,error)
!   call h5gn_members_f(file_id, GROUP_DATA, hdf5_nSol, error)
!   write(*,*) "Anzahl der Bloecke:" ,nBlock
!   write(*,*) "Anzahl der Loesungen:" ,hdf5_nSol
!   allocate(solution_name(hdf5_nSol))
!
!   do nd = 1, hdf5_nSol
!      call h5gget_obj_info_idx_f(file_id, GROUP_DATA, nd-1,solution_name(nd), solution_type, error)
!      call h5gopen_f(group_id,solution_name(nd),group_id1,error) ! OPEN TIMESTEP GROUP
!      if (nd == 1) then
!         do nb = 1, nBlock
!            write(block_group,'(A,I0)') GROUP_BLOCK,nb
!            call h5gopen_f(group_id1,block_group,group_id2,error)   ! open block Group
!            call h5gn_members_f(group_id1, block_group, nVar_in, error)
!            if (nd == 1 .and. nb == 1) then
!               allocate(varName_in(nVar_in))
!               do var = 1, nVar_in
!                  call h5gget_obj_info_idx_f(group_id1, block_group, var-1,varName_in(var), var_type, error)
!                  call h5dopen_f(group_id2, varName_in(var), dset_id, error)
!                  call h5dget_space_f(dset_id,dspace_id,error)
!                  
!                  call h5sget_simple_extent_dims_f(dspace_id,dims_data(:,nb),maxdims,error)
!
!                  call h5dclose_f(dset_id, error)
!                  
!               end do
!            end if
!            call h5gclose_f(group_id2, error)
!         end do
!      end if
!      call h5gclose_f(group_id1, error) !CLOSE TIMESTEP GROUP
!   end do
!   call h5gclose_f(group_id, error)  !CLOSE DATA GROUP

   ! close the file.
   call h5fclose_f(file_id, error)
   ! close fortran interface.
   call h5close_f(error)

   write(fu,'(A)') '<?xml version="1.0" ?>'
   write(fu,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
   write(fu,'(A)') '<Xdmf Version="2.0">'
   write(fu,'(A)') ' <Domain>'
!   if (hdf5_nSol > 1) then
!      write(fu,'(A)') '   <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
!      write(fu,'(A)') '   <Time TimeType="List">'
!      write(fu,'(A,I0,A)',advance="no") '     <DataItem Format="XML" Dimensions="',hdf5_nSol,'">'
!      do ns = 1, hdf5_nSol
!         write(fu,'(1X,A)',advance="no") trim(solution_name(ns))
!      end do
!      write(fu,'(A)')  '</DataItem>'
!      write(fu,'(A)') '   </Time>'
!   end if
   !do ns = 1, hdf5_nSol
      write(fu,'(A)') '   <Grid Name="All" GridType="Tree">'
      do nb = 1, nBlock

         write(block_group,'(A,I0)') GROUP_BLOCK,nb
         write(fu,'(A,I0,A)') '   <Grid Name="Block',nb,'" GridType="Uniform">'
         if (dims(2,nb) ==1) then
            write(fu,'(A,I0,A)') '<Topology TopologyType="Quadrilateral" NumberOfElements="',dims(1,nb),'">'
            write(fu,'(A,I0,A)') '<DataItem Format="XML" DataType="Int" Dimensions="',dims(1,nb),' 1">'
            do d = 0,int(dims(1,nb)-1,4)
               write(fu,'(I0,1X)') d
            end do
            write(fu,'(A)') '</DataItem>'
            write(fu,'(A)') '</Topology>'
         else if (dims(3,nb) == 1) then
         write(fu,'(A,3(I0,1X),A)') '     <Topology TopologyType="2DSMesh" NumberOfElements="'&
                                   ,(dims(d,nb),d=RANK,1,-1) &
                                   ,'"/>'
         else
         write(fu,'(A,3(I0,1X),A)') '     <Topology TopologyType="3DSMesh" NumberOfElements="'&
                                   ,(dims(d,nb),d=RANK,1,-1) &
                                   ,'"/>'
         end if
         write(fu,'(A)') '     <Geometry GeometryType="X_Y_Z">'
         do nd = 1, RANK
            write(fu,'(A,3(I0,1X),A)') '       <DataItem Dimensions="'&
                                      ,(dims(d,nb),d=RANK,1,-1) &
                                      ,'" NumberType="Float" Precision="7" Format="HDF">'
            write(fu,'(8X,A,":",3("/",A))') trim(filename_sol),group_grid,trim(block_group),coord_name(nd)
            write(fu,'(A)') '       </DataItem>'
         end do
         write(fu,'(A)') '     </Geometry>'

!         do nd = 1, nVar_in
!            write(fu,'(3A)') '     <Attribute Name="',trim(varname_in(nd)),'" Center="Cell">'
!            write(fu,'(A,3(I0,1X),A)') '       <DataItem Dimensions="'&
!                                      ,(dims_data(d,nb),d=1,RANK) &
!                                      ,'" NumberType="Float" Precision="7" Format="HDF">'
!            write(fu,'(8X,A,":",4("/",A))') trim(filename_sol),group_data,trim(solution_name(ns))&
!                                           ,trim(block_group),trim(varname_in(nd))
!            write(fu,'(A)') '       </DataItem>'
!            write(fu,'(A)') '     </Attribute>'
!         end do
         write(fu,'(A)') '   </Grid>'
      end do
      write(fu,'(A)') '   </Grid>'
   !end do
   !if (hdf5_nSol > 1) then
      !write(fu,'(A)') '   </Grid>'
   !end if
   write(fu,'(A)') ' </Domain>'
   write(fu,'(A)') '</Xdmf>'

   close(fu)

   deallocate(dims)
   deallocate(dims_data)
end program
