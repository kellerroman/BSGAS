module config
   use const
   use config_parser
implicit none

character(len=100) :: filename_config = "config.cfg"

contains

subroutine read_config
   use control, only: nIter, res_out,res_out_start
   use spring, only: cell_inc,cell_parallel_inc, faktor_wall, faktor_strech, faktor_para, spring_max
   use structured_grid, only: filename_grid_in,filename_grid_out
   use unstr_output, only: output_intervall
   use unstr, only: wall_move_rest, point_weight

implicit none
   logical :: fexists
   integer :: fu
   integer :: io_stat
   integer :: pos
   character(len=100) :: line
   character(len=VARNAME_LENGTH) :: varname
   character(len=90) :: varvalue
   character(len=VARNAME_LENGTH),allocatable :: unset_paras(:)

call add_parameter("GRID_OUT"          ,filename_grid_out      ,"grid_out.h5"     )
call add_parameter("GRID_IN"           ,filename_grid_in       ,"grid.h5"         )
call add_parameter("SPRING_MAX"        ,spring_max             ,1.00E+04_REAL_KIND)
call add_parameter("FAKTOR_PARA"       ,faktor_para            ,1.00E-03_REAL_KIND)
call add_parameter("FAKTOR_STRECH"     ,faktor_strech          ,1.00E-05_REAL_KIND)
call add_parameter("FAKTOR_WALL"       ,faktor_wall            ,1.00E+00_REAL_KIND)
call add_parameter("POINT_WEIGHT"      ,point_weight           ,2.50E-00_REAL_KIND)
call add_parameter("WALL_MOVE_REST"    ,wall_move_rest                            )
call add_parameter("CELL_INC"          ,cell_inc               ,1.25E+00_REAL_KIND)
call add_parameter("CELL_PARALLEL_INC" ,cell_parallel_inc      ,1.25E+00_REAL_KIND)
call add_parameter("OUTPUT_INTERVALL"  ,output_intervall       ,10                )
call add_parameter("RES_OUT_START"     ,res_out_start          ,50                )
call add_parameter("RES_OUT"           ,res_out                ,100               )
call add_parameter("NITER"             ,niter                                     )

inquire(file=trim(filename_config),exist=fexists)
if (.not. fexists) then
   write(*,*) "Boundary file: '"//trim(filename_config)//"' not found!",__FILE__,__LINE__
   stop 1
end if
open(newunit=fu,file=trim(filename_config))

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
   !write(*,*) line
   pos = index(line,"=")

   if (pos < 1) then
      write(*,'("konnte line nicht verarbeiten, kein = gefunden:",a)') trim(line)
      stop 1
   end if
!   write(*,*) trim(line(1:pos-1)),trim(line(pos+1:))
   varname = trim(line(1:pos-1))
   varvalue = trim(adjustl(line(pos+1:)))
   select case(set_para(varname,varvalue))
   case(RETURN_CODE_NOT_IN_LIST)
      write(*,*) "'"//trim(varname)//"' is not in list"
      stop 1
   case(RETURN_CODE_NOT_A_INT)
      write(*,*) "value for '"//trim(varname)//"' is not an INTEGER ",varvalue
      stop 1
   case(RETURN_CODE_NOT_A_REAL)
      write(*,*) "value for '"//trim(varname)//"' is not an REAL ",varvalue
      stop 1
   case(RETURN_CODE_ALREADY_SET)
      write(*,*) "value for '"//trim(varname)//"' is already SET"
      stop 1
   end select
   
end do
close(fu)
call get_unset_paras(unset_paras)
if (ubound(unset_paras,1) > 0 ) then
   write(*,*) "There are unset Parameters which require Values"
   do pos = 1,ubound(unset_paras,1)
      write(*,*) pos, unset_paras(pos)
   end do
   stop 1
end if

call list_parameter()
call free_all()

! Parameter calculation
point_weight = 1.0E+00_REAL_KIND / point_weight
end subroutine read_config

end module config
