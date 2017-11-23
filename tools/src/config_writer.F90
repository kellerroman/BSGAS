program config_writer
   use config_parser
   use const
implicit none

character(len = *), parameter :: para_int = "(A20,"" = "",I0)"
character(len = *), parameter :: para_real= '(A20," = ",ES11.4)'
character(len = *), parameter :: para_str = '(A20," =  ",A)'

character(len=100)                        :: filename_config

integer                                   ::   nIter
integer                                   ::   res_out
integer                                   ::   res_out_start
integer                                   ::   output_intervall
integer                                   ::   solution_output
real(kind = 8)                            ::   cell_inc
real(kind = 8)                            ::   cell_parallel_inc
real(kind = 8)                            ::   faktor_wall
real(kind = 8)                            ::   faktor_strech
real(kind = 8)                            ::   faktor_para
real(kind = 8)                            ::   faktor_smooth
real(kind = 8)                            ::   spring_max
real(kind = 8)                            ::   relax_fkt
character(len=100)                        ::   filename_grid_in
character(len=100)                        ::   filename_grid_out

logical                                   :: fexists
logical                                   :: write_config
integer                                   :: fu
integer                                   :: io_stat
integer                                   :: pos
character(len=100)                        :: line
character(len=VARNAME_LENGTH)             :: varname
character(len=90)                         :: varvalue
character(len=VARNAME_LENGTH),allocatable :: unset_paras(:)
integer                                   :: i

filename_config = "config.cfg"
write_config = .FALSE.

call add_parameter("GRID_OUT"          ,filename_grid_out      ,"grid_out.h5"     )
call add_parameter("GRID_IN"           ,filename_grid_in       ,"grid.h5"         )

call add_parameter("RELAXATION_FACTOR" ,relax_fkt              ,9.00E-01_REAL_KIND)

call add_parameter("SPRING_MAX"        ,spring_max             ,1.00E+10_REAL_KIND)

call add_parameter("FAKTOR_PARA"       ,faktor_para            ,1.00E-03_REAL_KIND)
call add_parameter("FAKTOR_STRECH"     ,faktor_strech          ,1.00E-05_REAL_KIND)
call add_parameter("FAKTOR_WALL"       ,faktor_wall            ,1.00E+00_REAL_KIND)
call add_parameter("FAKTOR_SMOOTH"     ,faktor_smooth          ,1.00E-05_REAL_KIND)

call add_parameter("CELL_PARALLEL_INC" ,cell_parallel_inc      ,1.25E+00_REAL_KIND)
call add_parameter("CELL_INC"          ,cell_inc               ,1.25E+00_REAL_KIND)

call add_parameter("OUTPUT_INTERVALL"  ,output_intervall       ,10000             )
call add_parameter("SOLUTION_OUTPUT"   ,solution_output        ,10000             )

call add_parameter("RES_OUT_START"     ,res_out_start          ,50                )
call add_parameter("RES_OUT"           ,res_out                ,10                )
call add_parameter("NITER"             ,niter                  ,10000             )
i = 1
DO
   CALL get_command_argument(i, line)
   IF (LEN_TRIM(line) == 0) EXIT
   if (trim(line) == "-w") then
      write_config = .TRUE.
   else
      read(line,'(A)') filename_config
   end if   
   i = i+1
END DO
write(*,*) "Config file: '"//trim(filename_config)//"'"
if (trim(filename_config) == ".") stop
inquire(file=trim(filename_config),exist=fexists)
if (.not. fexists) then
   write(*,*) "Config file: '"//trim(filename_config)//"' not found!",__FILE__,__LINE__
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
      !write(*,*) "'"//trim(varname)//"' is not in list"
      !stop 1
   case(RETURN_CODE_NOT_A_INT)
      write(*,*) "value for '"//trim(varname)//"' is not an INTEGER ",varvalue
      !stop 1
   case(RETURN_CODE_NOT_A_REAL)
      write(*,*) "value for '"//trim(varname)//"' is not an REAL ",varvalue
      !stop 1
   case(RETURN_CODE_ALREADY_SET)
      !write(*,*) "value for '"//trim(varname)//"' is already SET"
      !stop 1
   end select
   
end do
close(fu)
call get_unset_paras(unset_paras)
if (ubound(unset_paras,1) > 0 ) then
   write(*,*) "There are unset Parameters which require Values"
   do pos = 1,ubound(unset_paras,1)
      write(*,*) pos, unset_paras(pos)
   end do
   !stop 1
end if

!call list_parameter()
call free_all()
if (write_config) then
   open(newunit=fu,file=trim(filename_config))
else
   fu = 6
end if

write(fu,'(A)') ""
write(fu,'(A)') "! **************************************************************************************************"
write(fu,'(A)') "! ***                                                                                            ***"
write(fu,'(A)') "! ***                          Configurationfile for BSGAS                                       ***"
write(fu,'(A)') "! ***                                                                                            ***"
write(fu,'(A)') "! **************************************************************************************************"
write(fu,'(A)') ""

line = "niter"             ;write(fu,para_int) line, niter
line = "res_out"           ;write(fu,para_int ) line, res_out           
line = "res_out_start"     ;write(fu,para_int ) line, res_out_start     
line = "output_intervall"  ;write(fu,para_int ) line, output_intervall  
line = "solution_output"   ;write(fu,para_int ) line, solution_output  
write(fu,'(A)') ""

line = "cell_inc"          ;write(fu,para_real) line, cell_inc          
line = "cell_parallel_inc" ;write(fu,para_real) line, cell_parallel_inc 
write(fu,'(A)') ""

line = "spring_max"        ;write(fu,para_real) line, spring_max        
line = "relaxation_factor" ;write(fu,para_real) line, relax_fkt         
write(fu,'(A)') ""

line = "faktor_wall"       ;write(fu,para_real) line, faktor_wall       
line = "faktor_strech"     ;write(fu,para_real) line, faktor_strech     
line = "faktor_para"       ;write(fu,para_real) line, faktor_para       
line = "faktor_smooth"     ;write(fu,para_real) line, faktor_smooth     
write(fu,'(A)') ""

line = "grid_in"           ;write(fu,para_str ) line, filename_grid_in  
line = "grid_out"          ;write(fu,para_str ) line, filename_grid_out 

if (write_config) close(fu)
end program
