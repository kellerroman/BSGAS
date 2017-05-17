!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BLOCK STRUCUTRED GRID ADAPTION SOLVER SCREEN IO MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    this module is intendet to manage all output to the screen and in future log file outputs
!!    ALL routines that can be called are named sw_*
!!    
!!
!! AUTHORS:             ROMAN KELLER(RK)
!!
!! START:               13.05.2017               
!! LAST CHANGE:         13.05.2017
!! VERSION:             V0.0.1
!!
!! CHANGELOG:
!!          13.05.2017,RK: Start of Project
!!
!! TODO:
!!          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module screen_io
use, intrinsic :: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT
use const
implicit none
integer(INT_KIND), parameter :: SCREEN_WIDTH = 130

character(len=*),parameter :: RED_START            = achar(27)//"[31m"
character(len=*),parameter :: RED_END              = achar(27)//"[0m"                    
character(len=*),parameter :: GREEN_START          = achar(27)//"[32m"
character(len=*),parameter :: GREEN_END            = achar(27)//"[0m"                    
contains
subroutine sw_program_start()
implicit none
call sw_full_line("")
call sw_full_line("")
call sw_full_line("START OF BSGAS")
call sw_full_line("")
call sw_full_line("")
end subroutine sw_program_start


subroutine sw_program_end()
implicit none
call sw_full_line("")
call sw_full_line("")
call sw_full_line("END OF BSGAS")
call sw_full_line("")
call sw_full_line("")
end subroutine sw_program_end


subroutine sw_full_line(str)
implicit none
character(len=*), intent(in) :: str

integer :: length
! if empty strin make an entire line of "="
if (len_trim(str) == 0) then
   write(stdout,'(A)',ADVANCE='NO') repeat("=",SCREEN_WIDTH)
else
   length = SCREEN_WIDTH - 2 - len_trim(str)
   write(stdout,'(A)'      ,ADVANCE="NO") repeat("=",length/2)
   write(stdout,'(1X,A,1X)',ADVANCE="NO") trim(str)
   write(stdout,'(A)'      ,ADVANCE="NO") repeat("=",length/2)
   if (mod(length,2) /= 0 ) write(stdout,'("=")',ADVANCE="NO")
end if
write(stdout,*)
end subroutine sw_full_line

subroutine sw_residual(iter,max_edge_f,max_point_f)
   use control, only:res_out_start, res_out
implicit none
integer, intent(in) :: iter
real(REAL_KIND), intent(in) :: max_point_f, max_edge_f

if (iter <= res_out_start .or. mod(iter,res_out) == 0) then
   write(*,*) iter, max_edge_f,max_point_f
end if

end subroutine sw_residual

subroutine sw_init_residual
implicit none
write(*,*) "ITERATION","F MAX EDGE","F MAX POINT"
end subroutine sw_init_residual
end module screen_io
