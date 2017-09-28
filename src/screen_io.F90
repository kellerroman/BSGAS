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

subroutine sw_grid_info(blocks)
use types
implicit none
type(t_block), intent(in) :: blocks(:)
integer :: nBlock,b

nBlock = ubound(blocks,1)
!write(*,*) "Number of Blocks:", nBlock
write(*,'(A3,9(1X,A5))') "B#", "NI","NJ","NK","WEST","EAST","SOUTH","NORTH","FRONT","BACK"
do b = 1, nBlock
   write(*,'(I3,9(1X,I5))') b, blocks(b) % nCells, blocks(b) % boundary_cond(:) % bc_type
end do
end subroutine sw_grid_info

subroutine sw_residual(iter,max_spring, min_spring &
                      ,max_edge_f,max_point_f, sum_point_f &
                      ,max_edge_len,min_edge_len)
   use control, only:res_out_start, res_out
implicit none
integer, intent(in) :: iter
real(REAL_KIND), intent(in) :: max_spring, min_spring &
                             , max_point_f, max_edge_f, sum_point_f &
                             , max_edge_len,min_edge_len

if (iter <= res_out_start .or. mod(iter,res_out) == 0) then
   write(*,'(I10,7(1X,ES10.3))') iter &
            , max_spring, min_spring & 
            , max_edge_f,max_point_f,sum_point_f &
            , max_edge_len,min_edge_len
end if

end subroutine sw_residual

subroutine sw_init_residual
implicit none
write(*,'(8(A10,1X))') "ITERATION","MAX SPRING", "MIN SPRING" &
                     , "F MAX EDGE","F MAX POINT","F AVG POINT" &
                     , "LEN MAX","LEN MIN"
end subroutine sw_init_residual

subroutine sw_edge_info(e)
   use unstr, only: git
   implicit none
   integer, intent(in) :: e
   integer :: ne
   write(*,'("Edge: ",I0)') e
   write(*,'(2X,"From",4(1X,I0)," to",4(1X,I0))') git % point_refs(:,git % edge_points(1,e)) &
                                              ,git % point_refs(:,git % edge_points(2,e))
   ne = git % edge_neighbor(1,e)
   write(*,'(2X,"Neighbors: ",I0)') ne
   write(*,'(4X,"From",4(1X,I0)," to",4(1X,I0))') git % point_refs(:,git % edge_points(1,ne)) &
                                              ,git % point_refs(:,git % edge_points(2,ne))
   ne = git % edge_neighbor(2,e)
   write(*,'(2X,"Neighbors: ",I0)') ne
   write(*,'(4X,"From",4(1X,I0)," to",4(1X,I0))') git % point_refs(:,git % edge_points(1,ne)) &
                                              ,git % point_refs(:,git % edge_points(2,ne))

   end subroutine
end module screen_io
