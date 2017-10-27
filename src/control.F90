module control
   use const
implicit none
logical :: end_adaption
logical :: output_grid
integer :: nIter = 5000
integer :: res_out = 10
integer :: res_out_start = 50
integer :: iter = 0
real(REAL_KIND) :: stop_con_res = 1.0E-12_REAL_KIND
contains

subroutine loop_control(max_res)
implicit none
real(REAL_KIND) , intent(in) :: max_res

iter = iter + 1

if (iter >= nIter) then
   end_adaption = .true.
   write(*,*) "Maximum Number of Iterations"
end if

if (max_res <= stop_con_res .and. iter >= 10) then
   end_adaption = .true.
   write(*,*) "Convergence Achieved", iter
end if

end subroutine loop_control

end module control
