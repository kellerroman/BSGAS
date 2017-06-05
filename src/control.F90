module control
   use const
implicit none
logical :: end_adaption
logical :: output_grid
integer :: nIter = 5000
integer :: res_out = 10
integer :: res_out_start = 50
integer :: iter = 0
contains

subroutine loop_control
implicit none

iter = iter + 1

if (iter >= nIter) then
   end_adaption = .true.
end if

end subroutine loop_control

end module control
