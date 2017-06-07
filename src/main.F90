!================================================================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BLOCK STRUCUTRED GRID ADAPTION SOLVER MAIN FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===                                                                                                                          ===
!! AUTHORS:             ROMAN KELLER(RK)
!===                                                                                                                          ===
!! START:               13.05.2017               
!! LAST CHANGE:         05.06.2017
!! VERSION:             V0.0.1
!===                                                                                                                          ===
!! CHANGELOG:
!!          13.05.2017,RK: Start of Project
!===                                                                                                                          ===
!! TODO:
!===                                                                                                                          ===
!================================================================================================================================

program BSGAS
use boundary, only: read_boundary, init_boundary, init_walledges
use control, only: loop_control, end_adaption, iter
use config, only: read_config
use screen_io, only: sw_program_start,sw_program_end, sw_init_residual, sw_residual, sw_grid_info
use structured_grid, only: read_grid, blocks
use unstr, only: strukt2unstr,calc_edge_length,calc_edge_forces, calc_point_forces, move_points, git
use unstr_output, only: write_output
use spring, only: init_springs,calc_edge_springs
use types
implicit none
real(REAL_KIND) :: max_point_f, max_edge_f

call sw_program_start()

call read_config

call read_grid
call read_boundary(blocks)
call sw_grid_info(blocks)

call strukt2unstr(blocks)
call init_boundary(git,blocks)
call init_walledges(git,blocks)

call calc_edge_length()
call init_springs
call sw_init_residual

call write_output(0)

end_adaption = .false.
ITER_LOOP: do while (.not. end_adaption)

   call loop_control()

   call calc_edge_length()

   call calc_edge_springs()

   call calc_edge_forces(max_edge_f)

   call calc_point_forces(max_point_f)

   call sw_residual(iter,max_edge_f,max_point_f)

   call move_points(max_edge_f)
   
   call write_output(iter)
end do ITER_LOOP

!call unstr2strukt

!call write_grid
call sw_program_end()
end program BSGAS
