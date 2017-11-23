program BSGAS
! **************************************************************************************************
! ***                          Block Structured Grid Adaption Solver Main File                   ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   13.05.2017
! Last changes: 03.10.2017
! Version:      V0.1.0
! --------------------------------------------------------------------------------------------------
! Description:
!   Transforms Block Structured Grids to obey cell-inc, wall-dist creiteria
!   Currently three reasons for cell size changes are imlemented: wall-dist, neighbor edge length and
!   parallel cell length differences
!   
! --------------------------------------------------------------------------------------------------
! Comments and Notes:
!   
! --------------------------------------------------------------------------------------------------
! References:
!
! --------------------------------------------------------------------------------------------------
! Author and Change History:
!   - 2017-05-13,RK : Started of Project
!   - 2017-09-29,RK : Added ifort as possible compiler
!   - 2017-10-03,RK : Added OpenMP support in gfortran
!   - 2017-10-18,RK : Added Implicit Solver usingn SuperLU
!
! **************************************************************************************************
use boundary, only: read_boundary, init_boundary, init_walledges
use control, only: loop_control, end_adaption, iter
use config, only: read_config
use screen_io, only: sw_program_start,sw_program_end, sw_init_residual, sw_residual, sw_grid_info
use structured_grid, only: read_grid, blocks, write_grid
use unstr, only: strukt2unstr,calc_edge_length, git, unstr2struct
use unstr_output, only: write_output
use spring, only: init_springs,calc_edge_springs
use types
use implicit, only: init_matrix, calc_matrix, solve_system
implicit none
real(REAL_KIND) :: max_edge_len,min_edge_len, avg_edge_len
real(REAL_KIND) :: max_walledge_len,min_walledge_len
real(REAL_KIND) :: max_spring, min_spring
real(REAL_KIND) :: max_res, sum_res
real(REAL_KIND) :: max_spring_res, sum_spring_res

call sw_program_start()

call read_config

call read_grid
call read_boundary(blocks)
call sw_grid_info(blocks)

call strukt2unstr(blocks)
call init_boundary(git,blocks)
call init_walledges(git,blocks)

call init_matrix(git)

call calc_edge_length(max_edge_len,min_edge_len,avg_edge_len,max_walledge_len,min_walledge_len)
call init_springs
!call calc_edge_springs(max_spring,min_spring)

call write_output(0)

call sw_init_residual

end_adaption = .false.
max_res = 1.0D0

ITER_LOOP: do while (.not. end_adaption)

   call loop_control(max_res)

   !call calc_edge_length(max_edge_len,min_edge_len,max_walledge_len,min_walledge_len)

   call calc_edge_springs(max_spring,min_spring,max_spring_res,sum_spring_res)

   call calc_matrix(git)

   call solve_system(max_res,sum_res,git % point_coords)

   call calc_edge_length(max_edge_len,min_edge_len,avg_edge_len,max_walledge_len,min_walledge_len)

   call sw_residual(iter                              &
                   ,max_res, sum_res                  &
                   ,max_spring_res, sum_spring_res    &
                   ,max_spring, min_spring            &
                   ,max_edge_len, min_edge_len, avg_edge_len  &
                   ,max_walledge_len, min_walledge_len)

   call write_output(iter)
   call unstr2struct(blocks,iter)
   call write_grid(iter)
end do ITER_LOOP

call unstr2struct(blocks)
call write_grid()

call sw_program_end()
end program BSGAS
