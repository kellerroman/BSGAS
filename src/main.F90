program BSGAS
! **************************************************************************************************
! ***                          Block Structured Grid Adaption Solver Main File                   ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   13.05.2017
! Last changes: 29.09.2017
! Version:      V0.1.0
! --------------------------------------------------------------------------------------------------
! Description:
!   Fuegt dem 16M Gitter von Markus Lempke weitere Zellen in der
!   sondern 5E-6 ist.
! --------------------------------------------------------------------------------------------------
! Comments and Notes:
!   
! --------------------------------------------------------------------------------------------------
! References:
!
! --------------------------------------------------------------------------------------------------
! Author and Change History:
!   - 2017-05-13: Started of Project by RK
!   - 2017-09-29: Added ifort as possible compiler
!
! **************************************************************************************************
use boundary, only: read_boundary, init_boundary, init_walledges
use control, only: loop_control, end_adaption, iter
use config, only: read_config
use screen_io, only: sw_program_start,sw_program_end, sw_init_residual, sw_residual, sw_grid_info,sw_edge_info
use structured_grid, only: read_grid, blocks, write_grid
use unstr, only: strukt2unstr,calc_edge_length,calc_edge_forces, calc_point_forces, move_points, git, unstr2struct
use unstr_output, only: write_output
use spring, only: init_springs,calc_edge_springs!, springs
use types
implicit none
real(REAL_KIND) :: max_point_f, max_edge_f, sum_point_f,max_edge_len,min_edge_len
real(REAL_KIND) :: max_spring, min_spring

call sw_program_start()

call read_config

call read_grid
call read_boundary(blocks)
call sw_grid_info(blocks)

call strukt2unstr(blocks)
call init_boundary(git,blocks)
call init_walledges(git,blocks)

call calc_edge_length(max_edge_len,min_edge_len)
call init_springs
call calc_edge_springs(max_spring,min_spring)

!write(*,*) blocks(3) % refs(1,1,1), git % point_edges(:,101)
!call sw_edge_info(301)
!write(*,*) blocks(3) % refs(2,1,1), git % point_edges(:,70402)
!call sw_edge_info(140902)

call write_output(0)


call sw_init_residual
end_adaption = .false.

ITER_LOOP: do while (.not. end_adaption)

   call loop_control()

   call calc_edge_length(max_edge_len,min_edge_len)

   call calc_edge_springs(max_spring,min_spring)

   call calc_edge_forces(max_edge_f)

   call calc_point_forces(max_point_f,sum_point_f)

   call sw_residual(iter,max_spring,min_spring,max_edge_f,max_point_f,sum_point_f,max_edge_len,min_edge_len)

   call move_points(max_edge_f)
   
   call write_output(iter)
end do ITER_LOOP

call unstr2struct(blocks)

call write_grid

call sw_program_end()
end program BSGAS
