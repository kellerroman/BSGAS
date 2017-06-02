!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BLOCK STRUCUTRED GRID ADAPTION SOLVER PROGRAM CONSTANTS MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    this module is intendet to contain all global constants for the program
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

module const
implicit none
public
integer, parameter :: REAL_KIND = selected_real_kind(15 , 307)
integer, parameter :: INT_KIND  = selected_int_kind(8)


enum, bind(C)
   enumerator :: WEST= 1, EAST, SOUTH, NORTH, FRONT, BACK
end enum

end module const
