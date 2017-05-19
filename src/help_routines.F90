module help_routines
use const
implicit none
private
public :: alloc
interface alloc
   module procedure alloc_real_1d
   module procedure alloc_real_2d
   module procedure alloc_int_1d
   module procedure alloc_int_2d
   module procedure alloc_int_3d
   module procedure alloc_int_3d_array
   module procedure alloc_log_1d
end interface alloc
contains

subroutine alloc_real_1d(array,d1)
implicit none
real(REAL_KIND), allocatable, intent(out) :: array(:)
integer(INT_KIND), intent(in)             :: d1

allocate (array(d1))

end subroutine alloc_real_1d

subroutine alloc_real_2d(array,d1,d2)
implicit none
real(REAL_KIND), allocatable, intent(out) :: array(:,:)
integer(INT_KIND), intent(in)             :: d1
integer(INT_KIND), intent(in)             :: d2

allocate (array(d1,d2))

end subroutine alloc_real_2d

subroutine alloc_int_1d(array,d1)
implicit none
integer(INT_KIND), allocatable, intent(out) :: array(:)
integer(INT_KIND), intent(in)             :: d1

allocate (array(d1))

end subroutine alloc_int_1d

subroutine alloc_int_2d(array,d1,d2)
implicit none
integer(INT_KIND), allocatable, intent(out) :: array(:,:)
integer(INT_KIND), intent(in)             :: d1
integer(INT_KIND), intent(in)             :: d2

allocate (array(d1,d2))

end subroutine alloc_int_2d

subroutine alloc_int_3d(array,d1,d2,d3)
implicit none
integer(INT_KIND), allocatable, intent(out) :: array(:,:,:)
integer(INT_KIND), intent(in)             :: d1
integer(INT_KIND), intent(in)             :: d2
integer(INT_KIND), intent(in)             :: d3

allocate (array(d1,d2,d3))

end subroutine alloc_int_3d

subroutine alloc_int_3d_array(array,d)
implicit none
integer(INT_KIND), allocatable, intent(out) :: array(:,:,:)
integer(INT_KIND), intent(in)             :: d(3)

call alloc_int_3d(array,d(1),d(2),d(3))

end subroutine alloc_int_3d_array

subroutine alloc_log_1d(array,d1)
implicit none
logical          , allocatable, intent(out) :: array(:)
integer(INT_KIND), intent(in)             :: d1

allocate (array(d1))

end subroutine alloc_log_1d
end module help_routines
