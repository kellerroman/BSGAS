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
end module help_routines
