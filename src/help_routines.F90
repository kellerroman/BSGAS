module help_routines
use const
implicit none
private
real(REAL_KIND), parameter :: EPS = 1.0E-9_REAL_KIND
public :: alloc, vec_common, vec_normalize, lower_case,cross_product, vec_same,scalar_product

interface alloc
   module procedure alloc_real_1d
   module procedure alloc_real_2d
   module procedure alloc_real_3d
   module procedure alloc_real_3d_array
   module procedure alloc_int_1d
   module procedure alloc_int_2d
   module procedure alloc_int_3d
   module procedure alloc_int_3d_array
   module procedure alloc_log_1d
   module procedure alloc_log_2d
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

subroutine alloc_real_3d(array,d1,d2,d3)
implicit none
real(REAL_KIND), allocatable, intent(out) :: array(:,:,:)
integer(INT_KIND), intent(in)             :: d1
integer(INT_KIND), intent(in)             :: d2
integer(INT_KIND), intent(in)             :: d3

allocate (array(d1,d2,d3))

end subroutine alloc_real_3d

subroutine alloc_real_3d_array(array,d)
implicit none
real(REAL_KIND), allocatable, intent(out) :: array(:,:,:)
integer(INT_KIND), intent(in)             :: d(3)

call alloc_real_3d(array,d(1),d(2),d(3))

end subroutine alloc_real_3d_array

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

subroutine alloc_log_2d(array,d1,d2)
implicit none
logical          , allocatable, intent(out) :: array(:,:)
integer(INT_KIND), intent(in)             :: d1, d2

allocate (array(d1,d2))

end subroutine alloc_log_2d

subroutine vec_normalize(vec)
implicit none
real(REAL_KIND), intent(inout) :: vec(3)

real(REAL_KIND) :: l

l = sqrt(vec(1) * vec(1) + vec(2) * vec(2) + vec(3) * vec(3) )
vec = vec / l

end subroutine vec_normalize

subroutine vec_common(v1,v2)
implicit none
real(REAL_KIND), intent(inout) :: v1(3)
real(REAL_KIND), intent(in)    :: v2(3)
real(REAL_KIND) :: sp
real(REAL_KIND) :: tv(3)

call vec_normalize(v1)
tv = v2
call vec_normalize(tv)

sp = scalar_product(v1,tv)

if (abs(abs(sp)-1.0D0) > EPS) then
   v1 = 0
end if

end subroutine vec_common

real(REAL_KIND) function scalar_product(v1,v2)
implicit none
real(REAL_KIND), intent(in) :: v1(3) , v2(3)

scalar_product = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

end function scalar_product

subroutine cross_product(v1,v2,v3)
implicit none
real(REAL_KIND), dimension(3), intent(out) :: v3
real(REAL_KIND), dimension(3), intent(in)  :: v1,v2
v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
call vec_normalize(v3)
end subroutine

function vec_same(v1,v2) result (are_same)
implicit none
real(REAL_KIND), intent(in)      :: v1(3)
real(REAL_KIND), intent(in)      :: v2(3)
logical                          :: are_same
if (abs(v1(1) - v2(1)) < EPS .and. &
    abs(v1(2) - v2(2)) < EPS .and. &
    abs(v1(3) - v2(3)) < EPS ) then
   are_same = .true.
else
   are_same = .false.
end if
end function vec_same

function lower_case( input_string ) result ( output_string )
   implicit none
   character( * ), intent( in )       :: input_string
   character( len( input_string ) )   :: output_string
   integer                            :: ii,ic,nlen
   nlen = len(input_string)
   do ii=1,nlen
      ic = ichar(input_string(ii:ii))
      if (ic >= 65 .and. ic < 90) then
         output_string(ii:ii) = char(ic+32)
      else
         output_string(ii:ii) = char(ic)
      end if
   end do
end function lower_case

end module help_routines
