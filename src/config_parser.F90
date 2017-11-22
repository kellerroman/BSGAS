!==================================================================================================!
!===     Ascii-Config Parser and Binary Config Writer                                           ===!
!===                                                                                            ===!
!===     author: Roman Keller                                                                   ===!
!===                                                                                            ===!
!===     date: 30.01.2017                                                                       ===!
!===                                                                                            ===!
!===                                                                                            ===!
!==================================================================================================!
module config_parser
implicit none
private
integer, parameter         :: CONFIG_FILE_VERSION     =  2 
integer, parameter, public :: VARNAME_LENGTH          = 20
integer, parameter         :: LENGTH_STRING           = 100
integer, parameter, public :: RETURN_CODE_GOOD        =  0
integer, parameter, public :: RETURN_CODE_NOT_IN_LIST = -1
integer, parameter, public :: RETURN_CODE_NOT_A_INT   = -2
integer, parameter, public :: RETURN_CODE_NOT_A_REAL  = -3
integer, parameter, public :: RETURN_CODE_ALREADY_SET = -4
integer, parameter         :: DATATYPE_INT            =  1
integer, parameter         :: DATATYPE_REAL           =  2
integer, parameter         :: DATATYPE_STRING         =  3
type :: t_datin_para
   character(len=VARNAME_LENGTH) :: varname
   logical :: required                          ! i
   logical :: is_set
   integer :: datatype
   integer, pointer :: int_var
   real(kind=8), pointer :: real_var
   character(len=LENGTH_STRING), pointer :: string_var
   type(t_datin_para), pointer :: next => null()
end type t_datin_para
type(t_datin_para), pointer :: first => null()

integer, protected :: n_datin_para = 0


public :: add_parameter                               &
        , list_parameter                              &
        , set_para                                    &
        , get_unset_para                              &
        , get_unset_paras                             &
        , free_all

interface add_parameter
   module procedure add_int
   module procedure add_real
   module procedure add_str
end interface add_parameter
contains
subroutine add_int(str,var_pointer,init_value)
implicit none
character(len=*), intent(in) :: str
integer, target, intent(in) :: var_pointer
integer, intent(in), optional :: init_value

type(t_datin_para), pointer :: new, tmp
! Speicher reservieren
allocate(new)
! Werte setzen
new % varname = lower_case(str)
new % datatype = DATATYPE_INT
new % int_var => var_pointer
new % is_set   = .false.
if ( present(init_value) ) then
   new % int_var = init_value
   new % required = .false.
else
   new % required = .true.
end if
! Am Beginn der Liste einf?gen
if (.not. associated(first)) then
   first => new
else
   tmp => first
   first => new
   first % next => tmp
end if
n_datin_para = n_datin_para + 1
end subroutine add_int

subroutine add_real(str,var_pointer,init_value)
implicit none
character(len=*), intent(in) :: str
real(kind=8), target, intent(in) :: var_pointer
real(kind=8), intent(in), optional :: init_value

type(t_datin_para), pointer :: new, tmp
! Speicher reservieren
allocate(new)
! Werte setzen
new % varname = lower_case(str)
new % datatype = DATATYPE_REAL
new % real_var => var_pointer
new % is_set   = .false.
if ( present(init_value) ) then
   new % real_var = init_value
   new % required = .false.
else
   new % required = .true.
end if
! Am Beginn der Liste einf?gen
if (.not. associated(first)) then
   first => new
else
   tmp => first
   first => new
   first % next => tmp
end if
n_datin_para = n_datin_para + 1
end subroutine add_real

subroutine add_str(str,var_pointer,init_value)
implicit none
character(len=*), intent(in) :: str
character(len=*), target, intent(in) :: var_pointer
character(len=*), intent(in), optional :: init_value

type(t_datin_para), pointer :: new, tmp
! Speicher reservieren
allocate(new)
! Werte setzen
new % varname = lower_case(str)
new % datatype = DATATYPE_STRING
new % string_var => var_pointer
new % is_set   = .false.
!write(*,*) str
if ( present(init_value) ) then
   new % string_var = init_value
   new % required = .false.
else
   new % required = .true.
end if
! Am Beginn der Liste einf?gen
if (.not. associated(first)) then
   first => new
else
   tmp => first
   first => new
   first % next => tmp
end if
n_datin_para = n_datin_para + 1
end subroutine add_str

subroutine list_parameter(io_unit)
use, intrinsic :: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT
implicit none
integer, optional :: io_unit
integer :: npara, io
type(t_datin_para), pointer :: tmp
npara = 0
tmp => first
if (present(io_unit)) then
   io = io_unit
else
   io = stdout
end if
do 
   if (.not. associated(tmp) ) exit
   npara = npara + 1
   if (tmp % is_set .or. .not. tmp % required) then 
      if ( tmp % datatype == DATATYPE_INT) then
         write(io,'(A20," =  ",I0)') tmp % varname, tmp % int_var
      else if ( tmp % datatype == DATATYPE_REAL) then
         write(io,'(A20," = ",ES11.4)') tmp % varname, tmp % real_var
      else if ( tmp % datatype == DATATYPE_STRING) then
         write(io,'(A20," =  ",A)') tmp % varname, trim(tmp % string_var)
      else
         write(stderr,*) "Error in List_parameter: DATAYPE unknown"
         stop 1
      end if
   else
      write(io,'(A20)') tmp % varname
   end if
   tmp =>tmp % next
end do

end subroutine list_parameter

function set_para(varname,varvalue) result (return_value)
implicit none
character(len=VARNAME_LENGTH), intent(in) :: varname
character(len=*), intent(in) :: varvalue
integer :: return_value
type(t_datin_para), pointer :: tmp

return_value = RETURN_CODE_GOOD
tmp => first
do 
   if (.not. associated(tmp) ) then
      return_value = RETURN_CODE_NOT_IN_LIST
      exit
   end if

   if (tmp % varname == lower_case(varname)) then
      if (tmp % is_set) then
         return_value = RETURN_CODE_ALREADY_SET
      end if
      if (tmp % datatype == DATATYPE_INT) then
         if (is_integer(varvalue)) then
            !write(*,*) "Changing ",trim(varname)," to value ",varvalue
            read(varvalue,*) tmp % int_var
            tmp % is_set = .true.
         else
            return_value = RETURN_CODE_NOT_A_INT
            exit
         end if
      else if (tmp % datatype == DATATYPE_REAL) then
         if (is_real(varvalue)) then
            !write(*,*) "Changing ",trim(varname)," to value ",varvalue
            read(varvalue,*) tmp % real_var
            tmp % is_set = .true.
         else
            return_value = RETURN_CODE_NOT_A_REAL
            exit
         end if
      else if (tmp % datatype == DATATYPE_STRING) then
            !write(*,*) "Changing ",trim(varname)," to value ",varvalue
            tmp % string_var = varvalue
            tmp % is_set = .true.
      else 
         write(*,*) "Error in SET_PARA",__FILE__,__LINE__
         stop 1
      end if
      exit
   end if
   
   tmp =>tmp % next
end do


end function set_para

function get_unset_para(last_para) result (para_name)
implicit none
character(len=VARNAME_LENGTH), intent(in) :: last_para
character(len=VARNAME_LENGTH) :: para_name
type(t_datin_para), pointer :: tmp
tmp => first
para_name = ""
do 
   if (.not. associated(tmp) ) exit
   if (    tmp % required              .and. &
      .not.tmp % is_set                .and. &
           tmp % varname /= last_para) then
      para_name = tmp % varname
      exit
   end if
   tmp =>tmp % next
end do

end function get_unset_para

subroutine get_unset_paras(para_name)
implicit none
character(len=VARNAME_LENGTH),allocatable, intent(out) :: para_name(:)
type(t_datin_para), pointer :: tmp
integer :: nunset
tmp => first
nunset = 0
do 
   if (.not. associated(tmp) ) exit
   if (    tmp % required              .and. &
      .not.tmp % is_set             ) then
      nunset = nunset + 1
   end if
   tmp =>tmp % next
end do
allocate(para_name(nunset))
nunset = 0
tmp => first
do 
   if (.not. associated(tmp) ) exit
   if (    tmp % required              .and. &
      .not.tmp % is_set             ) then
      nunset = nunset + 1
      para_name(nunset) = tmp % varname
   end if
   tmp =>tmp % next
end do
end subroutine get_unset_paras
subroutine free_all()
implicit none
type(t_datin_para), pointer :: tmp
do        
   tmp => first
   if (.not. associated(tmp)) exit
   first => first%next
   deallocate(tmp)
   end do                     
end subroutine free_all

function lower_case( input_string ) result ( output_string )
   ! -- argument and result
   implicit none
   character( * ), intent( in )       :: input_string
   character( len( input_string ) )   :: output_string
   integer                            :: ii,ic,nlen
   nlen = len(input_string)
   do ii=1,nlen
      ic = ichar(input_string(ii:ii))
      if (ic >= 65 .and. ic <= 90) then
         output_string(ii:ii) = char(ic+32)
      else
         output_string(ii:ii) = char(ic)
      end if
   end do
end function lower_case
logical function is_integer(string)
implicit none
integer :: ipos
character (len = *) :: string
is_integer =  .true.
do ipos = 1,len_trim(string)
   if ((ichar(string(ipos:ipos)) >57 .or. ichar(string(ipos:ipos)) <48 )&
      .and.string(ipos:ipos) /= "-" .and.string(ipos:ipos) /= "+") then
!            write(*,*) trim(string) , "ist keine integer-zahl"
      is_integer =  .false.
      return
   end if
end do

   end function

logical function is_real(string)
implicit none
integer :: ipos
integer :: elem_pos
character (len = *) :: string
is_real =  .true.
elem_pos = 1
!!!! aufteilung eines realen zahl   - 0 . 0 d + 1
! vorzeichen                  elem_pos = 1
! zahl vor dem komma          elem_pos = 2
! trennzeichen                elem_pos = 3
! zahl nach dem komma         elem_pos = 4
! exponent                    elem_pos = 5
! exponenten-vorzeichen       elem_pos = 6
! exponent                    elem_pos = 7
do ipos = 1,len_trim(string)
   if (ichar(string(ipos:ipos)) <=57 .and. ichar(string(ipos:ipos)) >=48) then
      select case(elem_pos)
         case (2)
         case (4)
         case (5)
            elem_pos = 7
         case (7)
         case default
            elem_pos = elem_pos + 1
      end select
   else if (string(ipos:ipos) == "-" .or. string(ipos:ipos) == "+") then
      select case(elem_pos)
         case (2:4)
            is_real = .false.
            return
         case (6:7)
            is_real = .false.
            return
         case default
            elem_pos = elem_pos + 1
      end select
   else if (string(ipos:ipos) == ".") then
      select case(elem_pos)
         case (1:2)
            elem_pos = elem_pos + 1
         case default
            is_real = .false.
            return
      end select
   else if (string(ipos:ipos) == "d" .or. string(ipos:ipos) == "e" .or.  &
            string(ipos:ipos) == "D" .or. string(ipos:ipos) == "E") then
      select case(elem_pos)
         case (2:4)
            elem_pos = 5
         case default
            is_real = .false.
            return
      end select
   end if
end do
if (elem_pos <= 2 .or. &
    elem_pos == 5 .or. &
    elem_pos == 6 ) then
   is_real = .false.
   return
end if
end function
end module config_parser
