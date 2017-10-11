program SuperLu_1D
implicit none
integer        , allocatable, dimension(:)   :: rowind, colptr
real(KIND = 8) , allocatable, dimension(:)   :: springs, lengths, res
real(KIND = 8) , allocatable, dimension(:)   :: values
real(KIND = 8) , allocatable, dimension(:,:,:) :: coord
integer        , allocatable, dimension(:,:) :: spring_pos,spring_points, point_index
real(kind = 8) , parameter                   :: length      = (5.080d-02+6.350d-03)/2.0D0
real(kind = 8) , parameter                   :: wall_length = 1.0d-6
real(kind = 8) , parameter                   :: cell_inc    = 1.2d0
real(kind = 8) , parameter                   :: spring_min  = 1.0d-2
integer                                      :: cell_i, cell_j
integer                                      :: i, j
integer                                      :: npkts, nnz, nrhs, ldb, info, iopt
integer                                      :: nspring
integer*8                                    ::  factors

integer :: iter , niter = 5000

integer :: pos, c3

cell_i = 3 
cell_j = 3 


!  [13]---(19) --- [14]---(21) ---[15] ---(23) --- [16]
!   |               |               |               |
!  (18)            (20)            (22)            (24)
!   |               |               |               |
!  [9] --- (12)---[10] --- (14)---[11] ---(16) --- [12]
!   |               |               |               |
!  (11)            (13)            (15)            (17)
!   |               |               |               |
!  [5] --- (5) --- [6] --- (7) --- [7] --- (9) --- [8]
!   |               |               |               |
!  (4)             (6)             (8)             (10)
!   |               |               |               |
!  [1] --- (1) --- [2] --- (2) --- [3] --- (3) --- [4]


npkts    = 2 * (cell_i + 1) * (cell_j + 1) - 2 * (cell_i + 1) - 2 * (cell_j + 1)
nspring  = (cell_i + 1) * (cell_j) + (cell_i) * (cell_j + 1)
nnz      = npkts * 3 - 2 * (cell_i + 1) - 2 * (cell_j + 1)

write(*,*) "Points:",npkts,"Springs:",nspring,"NonZero Elem.:",nnz

allocate(rowind(nnz))
allocate(colptr(npkts+1))

allocate(values(nnz+npkts))
allocate(res(npkts))

allocate(springs(nspring))
allocate(spring_pos(4,nspring))
allocate(spring_points(2,nspring))
allocate(lengths(npkts+1))
allocate(coord(2,cell_i+1,cell_j+1))
allocate(point_index(cell_i+1,cell_j+1))

springs = 1.0D0

pos = 0
do j = 1, cell_j+1
   do i = 1, cell_i+1
      !coord(1,i,j) = length * dble(i-1) / dble(npkts+1)
      coord(1,i,j) = dble(i-1)
      coord(2,i,j) = dble(j-1)
      pos = pos + 1
      point_index(i,j) = pos
   end do
end do

nrhs = 1
ldb = npkts

pos = 0
do j = 1, cell_j + 1
   do i = 1, cell_i + 1
      if (i > 1) then
         pos = pos + 1
         spring_points(1,pos) = point_index(i-1,j)
         spring_points(2,pos) = point_index(i,j)
      end if
      if (j > 1) then
         pos = pos + 1
         spring_points(1,pos) = point_index(i,j-1)
         spring_points(2,pos) = point_index(i,j)
      end if
   end do
end do

if (pos /= nspring) then
   write(*,*) "NSPRING wrongly approimated",nspring,pos
   stop 1
end if

do i = 1, nspring
   
end do

stop 

do 
   iter = iter + 1
   pos = 0
   c3 = 1
   do i = 1, nnz
      if (c3 == 1) then
         values(i) = - ( springs(pos+c3) + springs(pos+c3+1))
      else if (c3 == 0) then
         values(i) = springs(pos+c3+1)
      else
         values(i) = springs(pos+c3)
      end if
      c3 = c3 + 1
      if (c3 == 3) then
         c3 = 0
         pos = pos + 1
      end if
   end do

!   values(nnz+1) = -posX(1) * springs(1)
!   values(nnz+2:nnz+npkts-1) = 0.0d0
!   values(nnz+npkts) = -posX(npkts+2) * springs(npkts+1)

   ! First, factorize the matrix. The factors are stored in *factors* handle.
   iopt = 1
   call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr, values(nnz+1), ldb, factors, info )
   !
   if (info .eq. 0) then
   !         write (*,*) 'Factorization succeeded'
   else
      write(*,*) 'INFO from factorization = ', info
   endif
   !
   ! Second, solve the system using the existing factors.
   iopt = 2
   call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr, values(nnz+1), ldb, factors, info )
   !
   if (info .eq. 0) then
   !         write (*,*) 'Solve succeeded'
   else
      write(*,*) 'INFO from triangular solve = ', info
      stop 1
   endif


   
!   do i = 2,npkts+1
!      res(i-1) = abs(posX(i) - values(nnz+i-1))
!      posX(i) = values(nnz+i-1)
!   end do
!   do i = 1,npkts+1
!      lengths(i) = posX(i+1) -posX(i)
!   end do
!   !write(*,*) posX(:)
!   write(*,*) lengths(:)
   write(*,*) iter, maxval(res), sum(res)
   i = 1
   springs(i) = springs(i) * lengths(i) / wall_length
   do i = 2, npkts !+ 1
      if (i > npkts) then
         springs(i) = springs(i) * lengths(i) / lengths(i-1) / cell_inc
      else if (i == 1) then
         springs(i) = springs(i) * lengths(i) / lengths(i+1) / cell_inc
      else
         springs(i) = springs(i) * lengths(i) / min(lengths(i+1),lengths(i-1)) / cell_inc
      end if
      springs(i) = max(springs(i), spring_min)
   end do
   i = npkts + 1
   springs(i) = springs(i) * lengths(i) / wall_length

   if (maxval(res) < 1D-14) then
      write(*,*) "Convergent Solution Reached" 
      exit
   end if

   if (iter >= niter) then
      write(*,*) "Maximum Number of Iterations reached"
      stop 1
   end if
end do

!do i = 1,npkts+1
!   write(123,*) i,lengths(i), springs(i)
!end do
!
!do i = 1,npkts+2
!   write(1234,*) i,posX(i)
!end do

! Last, free the storage allocated inside SuperLU
iopt = 3
call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr, values(nnz+1), ldb, factors, info )
!
end program SuperLU_1d
