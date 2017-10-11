program SuperLu_1D
implicit none
integer        , allocatable, dimension(:) :: rowind, colptr
real(KIND = 8) , allocatable, dimension(:) :: posX, springs, lengths
real(KIND = 8) , allocatable, dimension(:) :: values, b
integer                                    :: npkts, nnz, nrhs, ldb, info, iopt,i 
integer*8                                  ::  factors

integer :: iter , niter = 50

integer :: pos, c3

npkts = 3
nnz = npkts * 3 - 2

allocate(rowind(nnz))
allocate(colptr(npkts+1))
allocate(values(nnz))
allocate(b(npkts))
allocate(springs(npkts+1))
allocate(lengths(npkts+1))
allocate(posX(npkts+2))
springs = 1.0D0
do i = 1,npkts+2
   posX(i) = dble(i-1)
end do
pos = 0
c3 = 1
do i = 1, nnz
   rowind(i) = pos + c3
   c3 = c3 + 1
   if (c3 == 3) then
      c3 = 0
      pos = pos + 1
   end if
end do

colptr(1) = 1
colptr(2) = 3
do i = 3, npkts
   colptr(i) = colptr(i-1) + 3
end do
colptr(npkts+1) = colptr(npkts) + 2

nrhs = 1
ldb = npkts

do iter = 1, niter

   values(1) = -(springs(2)+springs(1))
   values(2) = springs(2)
   values(3) = springs(2)
   values(4) = -(springs(2)+springs(3))
   values(5) = springs(3)
   values(6) = springs(3)
   values(7) = -(springs(4)+springs(3))
   b(1) = -posX(1) * springs(1)
   b(2) = 0.0d0
   b(3) = -posX(5) * springs(4)
   !b(:) = 1

   !
   ! First, factorize the matrix. The factors are stored in *factors* handle.
   iopt = 1
   call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr,&
                         b, ldb, factors, info )
   !
   if (info .eq. 0) then
   !         write (*,*) 'Factorization succeeded'
   else
      write(*,*) 'INFO from factorization = ', info
   endif
   !
   ! Second, solve the system using the existing factors.
   iopt = 2
   call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr,&
                        b, ldb, factors, info )
   !
   if (info .eq. 0) then
   !         write (*,*) 'Solve succeeded'
   else
      write(*,*) 'INFO from triangular solve = ', info
      stop 1
   endif


   !write (*,*) (b(i), i=1, min(n,10))
   do i = 2,npkts+1
      posX(i) = b(i-1)
   end do
   do i = 1,npkts+1
      lengths(i) = posX(i+1) -posX(i)
   end do
   !write(*,*) posX(:)
   write(*,*) lengths(:)
   i = 1
   springs(i) = springs(i) * lengths(i) / 0.1D0
   do i = 2, npkts 
      springs(i) = springs(i) * lengths(i) / lengths(i-1) / 1.3d0
   end do


end do

! Last, free the storage allocated inside SuperLU
iopt = 3
call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr,&
                     b, ldb, factors, info )
!
stop
end program SuperLU_1d
