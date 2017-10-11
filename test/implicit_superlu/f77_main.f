*
* Copyright (c) 2003, The Regents of the University of California, through
* Lawrence Berkeley National Laboratory (subject to receipt of any required 
* approvals from U.S. Dept. of Energy) 
* 
* All rights reserved. 
* 
* The source code is distributed under BSD license, see the file License.txt
* at the top-level directory.
*
      program f77_main
         implicit none
      integer maxn, maxnz
      parameter ( maxn = 10000, maxnz = 100000 )
      integer rowind(maxnz), colptr(maxn)
      real*8  values(maxnz), b(maxn)
      integer n, nnz, nrhs, ldb, info, iopt,i 
      integer*8 factors
      real*8 :: l1 = 1,l2 = 1,l3 = 1,l4 = 3
      real*8, parameter :: x0 = 0,x1 = 1,x2 = 2,x3 = 3,x4 = 4

      integer :: iter , niter = 50

      n = 3
      nnz = 7

      rowind(1) = 1
      rowind(2) = 2
      rowind(3) = 1
      rowind(4) = 2
      rowind(5) = 3
      rowind(6) = 2
      rowind(7) = 3

      colptr(1) = 1
      colptr(2) = 3
      colptr(3) = 6
      colptr(4) = 8

      nrhs = 1
      ldb = n

      do iter = 1, niter

      values(1) = -(l2+l1)
      values(2) = l2
      values(3) = l2
      values(4) = -(l2+l3)
      values(5) = l3
      values(6) = l3
      values(7) = -(l4+l3)
      b(1) = -x0 * l1
      b(2) = 0.0d0
      b(3) = -x4 * l4
      !b(:) = 1

*
* First, factorize the matrix. The factors are stored in *factors* handle.
      iopt = 1
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, rowind, colptr, 
     $                      b, ldb, factors, info )
*
      if (info .eq. 0) then
!         write (*,*) 'Factorization succeeded'
      else
         write(*,*) 'INFO from factorization = ', info
      endif
*
* Second, solve the system using the existing factors.
      iopt = 2
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, rowind, colptr, 
     $                      b, ldb, factors, info )
*
      if (info .eq. 0) then
!         write (*,*) 'Solve succeeded'
      else
         write(*,*) 'INFO from triangular solve = ', info
         stop 1
      endif


      !write (*,*) (b(i), i=1, min(n,10))
      write(*,*) b(1),b(2)-b(1),B(3) - B(2), x4 - b(3)

      l4 = l4 * (x4-b(3)) / 0.1d0
      l3 = l3 * (b(3) - b(2)) / (x4 - b(3)) / 1.3d0
      l2 = l2 * (b(2) - b(1)) / (b(3) - b(2)) / 1.3d0

      end do

* Last, free the storage allocated inside SuperLU
      iopt = 3
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, rowind, colptr, 
     $                      b, ldb, factors, info )
*
      stop
      end


