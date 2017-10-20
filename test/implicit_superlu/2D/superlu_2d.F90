program SuperLu_2D
implicit none
real(kind = 8) , parameter                   :: LENGTH      = (5.080d-02+6.350d-03)/2.0D0
real(kind = 8) , parameter                   :: WALL_LENGTH = 1.0d-6
real(kind = 8) , parameter                   :: CELL_INC    = 1.2d0
real(kind = 8) , parameter                   :: SPRING_MIN  = 1.0d-2

integer        , allocatable, dimension(:)   :: rowind, colptr        ! Matrix-Entries Position arrays
integer        , allocatable, dimension(:)   :: mat_dia               ! Indices of Diagonal Matrix elements
real(kind = 8) , allocatable, dimension(:)   :: values                ! Array-Entries
real(kind = 8) , allocatable, dimension(:)   :: rhs                   ! RHS-Entries
integer        , allocatable, dimension(:,:) :: mat_to_id             ! Matrix-Entries Position arrays

real(kind = 8) , allocatable, dimension(:)   :: springs, lengths, res

real(kind = 8) , allocatable, dimension(:,:,:) :: Pkt2D_coords
real(kind = 8) , allocatable, dimension(:,:) :: Pkt_id_coords

integer        , allocatable, dimension(:,:) :: spring_points         ! Pointer auf die 1D Nummer der Spring Punkte

integer        , allocatable, dimension(:,:) :: spring_pos_mat        ! Position der Springeintrage im Matrix array
integer        , allocatable, dimension(:)   :: spring_nmat_entry     ! Number of Entries of the spring in Matrix
integer        , allocatable, dimension(:,:) :: spring_pos_rhs        ! Position der Springeintrage im RHS array
real(kind = 8) , allocatable, dimension(:,:) :: spring_cell_rhs       ! Position der Springeintrage im RHS array
real(kind = 8) , allocatable, dimension(:,:) :: spring_dim_rhs        ! Position der Springeintrage im RHS array
integer        , allocatable, dimension(:)   :: spring_nrhs_entry     ! Number of Entries of the spring in RHS

integer        , allocatable, dimension(:)   :: spring_nneigh         ! Number of Entries of the spring in Matrix
integer        , allocatable, dimension(:,:) :: spring_neigh          ! Position der Springeintrage im RHS array

integer        , allocatable, dimension(:,:) :: Pkt2D_to_id           ! id des i,j- Punktes
integer        , allocatable, dimension(:,:) :: Pkt_id_to_ijk              ! 2D Position anhand der Pkt_id 
integer        , allocatable, dimension(:,:) :: point_pos_mat

integer        , allocatable, dimension(:,:) :: Pkt_springs           ! Ids of springs connected to this point
integer        , allocatable, dimension(:,:) :: Pkt_connected_Pkt      ! Ids of Pkt connected to this point
integer        , allocatable, dimension(:)   :: Pkt_nspring           ! Number of springs connected to this point
integer                                      :: cell_i, cell_j
integer                                      :: i, j
integer                                      :: npkts, nnz, nrhs, ldb, info, iopt
integer                                      :: nspring, ngrid_pkts
integer*8                                    ::  factors

integer                                      :: iter
integer        , parameter                   :: niter = 5000

integer                                      :: pos, c3, eq, found, found_rhs, d, ind
integer                                      :: p_eq, p_pos

integer                                      :: delta_index

logical                                      :: print_mat

real(kind = 8)                               :: temp
cell_i = 3
cell_j = 3 
 print_mat = .true.


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
!
!
!#MPKTS   001  002  003  004  005  006  007  008  009  010  011  012  013  014  015  016 
!Points  x002 x003 y005 x006 y006 x007 y007 y008 y009 x010 y010 x011 y011 y012 x014 x015 
!Index i    2    3    1    2    2    3    3    4    1    2    2    3    3    4    2    3 
!Index j    1    1    2    2    2    2    2    2    3    3    3    3    3    3    4    4

!        l1+l2 l2        l6                                                               l1*x1
!        +l6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         l2  l2+l8               l8                                                      l3*x4
!             +l3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  l4+l5     l5                  l11                                      l4*y1
!                  +l11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         l6           l5+l6      l7                   l13                                l5*x5
!                     +l7+l13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 l5       l5+l6      l7                   l13                            l6*y2
!                         +l7+l13

ngrid_pkts = (cell_i + 1) * (cell_j + 1)
npkts    = 2 * ngrid_pkts - 2 * (cell_i + 1) - 2 * (cell_j + 1)
nspring  = (cell_i + 1) * (cell_j) + (cell_i) * (cell_j + 1)
nnz      = npkts * 5 - 2 * (cell_i + 1) - 2 * (cell_j + 1)

write(*,*) "Points:",npkts,"Springs:",nspring,"NonZero Elem.:",nnz

allocate(res(npkts))

allocate(springs(nspring))
allocate(spring_points(2,nspring))
allocate(lengths(nspring))

allocate(spring_pos_mat(8,nspring))
allocate(spring_nmat_entry(nspring))

allocate(spring_pos_rhs(8,nspring))
allocate(spring_cell_rhs(8,nspring))
allocate(spring_dim_rhs(8,nspring))
allocate(spring_nrhs_entry(nspring))

allocate (spring_nneigh(nspring))
allocate (spring_neigh(2,nspring))

allocate(Pkt2D_coords (2,cell_i+1,cell_j+1))
allocate(Pkt2D_to_id(cell_i+1,cell_j+1))
allocate(Pkt_id_coords(2,(cell_i+1)*(cell_j+1)))
allocate(Pkt_id_to_ijk(2,ngrid_pkts))
allocate(Pkt_springs  (4,ngrid_pkts))
allocate(Pkt_connected_Pkt  (4,ngrid_pkts))
allocate(Pkt_nspring  (ngrid_pkts))

allocate(point_pos_mat(2,ngrid_pkts))

allocate(mat_to_id(2,npkts))

springs           = 1.0D0
point_pos_mat     = -1

spring_nneigh     = 0

spring_nmat_entry = 0
spring_nrhs_entry = 0

Pkt_nspring       = 0

c3 = 0
pos = 0
do j = 1, cell_j+1
   do i = 1, cell_i+1
      !Pkt2D_coords(1,i,j) = length * dble(i-1) / dble(npkts+1)
      Pkt2D_coords(1,i,j) = dble(i-1)/dble(cell_i)
      Pkt2D_coords(2,i,j) = dble(j-1)/dble(cell_j)
      pos = pos + 1
      Pkt2D_to_id(i,j) = pos
      Pkt_id_coords(:,pos) = Pkt2D_coords(:,i,j)
      Pkt_id_to_ijk(1,pos) = i
      Pkt_id_to_ijk(2,pos) = j
      if (i > 1 .and. i <= cell_i) then
         c3 = c3 + 1
         point_pos_mat(1,pos) = c3
         mat_to_id(1,c3) = 1
         mat_to_id(2,c3) = pos
      end if
      if (j > 1 .and. j <= cell_j) then
         c3 = c3 + 1
         point_pos_mat(2,pos) = c3
         mat_to_id(1,c3) = 2
         mat_to_id(2,c3) = pos
      end if
   end do
end do

if (npkts /= c3) then
   write(*,*) "NPKT wrongly approimated",npkts,c3
   stop 1
end if

if (print_mat) then
   do i = 1, ngrid_pkts
      if (point_pos_mat(1,i) /= -1) then
         write(*,'(" ",I3.3,1X)',ADVANCE="NO") point_pos_mat(1,i)
      end if
      if (point_pos_mat(2,i) /= -1) then
         write(*,'(" ",I3.3,1X)',ADVANCE="NO") point_pos_mat(2,i)
      end if
   end do
   write(*,*)

   do i = 1, ngrid_pkts
      if (point_pos_mat(1,i) /= -1) then
         write(*,'("x",I3.3,1X)',ADVANCE="NO") i
      end if
      if (point_pos_mat(2,i) /= -1) then
         write(*,'("y",I3.3,1X)',ADVANCE="NO") i
      end if
   end do
   write(*,*)

   do i = 1, ngrid_pkts
      if (point_pos_mat(1,i) /= -1) then
         write(*,'(" ",I3,1X)',ADVANCE="NO") Pkt_id_to_ijk(1,i)
      end if
      if (point_pos_mat(2,i) /= -1) then
         write(*,'(" ",I3,1X)',ADVANCE="NO") Pkt_id_to_ijk(1,i)
      end if
   end do
   write(*,*)
   do i = 1, ngrid_pkts
      if (point_pos_mat(1,i) /= -1) then
         write(*,'(" ",I3,1X)',ADVANCE="NO") Pkt_id_to_ijk(2,i)
      end if
      if (point_pos_mat(2,i) /= -1) then
         write(*,'(" ",I3,1X)',ADVANCE="NO") Pkt_id_to_ijk(2,i)
      end if
   end do
   write(*,*)
end if


nrhs = 1
ldb = npkts

pos = 0
do j = 1, cell_j + 1
   do i = 1, cell_i + 1
      if (i > 1) then
         pos = pos + 1
         spring_points(1,pos) = Pkt2D_to_id(i-1,j)
         spring_points(2,pos) = Pkt2D_to_id(i,j)
         do d = 1,2
            c3 = spring_points(d  ,pos)
            eq = spring_points(3-d,pos)
            ind = Pkt_nspring(c3) + 1
            Pkt_nspring(c3) = ind
            Pkt_springs(ind,c3) = pos
            Pkt_connected_Pkt(ind,c3) = eq
         end do
      end if
      if (i > 2) then
         if (j > 1 .and. j < cell_j+1) then
            c3 = 2
         else
            c3 = 1
         end if

         ind = spring_nneigh(pos) + 1
         spring_nneigh(pos) = ind
         spring_neigh(ind,pos) = pos - c3

         ind = spring_nneigh(pos-c3) + 1
         spring_nneigh(pos-c3) = ind
         spring_neigh(ind,pos - c3) = pos
      end if
      if (j > 1) then
         pos = pos + 1
         spring_points(1,pos) = Pkt2D_to_id(i,j-1)
         spring_points(2,pos) = Pkt2D_to_id(i,j)
         do d = 1,2
            c3 = spring_points(d,pos)
            eq = spring_points(3-d,pos)
            ind = Pkt_nspring(c3) + 1
            Pkt_nspring(c3) = ind
            Pkt_springs(ind,c3) = pos
            Pkt_connected_Pkt(ind,c3) = eq
         end do
      end if
   end do
end do

if (pos /= nspring) then
   write(*,*) "NSPRING wrongly approimated",nspring,pos
   stop 1
end if

delta_index = 0
do pos = 1, nspring
   delta_index = max(delta_index, abs(spring_points(2,pos) - spring_points(1,pos)))
end do
write(*,*) "Maxiumum Difference", delta_index
delta_index = delta_index * 2

iter = 0
do eq = 1, npkts
   d = mat_to_id(1,eq)
   p_eq = mat_to_id(2,eq)
   do pos = 1, npkts
   !do pos = max(1,eq-delta_index), min(npkts,eq+delta_index)
      found = 0
      if (d == mat_to_id(1,pos)) then
         p_pos = mat_to_id(2,pos)
         if (pos == eq) then 
            do c3 = 1, Pkt_nspring(p_eq)
               i = Pkt_springs(c3,p_eq)
               found = found + i
            end do
         else 
            do c3 = 1, Pkt_nspring(p_eq)
               if (p_pos == Pkt_connected_Pkt(c3,p_eq)) then
                  found = Pkt_springs(c3,p_eq)
                  exit
               end if
            end do
         end if
      end if
      if (found /= 0) then
         if (print_mat) &
            write(*,'(" ",I3.3,1X)',ADVANCE="NO") found
         iter = iter + 1
      else  
         if (print_mat) &
            write(*,'("    ",1X)',ADVANCE="NO")
      end if
   end do
   do c3 = 1, Pkt_nspring(p_eq)
      if (point_pos_mat(d,Pkt_connected_Pkt(c3,p_eq)) == -1) then
         if (print_mat) &
            write(*,'(" ",I3.3,1X)',ADVANCE="NO") Pkt_springs(c3,p_eq)
      end if
   end do
   if (print_mat) &
      write(*,'("    ",1X)')
end do

if (nnz /= iter) then
   write(*,*) "Number of Non-Zero Elemnents Wrongly Approximated",nnz,iter
   !stop 1
   nnz = iter
end if
stop
allocate(rowind(nnz))
allocate(colptr(npkts+1))

allocate(values(nnz))
allocate(rhs(npkts))
allocate(mat_dia(npkts))

iter = 0
do pos = 1, npkts
   colptr(pos) = iter + 1
   p_pos = mat_to_id(2,pos)
   d     = mat_to_id(1,pos)
   do eq = max(1,pos-delta_index), min(npkts,pos+delta_index)
      p_eq = mat_to_id(2,eq)
      if (d == mat_to_id(1,eq)) then
         if (pos == eq) then 
            iter = iter + 1
            rowind(iter) = eq
            mat_dia(eq) = iter
            do c3 = 1, Pkt_nspring(p_eq)
               i = Pkt_springs(c3,p_eq)
               ind = spring_nmat_entry(i) + 1
               spring_nmat_entry(i) = ind
               spring_pos_mat(ind,i) = iter
            end do
         else 
            do c3 = 1, Pkt_nspring(p_eq)
               if (p_pos == Pkt_connected_Pkt(c3,p_eq)) then
                  iter = iter + 1
                  rowind(iter) = eq
                  i = Pkt_springs(c3,p_eq)
                  ind = spring_nmat_entry(i) + 1
                  spring_nmat_entry(i) = ind
                  spring_pos_mat(ind,i) = iter
                  exit
               end if
            end do
         end if
      end if
   end do
end do
do eq = 1, npkts
   d = mat_to_id(1,eq)
   p_eq = mat_to_id(2,eq)
   do c3 = 1, Pkt_nspring(p_eq)
      if (point_pos_mat(d,Pkt_connected_Pkt(c3,p_eq)) == -1) then
         i = Pkt_springs(c3,p_eq)
         ind = spring_nrhs_entry(i) + 1
         spring_nrhs_entry(i) = ind
         spring_pos_rhs(ind,i) = eq
         spring_cell_rhs(ind,i) = Pkt_connected_Pkt(c3,p_eq)
         spring_dim_rhs(ind,i) = d
      end if
   end do
end do
colptr(npkts+1) = iter + 1
!write(*,*) rowind
!write(*,*) colptr
!write(*,*) mat_dia
!do i = 1, nspring
!   write(*,*) i, spring_nmat_entry(i), spring_nrhs_entry(i) &
!               , spring_pos_mat(1:spring_nmat_entry(i),i)   &
!               , spring_pos_rhs(1:spring_nrhs_entry(i),i)
!end do


do i = 1,nspring
   lengths(i) = sqrt ( (Pkt_id_coords(1,spring_points(1,i)) - Pkt_id_coords(1,spring_points(2,i))) &
                     * (Pkt_id_coords(1,spring_points(1,i)) - Pkt_id_coords(1,spring_points(2,i))) &
                     + (Pkt_id_coords(2,spring_points(1,i)) - Pkt_id_coords(2,spring_points(2,i))) &
                     * (Pkt_id_coords(2,spring_points(1,i)) - Pkt_id_coords(2,spring_points(2,i))) )
end do

i = 1
springs(i) = springs(i) * lengths(i) / wall_length
do i = 1, nspring
   if (spring_nneigh(i) > 0) then
      temp = 1.0D10
      do ind = 1, spring_nneigh(i)
         temp = min(temp,lengths(spring_neigh(ind,i)))
      end do
      springs(i) = springs(i) * lengths(i) / temp / cell_inc
   end if
end do
!do i = 2, npkts !+ 1
!   if (i > npkts) then
!      springs(i) = springs(i) * lengths(i) / lengths(i-1) / cell_inc
!   else if (i == 1) then
!      springs(i) = springs(i) * lengths(i) / lengths(i+1) / cell_inc
!   else
!      springs(i) = springs(i) * lengths(i) / min(lengths(i+1),lengths(i-1)) / cell_inc
!   end if
!   springs(i) = max(springs(i), spring_min)
!end do
!i = npkts + 1
!springs(i) = springs(i) * lengths(i) / wall_length
iter = 0
do 
   iter = iter + 1

   values = 0
   rhs    = 0
   do i = 1, nspring
      do c3 = 1, spring_nmat_entry(i)
         pos = spring_pos_mat(c3,i)
         values(pos) = values(pos) + springs(i)
      end do
      do c3 = 1, spring_nrhs_entry(i)
         pos = spring_pos_rhs(c3,i)
         ind = spring_cell_rhs(c3,i)
           d = spring_dim_rhs(c3,i)
         rhs(pos) = rhs(pos) - springs(i) * Pkt_id_coords(d,ind)
      end do
   end do
   do i = 1, npkts
      pos = mat_dia(i)
      values(pos) = -values(pos)
   end do

   ! First, factorize the matrix. The factors are stored in *factors* handle.
   iopt = 1
   call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr, rhs, ldb, factors, info )
   !
   if (info .eq. 0) then
   !         write (*,*) 'Factorization succeeded'
   else
      write(*,*) 'INFO from factorization = ', info
   endif
   !
   ! Second, solve the system using the existing factors.
   iopt = 2
   call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr, rhs, ldb, factors, info )
   !
   if (info .eq. 0) then
   !         write (*,*) 'Solve succeeded'
   else
      write(*,*) 'INFO from triangular solve = ', info
      stop 1
   endif
   
   do i = 1,npkts
      res(i) = abs(Pkt_id_coords(mat_to_id(1,i),mat_to_id(2,i)) - rhs(i))
      Pkt_id_coords(mat_to_id(1,i),mat_to_id(2,i)) = rhs(i)
   end do

   do i = 1,nspring
      lengths(i) = sqrt ( (Pkt_id_coords(1,spring_points(1,i)) - Pkt_id_coords(1,spring_points(2,i))) &
                        * (Pkt_id_coords(1,spring_points(1,i)) - Pkt_id_coords(1,spring_points(2,i))) &
                        + (Pkt_id_coords(2,spring_points(1,i)) - Pkt_id_coords(2,spring_points(2,i))) &
                        * (Pkt_id_coords(2,spring_points(1,i)) - Pkt_id_coords(2,spring_points(2,i))) )
   end do

   write(*,*) iter, maxval(res), sum(res), lengths(1:5), springs(1:5)

   i = 1
   springs(i) = springs(i) * lengths(i) / wall_length
   do i = 1, nspring
      if (spring_nneigh(i) > 0) then
         temp = 1.0D10
         do ind = 1, spring_nneigh(i)
            temp = min(temp,lengths(spring_neigh(ind,i)))
         end do
         if (lengths(i) / temp > cell_inc) then
            springs(i) = springs(i) * lengths(i) / temp / cell_inc
         end if
      end if
   end do

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
call c_fortran_dgssv( iopt, npkts, nnz, nrhs, values, rowind, colptr, rhs, ldb, factors, info )
!
end program SuperLU_2D
