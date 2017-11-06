program gridgen_o_grid
   use mod_gridgen
implicit none
integer, parameter :: NB = 23
integer, parameter :: nVar   = 0
real(kind = 8), parameter :: Pi = 3.1415927/180.D0

integer :: nk  = 32
integer :: ny  = 32
integer :: ny2 = 32
integer :: nj_wall = 32!64

real(kind = 8),parameter :: radius      = 3.0D-3
real(kind = 8),parameter :: radius_o    = radius * 0.3d0
real(kind = 8),parameter :: radius_out1 =13.0D-3
real(kind = 8),parameter :: radius_wall =15.0D-3

real(kind = 8) :: p1_x = 0.0D-3
real(kind = 8) :: p1_y = 0.0D-3

real(kind = 8) :: p2_x = 9.0D-3
real(kind = 8) :: p2_y = 0.0D-3

real(kind = 8) :: p3_x = 3.0D0 * radius * sin(pi * 30)
real(kind = 8) :: p3_y = 3.0D0 * radius * cos(pi * 30)

integer :: ni, nj

integer :: i, j, k, b
real(kind = 8) :: spktx,spkty,epkty,epktx,lx,ly, winkel,w2
!character(len = 1) :: arg
!integer :: testcase
WRITE(*,*) "GRIDGEN FÜR EIN EINFACHES O_GRID"

!testcase = 1
!call get_command_argument(1, arg)
!if (len_trim(arg) /= 0) then
!   read(arg,*) testcase
!end if
!if (testcase == 1) then
!   nk = 16
!   ny = 32
!else
!   nk = 64
!   ny = 128
!end if
call set_dimension(2)

call add_block(nk / 2, nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny2 )
call add_block(nk / 2, ny2 )

call add_block(nk / 2, nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)

call add_block(nk / 2, nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny - nk / 2)
call add_block(nk / 2, ny2)
call add_block(nk / 2, ny2)
call add_block(nk / 2, ny2)
call add_block(nk / 2, ny2)

call add_block(nk / 2, nj_wall)
call add_block(nk / 2, nj_wall)

call allocate_blocks(nVar)

write(*,'(2(I4,1X))') ((blocks(i) % nCells(j),J=1,2),I=1,NB)

b = 1
ni = blocks(b) % nPkts(1)
nj = blocks(b) % nPkts(2)
 k = 1
do j = 1, nj
   do i = 1,ni
      blocks(b) % xyzs(i,j,k,1) = dble(i-1)*RADIUS*0.3d0/dble(ni-1) + p1_x
      blocks(b) % xyzs(i,j,k,2) = dble(j-1)*RADIUS*0.3d0/dble(nj-1) + p1_y
   end do
end do

b = 2
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (45.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      spkty = dble(ni-i) * radius * 0.3d0 / dble(ni-1)
      spktx = radius * 0.3d0
      epktx = radius * sin(winkel)
      epkty = radius * cos(winkel)
      if (i == ni) then
         epktx = radius
         epkty = 0.0d0
      end if
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx + p1_x
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly + p1_y
   end do
end do

b = 3
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  45.0D0 * DBLE(I-1) /DBLE(NI-1) * pi
      spktx = DBLE(I-1)*RADIUS*0.3D0/DBLE(NI-1)
      spkty = RADIUS*0.3D0
      epktx = RADIUS * SIN(WINKEL)
      epkty = RADIUS * COS(WINKEL)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx + p1_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly + p1_y
   END DO
END DO

b = 4
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (45.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = radius * sin(winkel)
      spkty = radius * cos(winkel)
      epkty = dble(ni-i) * radius * 1.0d0 / dble(ni-1)
      epktx = radius * 1.5d0
      if (i == ni) then
         epktx = radius * 1.5d0
         epkty = 0.0d0
      end if
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx + p1_x
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly + p1_y
   end do
end do

b = 5
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  45.0D0 * DBLE(I-1) /DBLE(NI-1) * pi
      epktx = DBLE(I-1)*RADIUS*1.5D0/DBLE(NI-1)
      epkty = RADIUS*2.0D0- dble(i-1) * 1.0*radius/dble(ni-1)
      spktx = RADIUS * SIN(WINKEL)
      spkty = RADIUS * COS(WINKEL)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx + p1_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly + p1_y
   END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
b = 6
ni = blocks(b) % nPkts(1)
nj = blocks(b) % nPkts(2)
 k = 1
do j = 1, nj
   do i = 1,ni
      blocks(b) % xyzs(i,j,k,1) = -RADIUS * 0.3D0 + dble(i-1)*RADIUS*0.6d0/dble(ni-1) + p2_x
      blocks(b) % xyzs(i,j,k,2) = dble(j-1)*RADIUS*0.3d0/dble(nj-1) + p2_y
   end do
end do

b = 7
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (45.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      spkty = dble(ni-i) * radius * 0.3d0 / dble(ni-1)
      spktx = radius * 0.3d0
      epktx = radius * sin(winkel)
      epkty = radius * cos(winkel)
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx  + p2_x
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly  + p2_y
   end do
end do

b = 8
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (-45.0D0 + 90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      !winkel =  (90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      spktx = -RADIUS*0.3D0 + DBLE(I-1)*RADIUS*0.6D0/DBLE(NI-1)
      spkty = RADIUS*0.3D0
      epktx = RADIUS * SIN(WINKEL)
      epkty = RADIUS * COS(WINKEL)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx  + p2_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly  + p2_y
   END DO
END DO

b = 9
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (-90.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      spkty = dble(i-1) * radius * 0.3d0 / dble(ni-1)
      spktx = -radius * 0.3d0
      epktx = radius * sin(winkel)
      epkty = radius * cos(winkel)
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx  + p2_x
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly  + p2_y
   end do
end do
b = 10
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (45.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      w2 =  (60.0d0 + 30.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = radius * sin(winkel)
      spkty = radius * cos(winkel)
      epktx = radius_out1 * sin(w2) - p2_x
      epkty = radius_out1 * cos(w2) - p2_y
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx  + p2_x
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly  + p2_y
   end do
end do

b = 11
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (-45.0D0 + 90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      !winkel =  (90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      epktx = -RADIUS*1.5D0 + DBLE(I-1)*(radius_out1 * sin(pi * 60.0d0) - p2_x+ RADIUS*1.5d0) / DBLE(NI-1) 
      epkty = RADIUS + DBLE(I-1)*(radius_out1 * cos(60.0d0*pi)  - p2_y - RADIUS)/DBLE(NI-1)
      spktx = RADIUS * SIN(WINKEL)
      spkty = RADIUS * COS(WINKEL)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx  + p2_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly  + p2_y
   END DO
END DO

b = 12
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (-90.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      epkty = dble(i-1) * radius * 1.0d0 / dble(ni-1)
      epktx = -radius * 1.5d0
      spktx = radius * sin(winkel)
      spkty = radius * cos(winkel)
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx  + p2_x
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly  + p2_y
   end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
b = 13
ni = blocks(b) % nPkts(1)
nj = blocks(b) % nPkts(2)
 k = 1
do j = 1, nj
   do i = 1,ni
      blocks(b) % xyzs(i,j,k,1) = -radius_o*sqrt(2.0)+ dble(i-1)*sqrt(2.0)*radius_o/dble(ni-1) & 
                                                     + dble(j-1)*sqrt(2.0)*radius_o/dble(nj-1) + p3_x
      blocks(b) % xyzs(i,j,k,2) =           dble(j-1)*sqrt(2.0)*radius_o/dble(nj-1) & 
                                          - dble(i-1)*sqrt(2.0)*radius_o/dble(ni-1) + p3_y
   end do
end do

b = 14
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (+ 90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      !winkel =  (90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      spktx = DBLE(I-1)*sqrt(2.0)*radius_o/DBLE(NI-1)
      spkty = +radius_o* sqrt(2.0) - dble(i-1) * sqrt(2.0) * radius_o / dble(ni-1)
      epktx = RADIUS * SIN(WINKEL)
      epkty = RADIUS * COS(WINKEL)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx  + p3_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly  + p3_y
   END DO
END DO

b = 15
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (270.0d0 + 90.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = -radius_o* sqrt(2.0) + dble(i-1) * sqrt(2.0) * radius_o / dble(ni-1)
      spkty =  dble(i-1) * sqrt(2.0) * radius_o / dble(ni-1)
      epktx = radius * sin(winkel)
      epkty = radius * cos(winkel)
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx + p3_x 
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly + p3_y
   end do
end do

b = 16
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (180.0D0 + 90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      !winkel =  (90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      spktx = - dble(i-1) * sqrt(2.0) * radius_o / dble(ni-1)
      spkty = -radius_o* sqrt(2.0) + dble(i-1) * sqrt(2.0) * radius_o / dble(ni-1)
      epktx = RADIUS * SIN(WINKEL)
      epkty = RADIUS * COS(WINKEL)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx + p3_x 
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly + p3_y
   END DO
END DO

b = 17
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (90.0d0 + 90.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = +radius_o* sqrt(2.0) - dble(i-1) * sqrt(2.0) * radius_o / dble(ni-1)
      spkty = - dble(i-1) * sqrt(2.0) * radius_o / dble(ni-1)
      epktx = radius * sin(winkel)
      epkty = radius * cos(winkel)
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx + p3_x 
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly + p3_y
   end do
end do

b = 18
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (+ 90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      w2 =  (60.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = RADIUS * SIN(WINKEL)
      spkty = RADIUS * COS(WINKEL)
      epktx = radius_out1 * sin(w2) - p3_x
      epkty = radius_out1 * cos(w2) - p3_y
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx  + p3_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly  + p3_y
   END DO
END DO

b = 19
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (270.0d0 + 90.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = RADIUS * SIN(WINKEL)
      spkty = RADIUS * COS(WINKEL)
      epktx = -p3_x
      epkty = 2.0 * radius - p3_y + dble(i-1) * (radius_out1 - 2.0 * radius-0.0*p3_y) / dble(ni-1)
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx  + p3_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly  + p3_y
   END DO
END DO

b = 20
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (180.0D0 + 90.0D0 * DBLE(I-1) /DBLE(NI-1) )* pi
      spktx = RADIUS * SIN(WINKEL)
      spkty = RADIUS * COS(WINKEL)
      epktx = radius*1.5d0-DBLE(I-1)*RADIUS*1.5D0/DBLE(NI-1)-p3_x
      epkty = RADIUS*1.0D0+dble(i-1)*1.0*radius/dble(ni-1)-p3_y
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx + p3_x 
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly + p3_y
   END DO
END DO

b = 21
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (90.0d0 + 90.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = radius * sin(winkel)
      spkty = radius * cos(winkel)
      epktx = (radius_out1 * sin(pi * 60.0d0) - p3_x) - DBLE(I-1)*(radius_out1 * sin(pi * 60.0d0) - p3_x) / DBLE(NI-1) 
      epkty = radius_out1 * cos(60.0d0*pi) - p3_y - DBLE(I-1)*(radius_out1 * cos(60.0d0*pi)  - RADIUS)/DBLE(NI-1)
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx + p3_x 
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly + p3_y
   end do
end do

b = 22
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel =  (45.0d0 + 45.0d0 * dble(i-1) / dble(ni-1)) * pi
      w2 =  (60.0d0 + 30.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = radius_out1 * sin(w2) - p2_x
      spkty = radius_out1 * cos(w2) - p2_y
      epktx = radius_wall * sin(winkel) - p2_x
      epkty = radius_wall * cos(winkel) - p2_y
      lx = epktx-spktx
      ly = epkty-spkty

      blocks(b) % xyzs(i,j,k,1) = spktx + dble(j-1)/dble(nj-1) * lx  + p2_x
      blocks(b) % xyzs(i,j,k,2) = spkty + dble(j-1)/dble(nj-1) * ly  + p2_y
   end do
end do
b = 23
NI = blocks(b) % nPkts(1)
NJ = blocks(b) % nPkts(2)
do j = 1, nj
   do i = 1,ni
      winkel = (45.0d0 * dble(i-1) / dble(ni-1)) * pi
      w2 =  (60.0d0 * dble(i-1) / dble(ni-1)) * pi
      spktx = radius_out1 * sin(w2) - p3_x
      spkty = radius_out1 * cos(w2) - p3_y
      epktx = radius_wall * SIN(winkel) - p3_x
      epkty = radius_wall * COS(winkel) - p3_y
      lx = epktx-spktx
      ly = epkty-spkty
      blocks(b) % xyzs(i,j,k,1) = spktx + DBLE(J-1)/DBLE(NJ-1) * lx  + p3_x
      blocks(b) % xyzs(i,j,k,2) = spkty + DBLE(J-1)/DBLE(NJ-1) * ly  + p3_y
   END DO
END DO
call write_grid()
open(666,file="bc.cfg")
write(666,'(A)') "! Wall at South of all  and EAST OF BLOCK 1"
write(666,'(A)') "wall: 22N, 23N! Wall"
write(666,'(A)') "dn = 5E-5    ! Spacing of first Cell"
close(666)
write(*,*) "done"


end program
