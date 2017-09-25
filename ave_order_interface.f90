program ave_order_interface
implicit none

! calculates the orientational order parameter q
! q = 1 - 3/8 * sum[j=1,3]sum[k=j+1,4]{cos psi + 1/3}**2
! within both interfacial regions

! Modified, to read STDIN

! parameters
integer, parameter :: dp = selected_real_kind(p=12)

! simulation settings
integer, parameter :: nsteps = 900000
integer, parameter :: nw = 839
integer, parameter :: nc = 767
real(dp), parameter :: zmax = 24.0_DP
real(dp), parameter :: zmin = 11.0_DP
real(dp), parameter :: xbox = 31.92800_DP
real(dp), parameter :: ybox = 31.19539_DP
real(dp), parameter :: zbox = 31.94800_DP
integer, parameter :: natoms = nw*3 + nc
real(dp), parameter :: xbox2 = xbox/2.0_dp
real(dp), parameter :: ybox2 = ybox/2.0_dp
real(dp), parameter :: zbox2 = zbox/2.0_dp
real(dp), parameter :: center = 16.0_DP

! types
type neighbor
    integer :: id
    real(dp) :: dist
end type neighbor

! variables with physical meaning
real(dp), dimension(nw) :: x, y, z
real(dp) :: xx, yy, zz
real(dp) :: xx1, yy1, zz1
real(dp) :: xx2, yy2, zz2
real(dp) :: angle, q, dist
type(neighbor), dimension(5) :: neigh
logical, dimension(nw) :: lowerslab, upperslab

! auxilliary variables
integer :: it
integer :: iw
integer :: jw
integer :: ia
integer :: idummy
integer :: i, j, k, l
character(len=120) :: input, output_low, output_up

! *** Execution section ***

! get CLI arguments
call getarg(1, output_low)
call getarg(2, output_up)

! open files
open(5)
open(unit=200, status='new', action='write', file=output_low)
open(unit=300, status='new', action='write', file=output_up)

! read the trajectory (stepwise)

do it = 1, nsteps
    call read_oxygens(5, x, y, z)
    ! create a mask for the slabs
    lowerslab = .false.
    upperslab= .false.
    where (z >= zmin .and. z < center)
        lowerslab = .true.
    end where
    where (z <= zmax .and. z > center)
        upperslab = .true.
    end where
    ! find the 4 closest neighbors of each atom
    ! and calculate q
    do iw = 1, nw
        if (lowerslab(iw)) then
            neigh = neighbor(0, 10*zbox**2)
            do jw = 1, nw
                if (jw == iw) cycle
                if (z(jw) < center) then ! ignore the lower part
                    xx = x(jw) - x(iw)
                    yy = y(jw) - y(iw)
                    zz = z(jw) - z(iw)
                    ! PBC
                    if ( xx > xbox2 ) xx = xx - xbox
                    if ( yy > ybox2 ) yy = yy - ybox
                    if ( zz > zbox2 ) zz = zz - zbox
                    if ( xx < -xbox2 ) xx = xx + xbox
                    if ( yy < -ybox2 ) yy = yy + ybox
                    if ( zz < -zbox2 ) zz = zz + zbox
                    ! calc distance
                    dist = xx**2 + yy**2 + zz**2
                    ! check for nearest neighbor
                    if (dist < neigh(5)%dist) then
                        neigh(5) = neighbor(jw, dist)
                        ! sort neighborlist (bubble sort)
                        call sort(neigh, 5)
                    end if
                end if
            end do
            ! neigh(1:4) should now contain the 4 closest neighbors
            ! calculate q
            q = 0.0_dp
            do j = 1, 3
                do k = j+1, 4
                    l = neigh(j)%id
                    xx1 = x(l) - x(iw)
                    yy1 = y(l) - y(iw)
                    zz1 = z(l) - z(iw)
                    ! PBC
                    if ( xx1 > xbox2 ) xx1 = xx1 - xbox
                    if ( yy1 > ybox2 ) yy1 = yy1 - ybox
                    if ( zz1 > zbox2 ) zz1 = zz1 - zbox
                    if ( xx1 < -xbox2 ) xx1 = xx1 + xbox
                    if ( yy1 < -ybox2 ) yy1 = yy1 + ybox
                    if ( zz1 < -zbox2 ) zz1 = zz1 + zbox
                    l = neigh(k)%id
                    xx2 = x(l) - x(iw)
                    yy2 = y(l) - y(iw)
                    zz2 = z(l) - z(iw)
                    ! PBC
                    if ( xx2 > xbox2 ) xx2 = xx2 - xbox
                    if ( yy2 > ybox2 ) yy2 = yy2 - ybox
                    if ( zz2 > zbox2 ) zz2 = zz2 - zbox
                    if ( xx2 < -xbox2 ) xx2 = xx2 + xbox
                    if ( yy2 < -ybox2 ) yy2 = yy2 + ybox
                    if ( zz2 < -zbox2 ) zz2 = zz2 + zbox
                    angle = xx1*xx2+yy1*yy2+zz1*zz2
                    angle = angle/sqrt(xx1**2+yy1**2+zz1**2)
                    angle = angle/sqrt(xx2**2+yy2**2+zz2**2)
                    q = q + (angle + 1.0_dp/3.0_dp)**2
                end do
            end do
            q = 1.0_dp - 3.0_dp/8.0_dp*q
            write(200,*) q
        end if
        if (upperslab(iw)) then
            neigh = neighbor(0, 10*zbox**2)
            do jw = 1, nw
                if (jw == iw) cycle
                if (z(jw) > center) then ! ignore the lower part
                    xx = x(jw) - x(iw)
                    yy = y(jw) - y(iw)
                    zz = z(jw) - z(iw)
                    ! PBC
                    if ( xx > xbox2 ) xx = xx - xbox
                    if ( yy > ybox2 ) yy = yy - ybox
                    if ( zz > zbox2 ) zz = zz - zbox
                    if ( xx < -xbox2 ) xx = xx + xbox
                    if ( yy < -ybox2 ) yy = yy + ybox
                    if ( zz < -zbox2 ) zz = zz + zbox
                    ! calc distance
                    dist = xx**2 + yy**2 + zz**2
                    ! check for nearest neighbor
                    if (dist < neigh(5)%dist) then
                        neigh(5) = neighbor(jw, dist)
                        ! sort neighborlist (bubble sort)
                        call sort(neigh, 5)
                    end if
                end if
            end do
            ! neigh(1:4) should now contain the 4 closest neighbors
            ! calculate q
            q = 0.0_dp
            do j = 1, 3
                do k = j+1, 4
                    l = neigh(j)%id
                    xx1 = x(l) - x(iw)
                    yy1 = y(l) - y(iw)
                    zz1 = z(l) - z(iw)
                    ! PBC
                    if ( xx1 > xbox2 ) xx1 = xx1 - xbox
                    if ( yy1 > ybox2 ) yy1 = yy1 - ybox
                    if ( zz1 > zbox2 ) zz1 = zz1 - zbox
                    if ( xx1 < -xbox2 ) xx1 = xx1 + xbox
                    if ( yy1 < -ybox2 ) yy1 = yy1 + ybox
                    if ( zz1 < -zbox2 ) zz1 = zz1 + zbox
                    l = neigh(k)%id
                    xx2 = x(l) - x(iw)
                    yy2 = y(l) - y(iw)
                    zz2 = z(l) - z(iw)
                    ! PBC
                    if ( xx2 > xbox2 ) xx2 = xx2 - xbox
                    if ( yy2 > ybox2 ) yy2 = yy2 - ybox
                    if ( zz2 > zbox2 ) zz2 = zz2 - zbox
                    if ( xx2 < -xbox2 ) xx2 = xx2 + xbox
                    if ( yy2 < -ybox2 ) yy2 = yy2 + ybox
                    if ( zz2 < -zbox2 ) zz2 = zz2 + zbox
                    angle = xx1*xx2+yy1*yy2+zz1*zz2
                    angle = angle/sqrt(xx1**2+yy1**2+zz1**2)
                    angle = angle/sqrt(xx2**2+yy2**2+zz2**2)
                    q = q + (angle + 1.0_dp/3.0_dp)**2
                end do
            end do
            q = 1.0_dp - 3.0_dp/8.0_dp*q
            write(300,*) q
        end if
    end do
end do

contains

    subroutine sort(neigh, n)
        implicit none

        integer, intent(in) :: n
        type(neighbor), dimension(n), intent(inout) :: neigh

        type(neighbor) :: temp
        integer :: i, j

        do i = n, 1, -1
            do j = 2, i
                if (neigh(j-1)%dist > neigh(j)%dist) then
                    temp = neigh(j)
                    neigh(j) = neigh(j-1)
                    neigh(j-1) = temp
                end if
            end do
        end do

    end subroutine sort

    subroutine read_oxygens(iunit, x, y, z)
    implicit none
    ! reads in the oxygen positions from LAMMPS output

    ! calling parameters
    integer, intent(in) :: iunit
    real(kind = dp), intent(out), dimension(nw) :: x, y, z

    integer :: i, j, k
    real(kind = dp), dimension(3*nw) :: inputx, inputy, inputz

    do i = 1, 9
        read(iunit, *)
    end do
    do i = 1, 3*nw
        read(iunit, *) j, k, inputx(j - nc), inputy(j - nc), inputz(j - nc)
    end do
    do i = 1, nw
        j = 3 * i - 2
        x(i) = inputx(j) * xbox
        y(i) = inputy(j) * ybox
        z(i) = inputz(j) * zbox
    end do
    end subroutine read_oxygens

end program ave_order_interface
