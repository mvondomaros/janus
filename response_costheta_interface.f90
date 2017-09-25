program response_costheta_interface
implicit none

!! *** Calculates the response function of cos(theta) within both interfacial regions ***

integer, parameter :: dp = selected_real_kind(p=12)

integer, parameter :: stepmax = 5000

integer :: filemax
real(dp) :: ave_low, ave_up
character(len=250) :: output, buffy

integer :: ifile, istep
integer :: idummy
real(dp) :: rdummy
real(dp) :: dip_low, dip_up
real(dp) :: norm_low, norm_up

real(dp), dimension(stepmax) :: corr_low, corr_up

call getarg(1, buffy)
read(buffy, *) filemax
call getarg(2, buffy)
read(buffy, *) ave_low
call getarg(3, buffy)
read(buffy, *) ave_up
call getarg(4, output)

open(unit=5, action='read') ! stdin
corr_low = 0.0_dp
corr_up = 0.0_dp
do ifile = 1, filemax
    do istep = 1, stepmax
        read(5, *) idummy, rdummy, rdummy, dip_low, rdummy, rdummy, dip_up
        corr_low(istep) = corr_low(istep) + dip_low
        corr_up(istep) = corr_up(istep) + dip_up
    end do
end do
close(5)

corr_low = corr_low - filemax*ave_low
corr_up = corr_up - filemax*ave_up

norm_low = corr_low(1)
norm_up = corr_up(1)

corr_low = corr_low / norm_low
corr_up = corr_up / norm_up

open(unit=100, file=output, status='unknown', action='write')
do istep = 1, stepmax
    write(100, *) real(istep-1,dp)/100.0_dp, corr_low(istep), corr_up(istep)
end do

end program response_costheta_interface
