module hydrology_methods

use real_precision
use system_state
use site_parameters
use tuning_parameters
use pft_parameters
use input_file

implicit none

contains

!**********************************************************************!
!                                                                      !
!                       hydrology :: hydrology_methods                 !
!                       ------------------------------                 !
!                                                                      !
! subroutine hydrology(adp,sfc,sw,sswc,awl,kd,kx,eemm,etmm,pet2,prc,   !
! s1in,t,rlai,evap,tran,roff,interc,evbs,f2,f3,ft)                     !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Inclusion of century water model (which is a 'bucket model)
!! using doly evaporation and interception calculations.
!! @details lsfc is the field capacity for each layer
!! lsw is the wilting point for each layer
!! awl is the relatie root density for each layer 
!! ladp is the soil depth in mm for each layer
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine hydrology(adp,sfc,sw,sswc,awl,kd,kx,eemm,etmm,pet2,prc, &
 s1in,t,rlai,evap,tran,roff,interc,evbs,f2,f3,ft)
!**********************************************************************!
real(dp) :: adp(4),sfc(4),sw(4),sswc(4),awl(4),kd,kx,tf(13)
real(dp) :: ladp(4),lsfc(4),lsw(4),lsswc(4),s1in
real(dp) :: eemm,etmm,pet2,dp2,t,rlai,evap,tran,roff,f2,f3
real(dp) :: bst,kf,interc,ms,fs,rwc(4),w(4),rem,f1,ds,st,ws,sl,sf,sd
real(dp) :: fav,ev,ss,sums,evbs,pet3,ans1,tfscale,prc
integer :: lai,ft,i
!----------------------------------------------------------------------!


! dp2 at this point holds precip
dp2 = prc

!----------------------------------------------------------------------!
! SET parameter for THROUGHFALL.                                       !
!----------------------------------------------------------------------!
! tf holds the parameters for canopy interception at
! various lai values.These values are subtracted from 1 to get canopy
! interception so e.g. 1-tf(2)=0.05 would meen 5% of interception at
! lai=2
tf(1) = 1.0
tf(2) = 0.95
tf(3) = 0.935
tf(4) = 0.92
tf(5) = 0.905
tf(6) = 0.89
tf(7) = 0.875
tf(8) = 0.86
tf(9) = 0.845
tf(10)= 0.83
tf(11)= 0.815
tf(12)= 0.80
tf(13)= 0.785

! Scale for interception parameter
tfscale = 2.75
do i=1,13
  tf(i) = 1.0 - tfscale*(1.0 - tf(i))
enddo

! kd is the fraction of excess water flowing to deep storage
! and kf the one that doesn't.kd defined in wsparam
kf = 1.0 - kd

rem = rlai - int(rlai)
lai = int(rlai) + 1

! Why is pet2 multiplied by 3?pet2 holds the eemm
pet3 = pet2
pet2 = 3.0*pet2

!----------------------------------------------------------------------!
! Adjustment for the 'CITY' functional type.                           !
!----------------------------------------------------------------------!
if (pft(ft)%itag==2) then
  do i=1,4
    ladp(i) = adp(i)*tgp%p_city_dep/ssp%soil_depth/10.0
    lsfc(i) = sfc(i)*tgp%p_city_dep/ssp%soil_depth/10.0
    lsw(i) = sw(i)*tgp%p_city_dep/ssp%soil_depth/10.0
    lsswc(i) = sswc(i)*tgp%p_city_dep/ssp%soil_depth/10.0
  enddo
else
  ! Gets soil layer depth,field capacity,wilting point and saturation in mm
  do i=1,4
    ladp(i) = adp(i)
    lsfc(i) = sfc(i)
    lsw(i) = sw(i)
    lsswc(i) = sswc(i)
  enddo
endif

if (t<0.0) then
  s1in = 0.0
else
  s1in = dp2
endif

! If temp<0 then precip becomes snow and added to the snow layer
if (t<0.0) then
  ssv(ft)%snow = ssv(ft)%snow + dp2
  evap = 0.0
  interc = 0.0
else
  ! Canopy interception water loss (evap mm day-1).                             
  if (rlai>0) then
    ! Water interception by the canopy based on lai parameters
    interc = dp2*(1.0 - (tf(lai) + rem*(tf(lai+1) - tf(lai))))
    ! Intercepted water as a function of temperature,why?
    ! At 15C, factor becomes 1.O otherwise its always smaller
    interc = interc*min(1.0,(0.5 + t/32.0))
    ! evap holds canopy evaporation.At this point it is set equal to eemm
    evap = eemm
    ! This if is always false cause pet2=3*eemm 
    if (evap>pet2) evap = pet2
    ! Set canopy evap to intercept if it exceeds it.You can't have more canopy evap
    ! if the water is not on the canopy
    if (evap>interc)  evap = interc
    ! The remaining eemm when intercept evaporation is subtracted
    pet2 = pet2 - evap
    ! The remaining precip after intercept evaporation is subtracted
    ! This will move to the soil.It is assumed here that intercept is also moved
    ! to the soil as only intercept evaporation is subtracted.
    dp2 = dp2 - evap
    !interc = evap
  else
    evap = 0.0
    interc = 0.0
  endif

  ! If there is still snpw then add precip to liquid snow
  if (ssv(ft)%snow>0.0) then
    ssv(ft)%l_snow = ssv(ft)%l_snow + dp2
  ! otherwise to first soil layer
  else
    ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + dp2
  endif
endif

! Snow sublimation and snow melt
if (ssv(ft)%snow>0.0) then
  ! Snow sublimation as a function of remaining eemm
  sl = 0.85*pet2/50.0
  ! If it exceeds available snow
  if (sl>ssv(ft)%snow) then
    sl = ssv(ft)%snow
    pet2 = pet2 - sl
    ssv(ft)%snow = 0.0
  ! If snow still exists after sublimation,update remaining
  ! snow and eemm
  else
    ssv(ft)%snow = ssv(ft)%snow - sl
    pet2 = pet2 - sl
    ms = 0.0
    ! If temperature greater than 0,calculate snow melt
    if (t>0.0) then
      ms = t*40.0/real(30)
      if (ms>ssv(ft)%snow) then
        ms = ssv(ft)%snow
      endif
      ssv(ft)%snow = ssv(ft)%snow - ms
      ssv(ft)%l_snow = ssv(ft)%l_snow + ms
    endif
  endif
else
  sl = 0.0
endif


! precip has already been added to the first soil layer above.
! Here it adds snow melt
fs = 0.0
if ((ssv(ft)%l_snow>0.05*(ssv(ft)%l_snow + ssv(ft)%snow)).and.(t>0.0)) &
 then
  fs = ssv(ft)%l_snow - 0.05*(ssv(ft)%l_snow + ssv(ft)%snow)
  ssv(ft)%l_snow = 0.05*(ssv(ft)%l_snow + ssv(ft)%snow)
  ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + fs
  s1in = s1in + fs
endif

! Checks each soil layer.If the water holding capacity of the previous
! layer was exceeded,it sends remaining water to the layer below.
!----------------------------------------------------------------------!
! Soil water 150-300mm (ssv(ft)%soil_h2o(2)).                          !
!----------------------------------------------------------------------!
f1 = 0.0
if (ssv(ft)%soil_h2o(1)>lsfc(1)) then
  f1 = ssv(ft)%soil_h2o(1) - lsfc(1)
  ssv(ft)%soil_h2o(1) = lsfc(1)
endif
ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + f1

!----------------------------------------------------------------------!
! Soil water 150-300mm (ssv(ft)%soil_h2o(2)).                          !
!----------------------------------------------------------------------!
f2 = 0.0
if (ssv(ft)%soil_h2o(2)>lsfc(2)) then
  f2 = ssv(ft)%soil_h2o(2) - lsfc(2)
  ssv(ft)%soil_h2o(2) = lsfc(2)
endif
ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + f2

!----------------------------------------------------------------------!
! Soil water 300-450mm (ssv(ft)%soil_h2o(3)).                          !
!----------------------------------------------------------------------!
f3 = 0.0
if (ssv(ft)%soil_h2o(3)>lsfc(3)) then
  f3 = ssv(ft)%soil_h2o(3) - lsfc(3)
  ssv(ft)%soil_h2o(3) = lsfc(3)
endif
ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) + f3

!----------------------------------------------------------------------!
! Fill up to saturated water content from bottom up.                   !
!----------------------------------------------------------------------!
! Also calculates runoff (roff)

roff = 0.0
! Threshold over which I have saturation water for each layer
! Move the excess from bottom up until you reach the first layer where
! the excess becomes runoff
ans1 = lsfc(4) + tgp%p_roff2*(lsswc(4) - lsfc(4))
if (ssv(ft)%soil_h2o(4)>ans1) then
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4) - ans1
  ssv(ft)%soil_h2o(4) = ans1
  ans1 = lsfc(3) + tgp%p_roff2*(lsswc(3) - lsfc(3))
  if (ssv(ft)%soil_h2o(3)>ans1) then
    ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(3) - &
 ans1
    ssv(ft)%soil_h2o(3) = ans1
    ans1 = lsfc(2) + tgp%p_roff2*(lsswc(2) - lsfc(2))
    if (ssv(ft)%soil_h2o(2)>ans1) then
      ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + &
 ssv(ft)%soil_h2o(2) - ans1
      ssv(ft)%soil_h2o(2) = ans1
      ans1 = lsfc(1) + tgp%p_roff2*(lsswc(1) - lsfc(1))
      if (ssv(ft)%soil_h2o(1)>ans1) then
        roff = ssv(ft)%soil_h2o(1) - ans1
        ssv(ft)%soil_h2o(1) = ans1
      endif
    endif
  endif
endif

!----------------------------------------------------------------------!
! Deep storage (ds - mm day-1) and stream (st - mm day-1).             !
!----------------------------------------------------------------------!
sf = 0.0
sd = 0.0
if (ssv(ft)%soil_h2o(4)>lsfc(4)) then
  sf = kf*(ssv(ft)%soil_h2o(4) - lsfc(4))*tgp%p_roff
  sd = kd*(ssv(ft)%soil_h2o(4) - lsfc(4))*tgp%p_roff
  roff = roff + sf + sd
endif

ds = 0.0
st = 0.0

ss = ds*kx
ds = ds + sd - ss
st = st + sf + ss
ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) - sf - sd

!----------------------------------------------------------------------!
! Calculation of transpiration (tran - mm day-1).                      !
!----------------------------------------------------------------------!
! If the soil water content for each layer is greater than the wilting
! point then it calculates the relative water content for each layer
! fav is not used
! Calculates relative water content for each layer except the 1st,
! where we consider no water is removed from the plant for transpiration
fav = 0.0

if (ssv(ft)%soil_h2o(1)>lsw(1)) then
  if (inp%run%s070607) then
    fav = fav + ssv(ft)%soil_h2o(1) - lsw(1)
    rwc(1) = (ssv(ft)%soil_h2o(1) - lsw(1))/(lsfc(1)-lsw(1))
  else
    rwc(1) = 0.0
  endif
else
  rwc(1) = 0.0
endif

if (ssv(ft)%soil_h2o(2)>lsw(2)) then
  fav = fav + ssv(ft)%soil_h2o(2) - lsw(2)
  rwc(2) = (ssv(ft)%soil_h2o(2) - lsw(2))/(lsfc(2)-lsw(2))
else
  rwc(2) = 0.0
endif

if (ssv(ft)%soil_h2o(3)>lsw(3)) then
  fav = fav + ssv(ft)%soil_h2o(3) - lsw(3)
  rwc(3) = (ssv(ft)%soil_h2o(3) - lsw(3))/(lsfc(3)-lsw(3))
else
  rwc(3) = 0.0
endif

if (ssv(ft)%soil_h2o(4)>lsw(4)) then
  fav = fav + ssv(ft)%soil_h2o(4) - lsw(4)
  rwc(4) = (ssv(ft)%soil_h2o(4) - lsw(4))/(lsfc(4) - lsw(4))
else
  rwc(4) = 0.0
endif


! rwc is the relative water content calculated above
! and awl the relative root density for each layer
w(1) = rwc(1)*awl(1)*ladp(1)
w(2) = rwc(2)*awl(2)*ladp(2)
w(3) = rwc(3)*awl(3)*ladp(3)
w(4) = rwc(4)*awl(4)*ladp(4)
ws = w(1) + w(2) + w(3) + w(4)

! If the transpiration is greater than the potential
! then it sets transpiration to potential.It also subtract
! what is transpired from potential.

if (inp%run%s070607) then
! Think Euripedes added this, it wan't in my original but I've stuck 
! with it for the s070807 version.
  tran = 2*etmm
else
  tran = etmm
endif

if (tran>pet2)  tran = pet2
pet2 = pet2 - tran
! Removes water for each layer for transpiration in
! proportion to the relative water content calculated above
if (ws>1e-6) then
  ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) - tran*w(1)/ws
  ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) - tran*w(2)/ws
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) - tran*w(3)/ws
  ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) - tran*w(4)/ws
else
  sums = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4)
  if (sums>1e-6) then
     ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) - &
 ssv(ft)%soil_h2o(1)*tran/sums
     ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) - &
 ssv(ft)%soil_h2o(2)*tran/sums
     ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) - &
 ssv(ft)%soil_h2o(3)*tran/sums
     ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) - &
 ssv(ft)%soil_h2o(4)*tran/sums
! changed by Ghislain 6/10/03
  else
     ssv(ft)%soil_h2o(2) = 0.0
     ssv(ft)%soil_h2o(3) = 0.0
     ssv(ft)%soil_h2o(4) = 0.0
     tran = 0.0
  endif
endif

! Safeties so water content in each layer won't drop below zero
! It draws water from the layer below to balance it to zero
if (ssv(ft)%soil_h2o(1)<0.0) then
  ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(1)
  ssv(ft)%soil_h2o(1) = 0.0
endif
if (ssv(ft)%soil_h2o(2)<0.0) then
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(2)
  ssv(ft)%soil_h2o(2) = 0.0
endif
if (ssv(ft)%soil_h2o(3)<0.0) then
  ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) + ssv(ft)%soil_h2o(3)
  ssv(ft)%soil_h2o(3) = 0.0
endif

if (ssv(ft)%soil_h2o(4)<0.0) then
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4)
  ssv(ft)%soil_h2o(4) = 0.0
  if (ssv(ft)%soil_h2o(3)<0.0) then
    ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(3)
    ssv(ft)%soil_h2o(3) = 0.0
    if (ssv(ft)%soil_h2o(2)<0.0) then
      ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + ssv(ft)%soil_h2o(2)
      ssv(ft)%soil_h2o(2) = 0.0
      if (ssv(ft)%soil_h2o(1)<0.0) then
        tran = tran + ssv(ft)%soil_h2o(1)
        if (tran<0.0)  tran = 0.0
        ssv(ft)%soil_h2o(1) = 0.0
        write(*,*) 'No water at all.'
      endif
    endif
  endif        
endif

!----------------------------------------------------------------------!
! Calculation of bare soil evaporation 'evbs' (mm day-1).              !
!----------------------------------------------------------------------!
! evbs is a function of temperature and eemm
ev = (ssv(ft)%soil_h2o(2) - lsw(2))/(lsfc(2) - lsw(2))
if (ev<0.0)  ev = 0.0
bst = 0.0
if (t>0) bst = (t/16.0)
evbs = ev*tgp%p_bs*0.33*pet3*1.3*bst

! Subtracts from available eemm
if (evbs>pet2)  evbs = pet2
pet2 = pet2 - evbs

! Removes soil evaporation from 1st soil layer
if (evbs>ssv(ft)%soil_h2o(1)) then
  evbs = ssv(ft)%soil_h2o(1)
  ssv(ft)%soil_h2o(1) = 0.0
else
  ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) - evbs
endif

! Sums evaporation from all sources. Canopy, sublimation and soil
evap = evap + sl + evbs

! Update soil moisture trigger for budburst
do i=1,29
  ssv(ssp%cohort)%sm_trig(31-i) = ssv(ssp%cohort)%sm_trig(30-i)
enddo
ssv(ssp%cohort)%sm_trig(1) = (s1in - 0.0*evbs)


end subroutine hydrology



end module hydrology_methods

