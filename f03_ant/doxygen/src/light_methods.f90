module light_methods

use real_precision

contains

!**********************************************************************!
!                                                                      !
!                     dayl :: light_methods                            !
!                     ---------------------                            !
!                                                                      !
! real(dp) function dayl(lat,day)                                               !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Calculates daylength (h) for the gridcell lat and day of year
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
real(dp) function dayl(lat,day)
!**********************************************************************!
implicit none
real(dp) :: del !< solar declination
real(dp) :: tem !< cosine of the hour angle at sunrise in radians
real(dp) :: has !< hour angle in degrees   
real(dp) :: lat
integer :: day
real(dp), parameter :: conv = 1.74532925e-2
!----------------------------------------------------------------------!

del=-23.4*cos(conv*360.0*(day + 10.0)/365.0)
tem=-tan(lat*conv)*tan(del*conv)

if (abs(tem)<1.0) then
  has  = acos(tem)/conv
  dayl = 2.0*has/15.0
elseif (tem>0.0) then
  dayl =  0.0
else
  dayl = 24.0
endif

end function dayl





!**********************************************************************!
!                                                                      !
!                     pfd :: light_methods                             !
!                     --------------------                             !
!                                                                      !
! SUBROUTINE pfd(lat,day,hrs,cloud,direct,diffuse,total)               !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Calculates total,direct and diffuse radiation for the day
!! @details Units oof outputs are mol/m2/sec
!! Needs interpolated daily values of cloud cover instead
!! of the same one for all the days in the month
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine pfd(lat,day,hrs,cloud,direct,diffuse,total)
!**********************************************************************!
implicit none
real(dp) :: hrs !< day length in hours
real(dp) :: del !< solar declination
real(dp) :: toa !< solar flux density average over day
real(dp) :: rlat !< latitude in radians
real(dp) :: dawn_angle !< Hour angle in radians it will integrate for
real(dp) :: total !< Average direct and diffuse radiation in mole/m2/sec
real(dp) :: direct !< Average direct radiation in mole/m2/sec
real(dp) :: diffuse !< Average diffuse radiation in mole/m2/sec
real(dp) :: lat,cloud,diffprop,alpha,beta, &
 sigma,coscst,clearness
integer :: day
real(dp), parameter :: conv = 1.74532925e-2, pi = 3.1415927
!----------------------------------------------------------------------!

if (hrs>1e-6) then
  del  = -23.4*cos(conv*360.0*(day+10.0)/365.0)
  del  = del*conv
  rlat = lat*conv

!     solar radiation in watts/m2 average over the day 
!     (between sun rise en sun set)

  dawn_angle = (hrs/2.0)*(2.0*pi/24.0)

!Top of atmosphere irradiance, average over the day (daylength).
! 1370 = solar constant in W/m2 
  toa = 1370.0*24.0/(pi*hrs)*(sin(rlat)*sin(del)*dawn_angle + &
 cos(rlat)*cos(del)*sin(dawn_angle))

!Account for cloud cover

  coscst = 1.0-1.3614*cos(rlat)
  alpha = 18.0 - 64.884*coscst
  alpha = alpha*4.1842*1.0d4 ! convert from cal/cm to J/m2
  alpha = alpha/(24.0*3600.0) ! convert into w/m2
  
  beta = 0.682 - 0.3183*coscst

  sigma = 0.02*log(max(cloud,0.001))+0.03259

!Calculate the total irradiance
  total = toa*(beta-sigma*cloud*10.0)-alpha
  if (total<0.0) total=0.0

  clearness= total/toa

!Calculate diffuse irradiance (from Forest ETP WG)

  if(clearness<0.07) then
     diffprop = 1.0
  else if(clearness<0.35) then
     diffprop = 1.0-2.3*(clearness-0.07)**2
  else if(clearness<0.75) then
     diffprop = 1.33-1.46*clearness
  else
     diffprop = 0.23
  endif

  diffuse=total*diffprop

! PAR is about 48% of the total irradiance. A better formula could be
! found.
  total = 0.48*total
  diffuse = 0.48*diffuse ! not sure this is really right (Ghislain)

!Convert in mol/m2/sec 1W/m2=4.6e-6 mol/s/m2 for PAR radiation

  total = 4.6d-6*total
  diffuse = 4.6d-6*diffuse
  direct= total - diffuse
else
  total=0.0
  diffuse=0.0
  direct=0.0
endif

end subroutine pfd



!**********************************************************************!
!                                                                      !
!                     pfd_ant :: light_methods                         !
!                     ------------------------                         !
!                                                                      !
! SUBROUTINE pfd(lat,day,hrs,cloud,direct,diffuse,total)               !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Calculates total,direct and diffuse radiation for the day
!! @details Units oof outputs are mol/m2/sec
!! Needs interpolated daily values of cloud cover instead
!! of the same one for all the days in the month
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine pfd_ant(lat,day,hrs,cloud,direct,diffuse,total,swr, &
 read_par,subd_par,t,total_t,calc_zen,cos_zen)
!**********************************************************************!
REAL(dp) :: lat,hrs,del,rlat,toa,cloud,diffprop,alpha,beta
REAL(dp) :: sigma,dawn_angle,coscst,direct,diffuse,total,clearness
REAL(dp) :: swr,ts,h,toa_h,tem,dJ_R,dJ_K,cos_zen
INTEGER :: day,t,total_t
real(dp) :: conv = 1.74532925E-2
real(dp) :: pi = 3.1415927
real(dp) :: solar_const = 1370.0 ! solar constant (W/m2)
logical :: read_par,subd_par,calc_zen
!----------------------------------------------------------------------!

if (hrs.GT.1E-6) then
! declination angle
  del  = -23.4*cos(conv*360.0*(day+10.0)/365.0)
  del  = del*conv
  rlat = lat*conv

! calculate dawn angle
  dawn_angle = (hrs/2.0)*(2.0*pi/24.0)

! calculate daytime mean cos of the zenith angle  
  cos_zen = 1/dawn_angle*(sin(rlat)*sin(del)*dawn_angle + &
 cos(rlat)*cos(del)*sin(dawn_angle) )

! top of atmosphere irradiance
! - watts/m2 average over the day 
! - between sun rise and sun set
! toa = solar_const * cos_zen ! daytime mean irradiance (W/m2)
  toa = solar_const *24.0/(pi*hrs)*(sin(rlat)*sin(del)*dawn_angle + &
 cos(rlat)*cos(del)*sin(dawn_angle))

! calculate the mean daylight top of canopy (total) shortwave radiation
  if (read_par) then
!read from the input dataset as mean SWR over 24hrs
! - scale to mean over daylight hours only
    total = 24.0*swr/hrs   
  else
! calculate based on top of atmosphere value and cloud cover
! - also calculated from Forest ETP 
    coscst = 1.0-1.3614*cos(rlat)
    alpha  = 18.0 - 64.884*coscst
    alpha  = alpha*4.1842*1.0d4     ! convert from cal/cm to J/m2
    alpha  = alpha/(24.0*3600.0)  ! convert into w/m2  
    beta   = 0.682 - 0.3183*coscst
    sigma  = 0.02*log(max(cloud,0.001)) + 0.03259

    total  = toa*(beta-sigma*cloud*10.0) - alpha
  endif

  if (subd_par) then
! sub-daily variation in SWR 

! calculate using hour angle (h) 
! - scale by (instantaneous insolation at h)/(mean daylight insolation) law of cosines
! - when t = 0, h is 0 and is solar noon
! - when t = total_t, h is the dawn angle i.e. sunrise
    h       = dawn_angle * real(t)/real(total_t)
! calculate cos of the zenith angle
    cos_zen = sin(rlat)*sin(del) + cos(rlat)*cos(del)*cos(h)
! calculate instantaneous insolation at hour angle h
    toa_h   = solar_const * cos_zen 
! toa_h   = solar_const * ( sin(rlat)*sin(del) + 
! & cos(rlat)*cos(del)*cos(h) )
!          if(day.eq.160) print*, h,toa,toa_h,toa_h/toa
! scale SWR (total) from daytime mean to instantaneous 
! - scale by the ratio of toa instantaneous insolation to daytime mean
    if(toa.gt.0.0) then
      total = total*toa_h/toa
    else
      total = 0.0
    endif
! if(day.eq.160) print*, total 
! reset top of atmosphere to instantaneous value
    toa = toa_h

  endif

  if (total.LT.0.0) total=0.0
  if(toa.le.0.0) then
    clearness = 0.0
  else
    clearness = total/toa
  endif

! calculate diffuse irradiance (from Spitters etal 1986)
  IF(subd_par) THEN
! hourly data
    dJ_R = 0.847 - 1.61*cos_zen + 1.04*cos_zen**2
    dJ_K = (1.47-dJ_R)/1.66
    IF(clearness.LT.0.22) THEN
      diffprop = 1.0
    ELSE IF(clearness.LT.0.35) THEN
      diffprop = 1.0-6.4*(clearness-0.22)**2
    ELSE IF(clearness.LT.dJ_K) THEN
      diffprop = 1.47-1.66*clearness
    ELSE
      diffprop = dJ_R
    ENDIF
  ELSE
! daily data
    IF(clearness.LT.0.07) THEN
      diffprop = 1.0
    ELSE IF(clearness.LT.0.35) THEN
      diffprop = 1.0-2.3*(clearness-0.07)**2
    ELSE IF(clearness.LT.0.75) THEN
      diffprop = 1.33-1.46*clearness
    ELSE
      diffprop = 0.23
    ENDIF
  ENDIF

  diffuse = total*diffprop

! PAR is about 48% of the total irradiance. A better formula could be found.
  total   = 0.48*total
  diffuse = 0.48*diffuse ! not sure this is really right
                                 ! ... (Ghislain)
! convert to mol/m2/sec from W/m2
  total   = 4.6d-6*total
  diffuse = 4.6d-6*diffuse
  direct  = total - diffuse

else

  total   = 0.0
  diffuse = 0.0
  direct  = 0.0

endif

! set zenith angle to zero if requested (originaly SDGVM default)
if (.not.calc_zen) cos_zen = 1.0

end subroutine pfd_ant




!**********************************************************************!
!                                                                      !
!               pfd_without_cloud :: light_methods                     !
!               ----------------------------------                     !
!                                                                      !
! real(dp) pfd_without_cloud(lat,day,hrs)                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Photon Flux Density
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
real(dp) function pfd_without_cloud(lat,day,hrs)
!**********************************************************************!
real(dp) :: lat,hrs,del,rlat,r,sinbt,cloud,st,sdiff,trans
integer :: day
real(dp), parameter :: conv = 1.74532925e-2, pi = 3.1415927
!----------------------------------------------------------------------!

cloud = 0.50

if (hrs>1e-6) then
  del  = -23.4*cos(conv*360.0*(day+10.0)/365.0)
  del  = del*conv
  rlat = lat*conv

  trans = 0.3 + (0.35*(1.0 - cloud))

! solar radiation in watts/m2 at midday
  sinbt = sin(rlat)*sin(del) + cos(rlat)*cos(del)

  if (sinbt>1e-6) then
    r = 1370.0*trans**(1.0/sinbt)*sinbt
    sdiff = 0.3*(1.0 - trans**(1.0/sinbt))*1370.0*sinbt
    st = r + sdiff
    pfd_without_cloud = 2.0*st/pi*(0.5 + cos(rlat)/2.0)*2.1/1.0e+6
  else
   pfd_without_cloud = 0.0
  endif
else
  pfd_without_cloud = 0.0
endif

end function pfd_without_cloud





!**********************************************************************!
!                                                                      !
!                     pfd2 :: light_methods                            !
!                     ---------------------                            !
!                                                                      !
! real(dp) function pfd2(lat,day,hrs)                                           !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Photon Flux Density
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
real(dp) function pfd2(lat,day,hrs)
!**********************************************************************!
real(dp) :: lat,hrs,del,rlat,r,sinbt,rad
integer :: day
real(dp), parameter :: conv = 1.74532925e-2, pi = 3.1415927
!----------------------------------------------------------------------!

if (hrs>1e-6) then
  del  = -23.4*cos(conv*360.0*(day+10.0)/365.0)
  del  = del*conv
  rlat = lat*conv

! solar radiation in watts/m2 at midday
  sinbt = sin(rlat)*sin(del) + cos(rlat)*cos(del)
  r = 1360.0*0.7**1.32*sinbt

! convert to Langleys/day
  r = r/0.485

! account for cloud cover (assume 0.0)
  rad = r*(-0.21843*20.0 + 58.94408)/100.0

! solve for daily maximum (mol/m2/sec)
  pfd2 = (rad*41868.0*pi)/(hrs*2.0*3600.0)
  pfd2 = pfd2*0.5*4.255e-6
else
  pfd2 = 0.0
endif

end function pfd2





!**********************************************************************!
!                                                                      !
!                     pfds :: light_methods                            !
!                     ---------------------                            !
!                                                                      !
! real(dp) function pfds(tdss,hrs)                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Photon Flux Density
!! @details From Total Downward Shortwave at the Surface (TDSS) -
!! for use with GCM data.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
real(dp) function pfds(tdss,hrs)
!**********************************************************************!
real(dp) :: tdss,hrs
real(dp), parameter :: pi = 3.1415927
!----------------------------------------------------------------------!

if (hrs>1e-6) then
!Solve for daily maximum (mol/m2/sec)
  pfds = (tdss*pi*24.0)/(hrs*2.0)*4.255e-6
!  pfds = pfds*0.5*4.255E-6
else
  pfds = 0.0
endif

end function pfds



end module light_methods

