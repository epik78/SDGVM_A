!> @brief Collection of subroutines to read in data.
!! @details
!! @author Mark Lomas
!! @date July 2016

module daily_step

use real_precision
use phenology_methods
use productivity_methods
use misc_values
use site_parameters
use tuning_parameters
use system_state

implicit none

contains

!**********************************************************************!
!                                                                      !
!                          dailyStep :: daily_step                     !
!                          -----------------------                     !
!                                                                      !
!           (D)ynamic gl(O)ba(L) phytogeograph(Y) model                !
! SUBROUTINE DOLYDAY(tmp,prc,hum,ca,soilc,soiln,minn,adp,sfc,sw,sswc,  !
! awl,kd,kx,daygpp,resp_l,rlai,evap,tran,roff,interc,evbs,flow1,flow2, !
! year,mnth,day,pet,ht,thty_dys,ft,lmor_sc,nleaf,leaflitter,hrs,q,     !
! qdirect,qdiff,fpr,canga,gsn,rn,check_closure)                        !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! 
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine dailyStep(tmp,prc,hum,cld,ca,soilc,soiln,minn, &
 kd,kx,daygpp,resp_l,rlai,ht,ft,nleaf,hrs,q, &
 qdirect,qdiff,fpr,tleaf_n,tleaf_p,canga,gsn,rn,ce_light,ce_ci,ce_t, &
 ce_maxlight,ce_ga,ce_rh,check_closure,par_loops,lat,year,mnth,day, &
 thty_dys,gs_type,swr)
!**********************************************************************!
real(dp), parameter :: oi = 21000.0
real(dp) :: tmp,prc,hum,cld,ca, &
 t,suma,laiinc, &
 soilc,soiln,rh,tk,rd(12),rn,q,qdiff, &
 qdirect,hrs,canga,amx,gsn, &
 soilw,maxc,soil2g,can2a,daygpp, &
 gsum,ht, &
 amax,can2g,cangs,rlai, &
 p,canres,rem, &
 kd, &
 kx,minn,windspeed,resp_l, &
 npp_eff,nppstore, &
 fpr, &
 nleaf,old_total_carbon, &
 total_carbon,lat,tleaf_n,tleaf_p,swr
real(dp), dimension(30,12) :: ce_light,ce_ci,ce_ga,ce_maxlight
real(dp), dimension(30) :: ce_t,ce_rh
integer c3c4,thty_dys,ft, &
 mnth,i,lai,day,year, &
 co,par_loops,gs_type
logical veg,check_closure
!----------------------------------------------------------------------!

co = ssp%cohort

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if ((check_closure).and.(pft(co)%sla > 0.0)) then
  old_total_carbon = ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + &
 ssv(co)%nppstore(1) + ssv(co)%stem%tot(1) + ssv(co)%root%tot(1) + &
 ssv(co)%bio(1) + ssv(co)%bio(2)
endif
!----------------------------------------------------------------------!


if (pft(ft)%c3c4) then
  c3c4 = 1
else
  c3c4 = 0
endif 


p = 101325.0

! Real lai,the remainder and the whole plus one 
rlai = ssv(ssp%cohort)%lai%tot(1)
rem = rlai - int(rlai)
lai = int(rlai) + 1

!----------------------------------------------------------------------!
! INITIALISATIONS.                                                     !
!----------------------------------------------------------------------!
! Night-time respiration ? CHECK
rd = 0.82e-6

! maxc, moisture response parameter
! maxc is is the maximum influence of soil water on stomatal conductance.
! Goes into nppcalc
if (inp%run%s070607) then
  maxc = 690.0/622.6
else
  maxc = 1.0
endif


! initialisation a faire a chaque fois
if (pft(ft)%phen == 0) then
   veg = .false.
else
   veg = .true.
endif

! unsused anymore... so not important if set to this value every day. 
! Better: should be remove from nppcalc
amx  = -1000.0
amax = -1000.0

!----------------------------------------------------------------------!

soilw = ssv(ft)%soil_h2o(1) + ssv(ft)%soil_h2o(2) + &
 ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4)
soil2g = soilw/(ssp%soil_depth*10.0)

if (soil2g>1.0) soil2g = 1.0

!----------------------------------------------------------------------!
! Bare ground IF statement                                             !
!----------------------------------------------------------------------!
if (.not.(veg)) rlai = 0.0
!----------------------------------------------------------------------!

t = tmp
tk = 273.0 + t
rh = hum

if (rh>95.0)  rh=95.0
if (rh<30.0)  rh=30.0


t = tmp
! From Ziang et al., 2015 "Empirical estimation of daytine net radiation..."
! Net radiation as a function of incoming shortwave radiation,temperature and
! day of the year
! q is the PAR in mol/m2/sec calculated in PFD_ant.Here it becomes shortwave ratiation in W/m2
rn = 0.721*((1/0.48)*q/4.6d-6) + 0.777*t - 301.420*(1 + 0.033*cos(2*3.1415*no_day(year,mnth,day,thty_dys)/365)) + 296.842

! Old method calculating rn
!rn = 0.96*(q*1000000.0/4.0 + 208.0 + 6.0*t)
!rn = rn*0.52


!----------------------------------------------------------------------!
! canga=k^2 u / (log[(z-d)/z0])^2
! k=von Karman constant. k=0.41
! z=reference height
! d=zero plane displacement
! z0=roughness length
!----------------------------------------------------------------------!
! Boundary conductance (canga) derived from Eq.24 and Eq.25.
! Height (ht) is constant instead of being a function of lai as in Eq.25
! Also its a function of windspeed while in Eq.24 it's not although
! Eq.24 might be missing wind speed,check wording.
windspeed= 5.0 ! in m/s
canga = 0.168*windspeed/log((200.0 - 0.7*ht)/(0.1*ht))**2

!----------------------------------------------------------------------!

npp_eff = 0.0
!sum = 0.0
!eff_dec = 0.75
!eff_age = 90.0
!leafls   = pft(ft)%lls
!DO i=1,leafls
!  IF (i<=eff_age) THEN
!    eff = 1.0
!  ELSE
!    eff = eff_dec + real(leafls - i)*(1.0 - eff_dec)/(real(leafls) - &
! eff_age + 1.0)
!  ENDIF
!  npp_eff = npp_eff + leafv(i)*eff
!  sum = sum + leafv(i)
!ENDDO

!IF (abs(sum-rlai)>0.001) THEN
!  WRITE(*,*) 'leafv not = rlai ',SUM,RLAI,mnth,day,nppstore
!  STOP
!ENDIF
!IF (sum>0.0) THEN
!  npp_eff = npp_eff/sum
!ELSE
!  npp_eff = 0.0
!ENDIF
npp_eff = 1.0

if (veg) then
  call NPPCALC(npp_eff,c3c4,maxc,soilc,soiln,minn,soil2g,ssp%wilt, &
 ssp%field,rd,rlai,t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma,amx, &
 amax,gsum,hrs,canga/1.3,p,nleaf,fpr,tleaf_n,tleaf_p,ce_light,ce_ci,ce_t, &
 ce_maxlight,ce_ga,ce_rh,ft,par_loops,cld,lat,year,mnth,day,thty_dys,gs_type,swr)

!----------------------------------------------------------------------!
! Calculate the canopy conductance.                                    !
!----------------------------------------------------------------------!
  cangs = 8.3144*tk/p*can2g
  cangs = 1.6*cangs             ! convert CO2 gs into H2O gs
  gsn = cangs
!----------------------------------------------------------------------!
else
  can2a = 0.0
  can2g = 0.0
  canres = 0.
  cangs = 0.0
  gsn = 0.0
  suma = 0.0
endif

call suma_add(suma)


resp_l = canres*3600.0*(24.0-hrs)*12.0/1000000.0
daygpp = can2a*3600.0*hrs*12.0/1000000.0


!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if ((check_closure).and.(pft(co)%sla > 0.0)) then
  total_carbon = ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + &
 ssv(co)%nppstore(1) + ssv(co)%stem%tot(1) + ssv(co)%root%tot(1) + ssv(co)%bio(1) + &
 ssv(co)%bio(2)
  if (abs(total_carbon-old_total_carbon) >  1.0e-3) then
    write(*,*) 'Breach of carbon closure in DOLYDAY:', &
 total_carbon-old_total_carbon,' g/m^2.'
    stop
  endif
endif

end subroutine dailyStep





!**********************************************************************!
!                                                                      !
!                     evapotranspiration :: doly                       !
!                     --------------------------                       !
!                                                                      !
! SUBROUTINE evapotranspiration(t,rh,rn,canga,gsn,hrs,eemm,etmm)       !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Calculate evapotranspiration
!! @details etmm is transpiration in mm ,eemm evaporation in mm
!! @author Mark Lomas EPK
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine evapotranspiration(t,rh,rn,canga,gsn,hrs,eemm,etmm)
!**********************************************************************!
real(dp) :: t,rh,rn,canga,gsn,eemm,etmm,et,hrs,etwt,ee,ra
real(dp) :: Es_t,Delta,s_vp,vp_d
! Monteith Principles 4th Ed. Eq.2.27
real(dp), parameter :: lambda = 2.48*1000 ! Latent heat of vaporization J/g
real(dp), parameter :: Mw = 18.01 ! Molecular weight of water g/mol
real(dp), parameter :: R = 8.31 ! Gas constant J/K/mol
real(dp), parameter :: Es_t1 = 611 ! Saturation vapour pressure at 273 K in Pa
real(dp), parameter :: cp = 1012 ! Specific heat of air in J/Kg/K
real(dp), parameter :: Pa = 101325 ! Absolute Pressure in Pa for 1 atm
real(dp), parameter :: Rsp = 287.058 ! Specific gas constant in J/Kg/K
real(dp), parameter :: gamma = 67.0 ! Psychrometer constant Pa/K
real(dp), parameter :: c_k = 273.15 ! Conversion
!----------------------------------------------------------------------!
! Penman-Monteith equation for evapotranspiration.                     !
!----------------------------------------------------------------------!

! Saturation vapour pressure at temperature t in Pa.Similar to Tetens eq
! Monteith Principles 4th Ed. Eq.2.27
Es_t = Es_t1*exp(17.27*(t+c_k-273.15)/(t+c_k-36))
! Rate of change of saturation vapour pressure at temperature t in Pa/K
! Monteith Principles 4th Ed. Eq.2.28
Delta = lambda*Mw*Es_t/(R*(t+c_k)**2)
! Delta can also be calculated as 48.7*exp(0.0532*t)

!Dry air density in kg/m3 from wiki
ra = Pa/(Rsp*(t+c_k))
! Vapour pressure deficit in Pa
vp_d = (1.0 - rh/100.0)*Es_t

! rn in W/m2 (J/s/m2)
! et here in J/s/m2
if ((ssv(ssp%cohort)%lai%tot(1)>0.1).and.(msv%mv_soil2g>ssp%wilt)) then
  et = (Delta*rn + ra*cp*canga*vp_d)/(Delta + gamma*(1.0 + canga/gsn))
  
  if (et<0.0) then
    et = 0.0
  endif

  ! Units now become g/s/m2  
  etwt = et/lambda
  ! Multiplies by daylight hours in seconds and gets g/m2
  etwt = etwt*3600*hrs
  ! 1g of water is 10^-6 m3 so etwt is g/m2 but also 10^-6 m
  ! or 10^-3 mm
  etmm = etwt/1000.0

else
  et = 0.0
  etwt = 0.0
  etmm = 0.0
endif

! Evaporation,like evapotranspiration but without
! the canopy conductance resistance factor
ee = (Delta*rn + ra*cp*canga*vp_d)/(Delta + gamma)
! Evaporation in mm as above
eemm = (ee*3600.0*hrs)/(lambda*1000.0)

if (ee<0.0) then
  ee = 0.0
  eemm = 0.0
endif

end subroutine evapotranspiration



end module daily_step













