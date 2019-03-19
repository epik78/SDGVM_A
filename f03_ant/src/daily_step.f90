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
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
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
real(dp) :: tmp,prc,hum,cld,ca,maxevap, &
 t,suma,laiinc, &
 soilc,soiln,rh,tk,rd(12),rn,q,qdiff, &
 qdirect,hrs,vpd,canga,lam,rho,s,amx,gsn, &
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
!      abmin = tmp(mnth,day)*1.29772 - 19.5362
if (rh>95.0)  rh=95.0
if (rh<30.0)  rh=30.0


t = tmp
lam = 2500.0 - (2.367*t)
s = 48.7*exp(0.0532*t)
rn = 0.96*(q*1000000.0/4.0 + 208.0 + 6.0*t)
rn = rn*0.52


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
!! @details etmm is transpiration,eemm evaporation
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine evapotranspiration(t,rh,rn,canga,gsn,hrs,eemm,etmm)
!**********************************************************************!
real(dp) :: t,rh,rn,canga,gsn,eemm,etmm,et,svp,vpd,lam,rho,s,gam,hrs, &
 maxevap,etwt,ee

!----------------------------------------------------------------------!
! Penman-Monteith equation for evapotranspiration.                     !
! Units of ET = W/m2  (CHECK) N.B. Conductances in moles.              !
!----------------------------------------------------------------------!
svp = 6.108*exp((17.269*t)/(237.3 + t))*100.0
vpd = (1.0 - rh/100.0)*svp
lam = 2500.0 - 2.367*t
rho = 1288.4 - 4.103*t
s = 48.7*exp(0.0532*t)
gam = 101325.0*1.012/(0.622*lam)

! Transpiration since it uses the canopy conductance (gsn) Eq.38
! calculated in doly
if ((ssv(ssp%cohort)%lai%tot(1)>0.1).and.(msv%mv_soil2g>ssp%wilt)) then
  et = (s*rn + rho*1.012*canga*vpd)/(s + gam*(1.0 + canga/gsn))*tgp%p_et
 
  if (et<0.0) then
    et=0.0
  endif

  etwt = (et*3600.0*hrs)/lam
  etmm = etwt/1000.0
  if (etmm>0.0) then
    if (etmm/hrs>maxevap)  maxevap = etmm/hrs
  endif
else
  et = 0.0
  etwt = 0.0
  etmm = 0.0
endif

! Evaporation,it is similar with evapotranspiration but without
! the canopy conductance factor
ee = (s*rn + rho*1.012*canga*vpd)/(s + gam)
eemm = (ee*3600.0*hrs)/(lam*1000.0)

if (ee<0.0) then
  ee=0.0
  eemm=0.0
endif

end subroutine evapotranspiration



end module daily_step













