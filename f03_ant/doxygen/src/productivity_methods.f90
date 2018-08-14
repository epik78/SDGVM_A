module productivity_methods

use real_precision
use dims
use sunshade
use system_state
use pft_parameters
use site_parameters
use tuning_parameters
use input_methods
use luna_methods
use productivity_luna_methods
use func
use light_methods

implicit none

contains

!**********************************************************************!
!                                                                      !
!                   nppcalc :: productivity_methods                    !
!                   -------------------------------                    !
!                                                                      !
! SUBROUTINE NPPCALC(npp_eff,c3,maxc,soilc,soiln,minn,soil2g,wtwp,     !
! wtfc,rd,rlai,t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma,amx,!
! amax,gsum,hrs,canga,p,day,nleaf_sum,fpr,ft)                          !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine nppcalc(npp_eff,c3,maxc,soilc,soiln,minn,soil2g,wtwp, &
 wtfc,rd,rlai,t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma,amx, &
 amax,gsum,hrs,canga,p,nleaf_sum,fpr,tleaf_n,tleaf_p,ce_light,ce_ci,ce_t, &
 ce_maxlight,ce_ga,ce_rh,ft,par_loops,cld,lat,year,mnth,day,thty_dys,gs_type, &
 swr)
!**********************************************************************!
real(dp) :: suma,sumd,rlai,soilc,soil2g,wtfc,nmult,npp_eff,soiln,y1,y0,&
 x1,x0, mmult,minn,nup,t,tk,up,kc,ko,tau,nleaf_sum,rem,can(12),&
 vm(12),jmx(12),oi,canres,dresp(12),drespt(12),upt(12),nupw,vmx(12),rh,&
 c1,c2,canga,p,ga,q,qdiff,qdirect,cs,coeff,can2a,can2g,rn,wtwp,&
 kg,maxc,rht,tpav,tppcv,tpgsv,ca,rd(12),tpaj,tppcj,tpgsj,hrs,a(12),gs(12),&
 ci(12),tpac4,tppcc4,tpgsc4,xvmax,asunlit,ashade,pcsunlit,pcshade, &
 gssunlit,gsshade,jsunlit(12),jshade(12),qsunlit(12),qshade(12), &
 fsunlit(12),fshade(12),ax,amax,amx,gsum,lyr,qt,t154,xdvy,px,oitau1, &
 oitau2,gsmin,x,y,dv,kcoiko,qt1,qt2,vmxp1,jp1,fpr,apar,can_clump,cld,lat
! light-limited assimilation rate (j) and irradiance (q) for sunlit
! and shade
! fraction of sunlit and shade
real(dp), dimension(12) :: a_sd,ci_sd,gs_sd,fsunlit_sd,fshade_sd,qsunlit_sd, &
 qshade_sd,jsunlit_sd,jshade_sd,maxlight,vcmax,jmax,nleaf,pleaf,jm

real(dp) :: qdirect_sd,qdiff_sd,q_sd,swr,tleaf_n,tleaf_p

real(dp) :: sd_scale(20),sd_scale2,cos_zen
integer :: calc_zen,year,mnth,thty_dys,subd_par,read_par,gs_type
logical :: jfv
integer :: i,lai,c3,day,ft,luna_calc_days,par_loops,ii
real(dp) :: soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo,&
 albedodiff,albedobeam,can_sum,k
real(dp), dimension(30,12) :: ce_light,ce_ci,ce_ga,ce_maxlight
real(dp), dimension(30) :: ce_t,ce_rh
real(dp) :: max_dpchg,canrd
!topleaf values
real(dp) :: tassim,tgs,tci,p_rd,leaf_nit,ncan,pcan

real(dp) :: ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,ftjva,ftjvb,ftvna,ftvnb,ftjvna,ftjvnb,ftg0,ftg1,env_vcmax,env_jmax
real(dp) :: sla
real(dp), dimension(12) :: pnlc,enzs

integer :: s070607
logical :: gold,hw_j
!----------------------------------------------------------------------!

      if ((day.eq.100).and.(mnth.eq.1).and.(year.eq.1902)) stop
!      print*
!      print'(3i5,2f10.5)',year,mnth,day,soil2g,rlai
!      print*,'age in nppstore',ssv(ssp%cohort)%age

pnlc = 0.0

ftToptV = 0.0
ftHaV = 0.0
ftHdV = 0.0
ftToptJ = 0.0
ftHaJ = 0.0
ftHdJ = 0.0
ftvna = 0.0
ftvnb = 0.0
ftjva = 0.0
ftjvb = 0.0
ftjvna = 0.0
ftjvnb = 0.0
env_vcmax = 0.0
env_jmax = 0.0
enzs = 0.0

p_rd = 2.0e-2

ftg0 = pft(ssp%cohort)%g0
ftg1 = pft(ssp%cohort)%g1
sla = pft(ssp%cohort)%sla

! Use Harley 1992 J to Jmax relationship

gold = .true.

if (inp%run%s070607) then
  s070607 = 1
  hw_j = .false.
else
  s070607 = 0
  hw_j = .true.
endif

apar=0.0

rem = rlai - int(rlai)
lai = int(rlai) + 1
tk = 273.0 + t

suma = 0.0
sumd = 0.0
if(inp%run%vcmax_type.ne.9) then
  vcmax = 0.0
  jmax  = 0.0
endif
leaf_nit  = 0.0
nleaf_sum = 0.0
canres    = 0.0
nleaf     = 0.0

!----------------------------------------------------------------------!
! Set some parameters used in GOUDRIANS LAW routine (quicker to do it
! here).
!----------------------------------------------------------------------!
! soilalbedo: value from Wang
! leafalbedo: leaf transmittance and reflectance are equal. Then albedo
! is the double kbeam: extinction coefficent of the direct beam
! kdiff: extinction coefficient of diffuse radiations
! I found somewhere an approximate relation: kbeam=0.75*kdiff
! ... hum, I don't know if albedodiff is realy different or not...
! so equal
!----------------------------------------------------------------------!
soilalbedo = 0.15
leafalbedo = 0.15
kbeam = 0.5
kdiff = 0.66
m = sqrt(1.0-leafalbedo)
kbeamstar = m*kbeam
canopyalbedo = 2.0*kbeam/(kbeam + kdiff)*(1.0 - m)/(1.0 + m)
albedobeam = canopyalbedo + (soilalbedo - canopyalbedo)*exp(-kbeam*rlai)
albedodiff = albedobeam
!----------------------------------------------------------------------!
qt = 2.3
qt1 = qt**(t/10.0)/qt**(2.5)
qt2 = qt**(30.0/10.0)/qt**(2.5)

suma = 0.0
sumd = 0.0
nleaf_sum = 0.0

! parameter - interval between luna calculations of Vcmax (in days)
luna_calc_days = 10

! below are the kinetic parameters from Farquhar etal 1980
kc  = exp(35.8 - 80.5/(0.00831*tk))
ko  = exp(9.6 - 14.51/(0.00831*tk))*1000.0
tau = exp(-3.949 + 28.99/(0.00831*tk))

! changed by Ghislain 08/10/03      IF (rlai.GT.0.1) THEN
q = qdiff + qdirect

! if there are light and leaves - calculate canopy properties 
if ((rlai>0.1).and.(q>0.0)) then
! Influence of soil water on stomatal conductance (kg).
! From Eq.21 slightly changed
  if (soil2g>wtwp) then
    kg = maxc*((soil2g - wtwp)/(wtfc - wtwp))**tgp%p_kgw
  else
    kg = 0.0
  endif
  if (kg>maxc)  kg = maxc
  vmxp1 = npp_eff*kg**tgp%p_nu3

! Nitrogen uptake.Very different from Eq.31-Eq.36
! Here it is proportional to soilc and kg
! Then multiplied  by soiln and some constants
! Nitrogen uptake CHECK.
! nup = 120.0*(exp(-0.00008*soilc))

  if (inp%run%ncalc_type==0) then
    nupw = tgp%p_nu1*(exp(tgp%p_nu2*soilc))*kg**tgp%p_nu3
  elseif (inp%run%ncalc_type==1) then
    nupw = tgp%p_nu1*(exp(tgp%p_nu2*soilc))
  else
    write(*,*) 'ncalc_type must be 0 or 1. ncalc_type = ',inp%run%ncalc_type
    stop
  endif
  if (nupw<0.0) nupw = 0.0

! Nitrogen multiplier.
  IF(pft(ft)%phen.NE.3) THEN
    nmult = soiln*tgp%p_nu4
    if (nmult>=1.0)  nmult = 1.0
  ELSE
! nmult = soiln*tgp%p_nu4
    nmult=pft(ft)%fert(6)+pft(ft)%fert(1)/pft(ft)%fert(4)
  ENDIF

! All the lines that lead to the calculation of mmult are not
! needed since its not used,variable nup is overwritten in the
! next line
  y1 = 1.3
  y0 = 1.0
  x1 = 50.0
  x0 = 10.0
  mmult = minn*(y1 - y0)/(x1 - x0) + (y0*x1 - x0*y1)/(x1 - x0)
  if (mmult>y1)  mmult = y1
  if (mmult<y0)  mmult = y0
  
  nup = nupw*nmult*mmult
  nup = nupw*nmult
  up = nup
  if (up<0.0)  up = 0.0

!----------------------------------------------------------------------!
! canopy N scaling
!----------------------------------------------------------------------!
  k = 0.5
  if (inp%run%cstype.eq.1)  k = 0.5*pft(ssp%cohort)%can_clump
! Total the Beer's Law values for each lai up to the current one.      !
  if (inp%run%s070607) then
    coeff =-0.5
    can_sum = 0.0 
    do i=1,lai-1
      k = i - 1
      can_sum = can_sum + exp(coeff*real(k))
    enddo
    can_sum = can_sum + exp(coeff*real(lai - 1))*rem
  else
! Total the Beer's Law values for canopy
    can_sum = 0.0
    do i=1,lai-1
      if(inp%run%cstype.lt.2) then
        can_sum = can_sum + exp(-k*real(i))
      elseif (inp%run%cstype.eq.2) then
! use light proportion to scale
        can_sum = can_sum + sum(ce_light(:,i))
      else
        write(*,*) 'cstype ',inp%run%cstype,' undefined. set to a value<3 in &
&<input.dat> or define additional canopy scaling method'
        stop
      endif
    enddo
    if(inp%run%cstype.lt.2) can_sum = can_sum + exp(-k*real(lai))*rem
    if(inp%run%cstype.eq.2) can_sum = can_sum + sum(ce_light(:,lai))*rem
  endif

!calculate maximum daily change in vcmax for LUNA model after Ali, Xu, et al 2015
  if (inp%run%vcmax_type.eq.9) then
    max_dpchg = max_daily_pchg(sum(ce_t(:))/30.0)  
  endif

!----------------------------------------------------------------------!
! Canopy scaling of variables that do not vary with light
!----------------------------------------------------------------------!
  do i=1,lai
    lyr = real(i)

    if (inp%run%s070607) then
      can(i) = exp(coeff*(lyr - 1.0))/can_sum
    else
! Proportion of canopy N in the LAI layer 'lyr'
! - this is supposed to be proportional to the light in layer 'lyr' 
! but this is inconsistent with the rad scheme used by SDGVM
      if(inp%run%cstype.lt.2) then 
        can(i) = exp(-k*(lyr))/can_sum
      elseif(inp%run%cstype.eq.2) then
! use light proportion to scale
        can(i) = sum(ce_light(:,i))/can_sum
      endif
    endif
! Leaf nitrogen
!    nleaf(i) = (upt(i))*tgp%p_nleaf
    if((i.eq.lai).and.(.not.(inp%run%s070607))) can(i) = can(i)*rem

    upt(i) = up*can(i)

!----------------------------------------------------------------------!
! calculate leaf N in canopy layer i
!----------------------------------------------------------------------!
    IF (inp%run%ncalc_type.le.1) THEN 
! default topleaf N 
      nleaf(i) = up*can(i)*tgp%p_nleaf
    ELSEIF (inp%run%ncalc_type.eq.2) THEN
! based on trait regressions/specified top leaf N
! calculate total canopy n
       ncan     = tleaf_n*can_sum
       nleaf(i) = can(i)*ncan
    ELSE
      write(*,*) 'ncalc_type ',inp%run%ncalc_type,' undefined. set to a value &
     &<3 in <input.dat> or define your own N calculation method'
       STOP
    ENDIF

! this scales the bottom layer from a fractional layer to a full layer
! - unless total LAI is <1, in which case leaf N is assumed for a full leaf layer 
! calculations below this assume a full leaf layer and then scale all fluxes in lowest lai layer by rem 
    if (.not.(inp%run%s070607)) then 
      if ((i.eq.lai).and.(lai.ne.1)) nleaf(i) = nleaf(i)/rem
    endif

!----------------------------------------------------------------------!
! Dark respiration proportional to nitrogen uptake
! Different from Eq.27,not controlled by temperature
!----------------------------------------------------------------------!
    dresp(i) = nleaf(i)*tgp%p_dresp
    drespt(i) = dresp(i)
!----------------------------------------------------------------------!

! If soil water not completely limiting calculate photosynthesis
    if (kg>1.0e-10) then

!----------------------------------------------------------------------!
! calculate Vcmax@25oC (vm) 
!----------------------------------------------------------------------!
      IF(inp%run%vcmax_type.eq.0) THEN
! default 070607 SDGVM
        vm(i) = nleaf(i)*tgp%p_vm         

      ELSEIF(inp%run%vcmax_type.eq.1) THEN
! from Walker et al - N only
        vm(i) = exp(3.712+0.65*log(nleaf(i)))

! from Walker et al - N & P

        IF(inp%run%soilp_map.eq.1) THEN
! assume same change in leaf p through the canopy as leaf n
          pcan = tleaf_p*can_sum
          pleaf(i) = can(i)*pcan
          vm(i) = exp(3.946 + 0.921*log(nleaf(i)) + & 
 0.121*log(pleaf(i))+0.282*log(nleaf(i))*log(pleaf(i)))
        ENDIF

      ELSEIF ((inp%run%vcmax_type.eq.2).OR.(inp%run%vcmax_type.eq.3).or.(inp%run%vcmax_type.eq.5).or.(inp%run%vcmax_type.eq.6)) THEN
! vcmax based on environment 
        vm(i) = env_vcmax * can(i)/can(1)

      ELSEIF (inp%run%vcmax_type.eq.4) THEN
! specified as a PFT parameter
        vm(i) = ftvna + ftvnb*nleaf(i)
            
! if specified as a constant, i.e. ftvnb = 0 then scale using can
        IF(ftvnb.lt.1d-5) vm(i) = vm(i)*can(i)/can(1)

      ELSEIF (inp%run%vcmax_type.eq.7) THEN
! after Maire et al 2012
            
!       if(i.eq.1) print'(30(f6.1,1x))', ce_light(:,1) *1d6
!       if(i.eq.1) print'(f6.1)',   sum(ce_light(:,1))/30.d0 *1d6
           
! Jmax is simulated as a function of Vcmax 
        jfv    = .TRUE.
! calculate Vcmax at the mean temperature of the last month in mol m-2s-1
        vm(i) = brent_solver(0,1.0e-6_dp,5.0e-4_dp,oi,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
 sum(ce_t(:))/30.0,0.0_dp,1,0.0_dp,0.0_dp,0.0_dp,i, &
 sum(ce_ci(:,i))/30.0,sum(ce_light(:,i))/30.0,inp%run%ttype, &
 ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.0,jfv)

! convert back to value at 25oC 
! invert temp correction scalar to get values at 25oC
        vm(i) = vm(i)/T_SCALAR(sum(ce_t(:))/30.0,'v',inp%run%ttype, &
 ftToptV,ftHaV,ftHdV,sum(ce_t(:))/30.0,jfv) 
        jm(i) = jm(i)/T_SCALAR(sum(ce_t(:))/30.0,'j',inp%run%ttype, &
 ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.0,jfv)

! convert to umol m-2s-1
        vm(i) = vm(i)*1e6

      ELSEIF (inp%run%vcmax_type.eq.8) THEN
! specified as a PFT parameter but jmax specified as the Walker function of Vcmax
        vm(i) = ftvna + ftvnb*nleaf(i)
             
! if specified as a constant, i.e. ftvnb = 0 then scale using can
        IF (ftvnb.lt.1e-5) vm(i) = vm(i)*can(i)/can(1)

      ELSEIF (inp%run%vcmax_type.eq.9) THEN
! use LUNA model after Ali, Xu, et al 2015

        vm(i) = vcmax(i)
        jm(i) = jmax(i)

        IF (mod(day,luna_calc_days).eq.0) THEN
          call LUNA(i,vm(i),jm(i),PNlc(i),enzs(i),sum(ce_ga(:,i))/30.0, &
 sla,nleaf(i),sum(ce_light(:,i))/30.0,sum(ce_maxlight(:,i))/30.0, &
 ca,oi,hrs,sum(ce_rh(:))/30.0,sum(ce_t(:))/30.0, &
 gs_type,ftg0,ftg1,max_dpchg*luna_calc_days,inp%run%ttype,ftToptV,ftHaV, &
 ftHdV,ftToptJ,ftHaJ,ftHdJ)  
        ENDIF

      ELSE
        write(*,*) 'vcmax_type:',inp%run%vcmax_type,'undefined. set to a value &
 &1-9 in <input.dat> or define your own vcmax calculation method'
        STOP
      ENDIF

!----------------------------------------------------------------------!
! calculate Jmax@25oC (jm) 
!----------------------------------------------------------------------!   
      IF ((inp%run%vcmax_type.eq.0).or.(inp%run%vcmax_type.eq.5).or.(inp%run%vcmax_type.eq.6)) THEN
! 070607 default - Wullschleger 1993 
        jm(i)  = (29.1 + 1.64*vm(i))
        jfv    = .TRUE.
      ELSEIF((inp%run%vcmax_type.eq.1).or.(inp%run%vcmax_type.eq.7).or.(inp%run%vcmax_type.eq.8)) THEN
! Walker -  only vcmax
        jm(i)  = exp(1.0+0.890*log(vm(i)))
        jfv    = .TRUE.
      ELSEIF((inp%run%vcmax_type.eq.2).OR.(inp%run%vcmax_type.eq.3)) THEN
! from van Bodegom for TERRABITES project
! - not a function of vcmax 
        jm(i)  = env_jmax * can(i)/can(1)
        jfv    = .FALSE.
      ELSEIF(inp%run%vcmax_type.eq.4) THEN
! specified as a PFT parameter
        jm(i)  = ftjva + ftjvb*vm(i)
        jfv    = .TRUE.
      ELSEIF(inp%run%vcmax_type.eq.9) THEN
! do nothing - jmax already defined above            
        jfv    = .FALSE.
      ELSE
        write(*,*) 'vcmax_type ',inp%run%vcmax_type,' undefined. set to a value<1-9 in <input.dat>'
        STOP
      ENDIF

!----------------------------------------------------------------------!
! Always set C4 vcmax and PEPC (hijacking the jmax variable) from van Bodegom trait environment
! relationships or mean. 
! currently assumes same canopy and temp scaling and age reduction as C3 
!----------------------------------------------------------------------!
      IF (c3.NE.1) THEN
        vm(i) = env_vcmax * can(i)/can(1)
        jm(i) = env_jmax  * can(i)/can(1)
        jfv   = .FALSE.
      ENDIF

! temperature scale Vcmax & Jmax
!     if(inp%run%vcmax_type.le.6) then
      vmx(i) = vm(i)*T_SCALAR(t,'v',inp%run%ttype,ftToptV,ftHaV,ftHdV,sum(ce_t(:))/30.0,jfv)
      jmx(i) = jm(i)*T_SCALAR(t,'j',inp%run%ttype,ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.0,jfv)

      vmx(i) = vm(i)
      jmx(i) = jm(i)

! leaf age scale Vcmax & Jmax
      vmx(i) = vmx(i)*npp_eff
      jmx(i) = jmx(i)*npp_eff

! water limitation scale Vcmax & Jmax
      vmx(i) = vmx(i)*kg**tgp%p_nu3
      jmx(i) = jmx(i)*kg**tgp%p_nu3

      if (inp%run%s070607) then
        vm(i) = nleaf(i)*tgp%p_vm
        vmx(i)=vm(i)*qt1
        if (t>30.0) vmx(i)=qt2
        vmx(i) = vmx(i)*vmxp1

! jmx is light-saturated rate of electron transport Eq.13
        jmx(i) = (29.1 + 1.64*vmx(i))
! vmx is the maximum rate of carboxulation Vmax obtained from
! leaf nitrogen and temperature
      endif

!convert from umol to mol
      vmx(i) = vmx(i) * 1.0e-6
      jmx(i) = jmx(i) * 1.0e-6

! dark respiration (daytime) mol m-2 s-1
! assumes rd is proportional to soil water lim
      if (inp%run%s070607) then
        rd(:) = vmx(i)*0.02
      else
! SDGVM previously calculated rd but wasn't layer specific
        rd(i) = vmx(i)*p_rd
      endif

! if subd_par = 0 then calculate light here  
      if (.not.inp%run%subdaily) then
! calculate incident light in canopy layer
!        CALL GOUDRIAANSLAW_ANT(lyr-0.5,rlai,qdirect,qdiff, &
! fsunlit(i),qsunlit(i),fshade(i),qshade(i),can_clump,cos_zen, &
! inp%run%s070607)

        call GOUDRIAANSLAW2(lyr-0.5,rlai,qdirect,qdiff,fsunlit(i),&
 qsunlit(i),fshade(i),qshade(i),soilalbedo,leafalbedo,kbeam,&
 kdiff,m,kbeamstar,canopyalbedo,albedobeam,albedodiff)

! calculate the sunlit and shaded light limited electron transport rate
        call JCALC(qsunlit(i),jmx(i),hw_j,jsunlit(i))
        call JCALC(qshade(i),jmx(i),hw_j,jshade(i))

! switch 1 endif
      endif

!----------------------------------------------------------------------!
! water limitation IF
!----------------------------------------------------------------------!
    else

      if (subd_par.eq.0) then
        CALL GOUDRIAANSLAW_ANT(lyr-0.5,rlai,qdirect,qdiff, &
 fsunlit(i),qsunlit(i), fshade(i),qshade(i),can_clump,cos_zen,inp%run%s070607)

        vmx(i)     = 0.0
        jsunlit(i) = 0.0
        jshade(i)  = 0.0
      endif

! LUNA model has a term that reduces Vcmax during drought and periods of stress 
      if (inp%run%vcmax_type.eq.9) then
        vm(i) = vcmax(i)
        jm(i) = jmax(i)

        call LUNA_nogrowth(max_dpchg,vm(i),jm(i),enzs(i))  
! temperature scale Vcmax & Jmax
        vmx(i) = vm(i)*T_SCALAR(t,'v',inp%run%ttype,ftToptV,ftHaV,ftHdV,sum(ce_t(:))/30.0,.FALSE.)
        jmx(i) = jm(i)*T_SCALAR(t,'j',inp%run%ttype,ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.0,.FALSE.)

! leaf age scale Vcmax & Jmax
        vmx(i) = vmx(i)*npp_eff
        jmx(i) = jmx(i)*npp_eff
          
! water limitation scale Vcmax & Jmax
        vmx(i) = vmx(i)*kg**tgp%p_nu3
        jmx(i) = jmx(i)*kg**tgp%p_nu3

! convert from umol to mol
        vmx(i) = vmx(i) * 1.0e-6
        jmx(i) = jmx(i) * 1.0e-6

      endif

! water limitation IF
    endif

! leaf respiration - nighttime 
    if ((inp%run%vcmax_type.eq.2).or.(inp%run%vcmax_type.eq.3).or.(inp%run%vcmax_type.eq.5).or. &
 (inp%run%vcmax_type.eq.6).or.(inp%run%vcmax_type.eq.7) ) then
! these vcmax types assume no leaf N and that day AND night leaf resp = f(vcmax)
      nleaf(i) = 0.0
      nleaf_sum= 0.0
      dresp(i) = vmx(i)*p_rd ! * nightr:dayr scalar
    else
      dresp(i)  = nleaf(i)*tgp%p_dresp
    endif 
    drespt(i) = dresp(i)         ! I'm not sure what drespt does - nothing I think

! sum to get canopy N & Respiration
    IF (i.LT.lai) THEN
      nleaf_sum = nleaf_sum + nleaf(i)
      canres    = canres + dresp(i)
    ELSE
      nleaf_sum = nleaf_sum + nleaf(i)*rem
      canres    = canres + dresp(i)*rem
    ENDIF

! initial LAI loop
  enddo

! else there is no light or leaves
else

  if (inp%run%vcmax_type.eq.9) then
! LUNA model retains the final value of vcmax pre total leaf-loss to initialise the following year
    vm(:)     = vcmax(:) 
    jm(:)     = jmax(:)
  else
    vm(:)     = 0.0
    jm(:)     = 0.0
  endif

  can(:)     = 0.0
  upt(:)     = 0.0
  nleaf(:)   = 0.0
  dresp(:)   = 0.0
  drespt(:)  = 0.0
  jmx(:)     = 0.0
  vmx(:)     = 0.0
  jshade(:)  = 0.0
  jsunlit(:) = 0.0
  fsunlit(:) = 0.0
  fshade(:)  = 0.0
  qsunlit(:) = 0.0
  qshade(:)  = 0.0
        
  canres   = 0.0
  !vcmax    = 0.0
  !jmax     = 0.0
  leaf_nit = 0.0        
 
! end light and leaves loop
endif

!----------------------------------------------------------------------!
! Assimilation calculations using subroutines ASSVMAX and ASSJ.        !
!----------------------------------------------------------------------!
can2a = 0.0
can2g = 0.0

!----------------------------------------------------------------------!
! Calculation of constant multipliers used in assimilation calcs.      !
! Used to reduce the number of floating point operations and so tov    !
! speed up the model.                                                  !
!----------------------------------------------------------------------!
t154 = 1.54*t
px = 1.0e-6*p
oitau1 = oi/tau
oitau2 = 0.5*oi/tau
gsmin = 0.005
x = 3.0 + 0.2*t
if (x<0.25) x = 0.25
y = 1.5
dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)
xdvy = x/(1.0 + dv/y)
kcoiko = kc*(1.0 + oi/ko)
!----------------------------------------------------------------------!

! zero layer arrays
a = 0.0
ci = 0.0
gs = 0.0

if (rlai>0.1) then

!  if (rlai>0.1) then
! This should transform the boundary conductance to molar 
! (Ant) leaf boundary layer conductance 
! this needs attention, this is not right to use canga divided by lai
    ga = canga/(8.3144*tk/p)/rlai

!------------------------------------------------------------------------!
!subd_par if statement mean daily PAR or downscaled sub-daily PAR 
! - if a daily mean PAR is used (SDGVM default) use above calculated PAR values
!------------------------------------------------------------------------!
  if (.not.(inp%run%subdaily)) then

    do i=1,lai

      ssp%lai = i
      dresp(i) = drespt(i)

! assimilation LAI loop - fshade, fsunlit, jshade, jsunlit have been calculated above 
      call assimilation_calc(c3,t,rn,soil2g,wtwp,fshade(i),fsunlit(i),qshade(i),qsunlit(i), &
 jshade(i),jsunlit(i),canga,tk,p, &
 rlai,rh,oi,ca,vmx(i),kg,ko,kc,tau,rd(i),ga,t154,xdvy,dv,px,oitau1,oitau2, &
 gsmin,kcoiko,upt(i),a(i),ci(i),gs(i),gs_type,1)

    enddo

  else
!------------------------------------------------------------------------!
! scaling factor for integration of subdaily variability 
!------------------------------------------------------------------------!
    sd_scale(:) = 1.0
    sd_scale(1) = 2.0
    sd_scale(par_loops+1) = 2.0
    sd_scale2   = 1.0 / real(par_loops)

! zero arrays
    a  = 0.0
    ci = 0.0
    gs = 0.0
    fsunlit = 0.0
    fshade  = 0.0
    qsunlit = 0.0
    qshade  = 0.0
    qdirect = 0.0
    qdiff   = 0.0        

    a_sd  = 0.0
    ci_sd = 0.0
    gs_sd = 0.0
    fsunlit_sd = 0.0
    fshade_sd  = 0.0
    qsunlit_sd = 0.0
    qshade_sd  = 0.0
    qdirect_sd = 0.0 
    qdiff_sd   = 0.0

    do ii=1,par_loops+1
!    do ii=par_loops,par_loops

      call pfd_ant(lat,no_day(year,mnth,day,thty_dys),hrs,cld, &
 qdirect_sd,qdiff_sd,q_sd,swr,inp%run%read_par,inp%run%subdaily, &
 (ii-1),par_loops,inp%run%calc_zen,cos_zen)

! At dawn assume diffuse light is 5 umol/m2/s and cos_zen a fraction above zero
      if (qdirect_sd + qdiff_sd<5.0e-6)  qdiff_sd = 5.0e-6  
      if (cos_zen.le.0.0)  cos_zen  = cos(3.141/2)

      do i=1,lai
        ssp%lai = i

! In the average PAR version of the model this is contained within a water 
! limitation if statement
        call goudriaanslaw_ant(real(i,dp)-0.5_dp,rlai,qdirect_sd,qdiff_sd, &
 fsunlit_sd(i),qsunlit_sd(i),fshade_sd(i),qshade_sd(i),pft(ssp%cohort)%can_clump,cos_zen, &
 inp%run%s070607)

        call jcalc(qsunlit_sd(i),jmx(i),inp%run%s070607,jsunlit_sd(i))
        call jcalc(qshade_sd(i),jmx(i),inp%run%s070607,jshade_sd(i))

        call assimilation_calc(c3,t,rn,soil2g,wtwp,fshade_sd(i),fsunlit_sd(i), &
 qshade_sd(i),qsunlit_sd(i),jshade_sd(i),jsunlit_sd(i),canga,tk,p, &
 rlai,rh,oi,ca,vmx(i),kg,ko,kc,tau,rd(i),ga,t154,xdvy,dv,px,oitau1,oitau2, &
 gsmin,kcoiko,upt(i),a_sd(i),ci_sd(i),gs_sd(i),gs_type,ii)

      enddo

! Numerically integrate canopy values to get daily mean 
! (mean is necessary because per second values are scaled to daytime values in doly)   
      a(:)       =  a(:)       + a_sd(:)       / sd_scale(ii) 
      ci(:)      =  ci(:)      + ci_sd(:)      / sd_scale(ii)  
      gs(:)      =  gs(:)      + gs_sd(:)      / sd_scale(ii)  
      fsunlit(:) =  fsunlit(:) + fsunlit_sd(:) / sd_scale(ii)  
      fshade(:)  =  fshade(:)  + fshade_sd(:)  / sd_scale(ii)  
      qsunlit(:) =  qsunlit(:) + qsunlit_sd(:) / sd_scale(ii)  
      qshade(:)  =  qshade(:)  + qshade_sd(:)  / sd_scale(ii)  
      qdirect    =  qdirect    + qdirect_sd    / sd_scale(ii)  
      qdiff      =  qdiff      + qdiff_sd      / sd_scale(ii)  

      if (ii==1) maxlight(:) = fsunlit_sd(:)*qsunlit_sd(:) + fshade_sd(:)*qshade_sd(:)
!------------------------------------------------------------------------!
! End of sub-daily loop
!------------------------------------------------------------------------!
    enddo

    a       =  a        * sd_scale2  
    ci      =  ci       * sd_scale2  
    gs      =  gs       * sd_scale2  
    fsunlit =  fsunlit  * sd_scale2
    fshade  =  fshade   * sd_scale2 
    qsunlit =  qsunlit  * sd_scale2
    qshade  =  qshade   * sd_scale2
    qdirect =  qdirect  * sd_scale2
    qdiff   =  qdiff    * sd_scale2 

  endif

! calculate environment variables for Maire & LUNA vcmax calc
  ce_t(1:29)  = ce_t(2:30)
  ce_rh(1:29) = ce_rh(2:30)
  do i=1,lai
    ce_ci(1:29,i)       = ce_ci(2:30,i) 
    ce_ga(1:29,i)       = ce_ga(2:30,i) 
    ce_light(1:29,i)    = ce_light(2:30,i)
    ce_maxlight(1:29,i) = ce_maxlight(2:30,i)
  enddo

  ce_t(30)       = t
  ce_rh(30)      = rh
  ce_ci(30,:)    = ci(:)
  ce_ga(30,:)    = ga
  ce_light(30,:) = fsunlit(:)*qsunlit(:)+fshade(:)*qshade(:)    
  ce_maxlight(30,:) = maxlight 

!----------------------------------------------------------------------!
! End of LAI loop.                                                     !
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
! End of 't & rn > 0' 'IF' statement.                                  !
!----------------------------------------------------------------------!
else
  a = 0.0
  gs = 0.0
  dresp = 0.0
  ci = ca
endif

!----------------------------------------------------------------------!
! Sum canopy layers to get canopy integrated values
!----------------------------------------------------------------------!
can2a = 0.0
can2g = 0.0
canrd = 0.0
apar  = 0.0
DO i=1,lai
  IF (i.LT.lai) THEN
    can2a = can2a + a(i)
    can2g = can2g + gs(i)
    canrd = canrd + rd(i)
    apar  = apar  + fsunlit(i)*qsunlit(i)+fshade(i)*qshade(i)
  ELSE
    can2a = can2a + a(i)*rem
    can2g = can2g + gs(i)*rem
    canrd = canrd + rd(i)*rem
    apar  = apar + rem*(fsunlit(i)*qsunlit(i)+fshade(i)*qshade(i))
  ENDIF
ENDDO

! calculate fapar
fpr    = apar/(qdiff+qdirect)

! topleaf values
tassim =  a(1)
tgs    = gs(1)
tci    = ci(1)

!----------------------------------------------------------------------!
! Cumulate daily assimilation of the last lai layer, 'suma'.           !
!----------------------------------------------------------------------!
IF (lai.GT.1) THEN
  suma = suma + (rem*a(lai) + (1.0 - rem)*a(lai - 1))* &
 3600.0*hrs - (rem*dresp(lai) + (1.0 - rem)*dresp(lai - 1))*3600.0*(24.0 - hrs)
  sumd = sumd + (rem*dresp(lai) + (1.0 - rem)*dresp(lai - 1))*3600.0*(24.0 - hrs)
ELSE
! suma = suma + rem*a(lai)*3600.0*hrs - rem*dresp(lai)* &
!3600.0*(24.0 - hrs)
  sumd = sumd + rem*dresp(lai)*3600.0*(24.0 - hrs)
  suma = 0.0
ENDIF


if (lai>1) then
  gsum = gs(lai - 1)*(1.0 - rem) + gs(lai)*rem
else
  gsum = gs(1)
endif

if ((rlai>0.0).and.(t>0.0).and.(rn>0.0)) then
  if (mod(day-1,1)==0) then
! CALL XAMAX(ax,t,oi,xvmx,rd,xj,xdresp,ci(1))
! Temporary change by Ghislain
     ax=0
! end if change
    if (rlai<-1.0) then
      if (ax*rem>amax)  then
        amax = a(1)*rem
      endif
      if (ax*rem>amx)  amx = ax*rem
    else
      if (ax>amax)  then
        amax = a(1)
      endif
      if (ax>amx)  amx = ax
    endif
  endif
endif

suma = suma*1.0e-6


end subroutine nppcalc





!**********************************************************************!
!                                                                      !
!                   assj2 :: productivity_methods                      !
!                   -----------------------------                      !
!                                                                      !
! ASSJ computes the assimilation rate by forming a cubic equation      !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,j,g0,g1,rh,kg,tau,rd                                    !
!     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  !
!                                                                      !
! SUBROUTINE assj2(a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck)        !
!                                                                      !
!**********************************************************************!
subroutine assj2(a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck)
!**********************************************************************!
real(dp) :: a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck,r2q3,arg1, &
 ax,bx,cx,dx,ex,fx,gx,hx,p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc
!----------------------------------------------------------------------!

pi = atan(1.0)*4.0
sc = 1000000.0

!----------------------------------------------------------------------!
! 'a2,b2,c2,d2' are such that 'w=(a2a + b2)/(c2a + d2)'.               !
!----------------------------------------------------------------------!
ax = j*tau*pa*(kg*g1*rh - 160.0)
bx = j*tau*pa*kg*g0*pa
cx = 4.0*pa*kg*tau*g1*rh - 640.0*pa*tau + 4.0*po*kg*g1*rh
dx = 4.0*pa**2.0*kg*tau*g0 + 4.0*po*kg*g0*pa

!----------------------------------------------------------------------!
! 'ex,fx,gx,hx' are such that '0.5po/(taupc)-1=(exa+fx)/(gxa+hx)'.     *
!----------------------------------------------------------------------!
ex =-tau*pa*kg*g1*rh + 160.0*tau*pa + 0.5*po*kg*g1*rh
fx =-tau*pa**2.0*kg*g0 + 0.5*po*kg*g0*pa
gx = tau*pa*(kg*g1*rh - 160.0)
hx = tau*pa**2.0*kg*g0

!----------------------------------------------------------------------!
! Cubic coefficients 'p,q,r'.                                          !
!----------------------------------------------------------------------!
p = sc*rd + (sc*ax*ex + cx*hx + dx*gx)/(cx*gx)
q = (sc*rd*cx*hx + sc*rd*dx*gx + sc*ax*fx + sc*bx*ex + dx*hx)/(cx*gx)
r = (sc*rd*dx*hx + sc*bx*fx)/(cx*gx)

!----------------------------------------------------------------------!
! Sove the cubic equaiton to find 'a'.                                 !
!----------------------------------------------------------------------!
Qt = (p**2 - 3.0*q)/9.0
Rt = (2.0*p**3 - 9.0*p*q + 27.0*r)/54.0

r2q3 = Rt**2 - Qt**3

if (r2q3<0.0) then
  arg1 = Rt/Qt**1.5
  if (arg1>1.0)  arg1 = 1.0
  if (arg1<-1.0)  arg1 = -1.0
  th = acos(arg1)
  a1 =-2.0*Qt**0.5*cos(th/3.0) - p/3.0
  a2 =-2.0*Qt**0.5*cos((th + 2.0*pi)/3.0) - p/3.0
  a3 =-2.0*Qt**0.5*cos((th + 4.0*pi)/3.0) - p/3.0
  a = a1
  if (a2>a)  a = a2
  if (a3>a)  a = a3
else
  at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5)**(1.0/3.0)
  if (abs(at)>1e-6) then
    bt = Qt/at
  else
    bt = 0.0
  endif
  a1 = at + bt - p/3.0
!        a2 =-0.5*(at + bt) - p/3.0 + i*3.0**0.5*(at - bt)/2.0
!        a3 =-0.5*(at + bt) - p/3.0 - i*3.0**0.5*(at - bt)/2.0
  a = a1
endif

!----------------------------------------------------------------------!
! Compute 'gs,pc,w' corresponding to 'a'.                              !
!----------------------------------------------------------------------!
gs = (g0 + g1*a*rh/pa)*kg
pc = pa - a*160.0/gs
w = j*pc/4.0/(pc + po/tau)
acheck = (w*(1.0 - 0.5*po/(tau*pc)) - rd)*sc

!**********************************************************************!
end subroutine assj2
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   xamax :: productivity_methods                      !
!                   -----------------------------                      !
!                                                                      !
! This subroutine calculates the value of 'amax' at the maximum        !
! iradience as opposed to the average irradiance in the main program.  !
!                                                                      !
! SUBROUTINE xamax(a,xt,oi,vmx,rd,j,dresp,pc)                          !
!                                                                      !
!**********************************************************************!
subroutine xamax(a,xt,oi,vmx,rd,j,dresp,pc)
!**********************************************************************!
real(dp) :: a,xt,oi,vmx,rd,j,dresp,pc,t,tk,c1,c2,kc,ko,tau,tpav,tpwv,&
 tpaj,tpwj,sc
!       tppcj,tppcv,tpgsv,tpgsj,acheck
!----------------------------------------------------------------------!

t = 2.64 + 1.14*xt
tk = 273.0 + t
sc = 1000000.0

c1 = 142.4 - 4.8*t
if (c1>80.0)  c1 = 80.0
if (c1<8.0)  c1 = 8.0
c2 = 12.7 - 0.207*t
if (c2>10.0)  c2 = 10.0
if (c2<6.9)  c2 = 6.9
kc = exp(35.8 - 80.5/(0.00831*tk))
ko = exp(9.6 - 14.51/(0.00831*tk))*1000.0
tau = exp(-3.949 + 28.99/(0.00831*tk))

!      rht = rh
!      IF (kg*c2*rh.LT.170.0)  rht = 170.0/(kg*c2)
!      CALL ASSVMAX(tpav,tpwv,tppcv,tpgsv,oi,ca,vmx,c1,c2,
!     &rht,ko,kc,kg,tau,rd,acheck)
!      CALL ASSJ(tpaj,tpwj,tppcj,tpgsj,oi,ca,j,c1,c2,rht,
!     &kg,tau,rd,acheck)

tpwv = vmx*pc/(pc + kc*(1.0 + oi/ko))
tpav = (tpwv*(1.0 - 0.5*oi/(tau*pc)) - rd)*sc

tpwj = j*pc/(4.0*(pc + oi/tau))
tpaj = (tpwj*(1.0 - 0.5*oi/(tau*pc)) - rd)*sc

if (tpav<tpaj) then
  A = tpav
else
  A = tpaj
endif
if (A<-dresp)  A =-dresp

!**********************************************************************!
end subroutine xamax
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assc42 :: productivity_methods                     !
!                   ------------------------------                     !
!                                                                      !
! ASS4 computes the assimilation rate by forming a cubic equation      !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. This is then       !
! compared with assimilation limited to light dependant and vmax       !
! dependant assimilation and the minimum one is taken. The variables   !
! here correspond to the Doly paper, in the main program in the        !
! folowing mannor.                                                     !
!        po,pa,j,g0,g1,rh,kg,tau,rd                                    !
!     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  !
!                                                                      !
! SUBROUTINE assc42(a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,      !
! acheck,up)                                                           !
!                                                                      !
!**********************************************************************!
subroutine assc42(a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,acheck,&
 up)
!**********************************************************************!
real(dp) :: a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,acheck,up, &
 p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc,kp,kmq,pg,l,vmq,absx, &
 qrub,f,a0,b0,b1,c1,d1,e1,f1,b2,c2,d2,b3,c3,d3,wq,wv, &
 a4,b4,c4,d4,f2,aq,av,a2t,r2q3,arg1
!----------------------------------------------------------------------!

pi = atan(1.0)*4.0
sc = 1000000.0

kp = 0.7
kmq = 2.0
pg = 101325.0
l = 0.000005
vmq = 2.0
absx = 0.8
qrub = 0.11
f = 0.6

vmax = up*360.0*vmq**((t - 25.0)/10.0)/(360.0 + up)
vmax = vmax/((1 + exp(0.3*(10.0 - t)))*(0.8 + exp(0.14*(t - 36.0))))
vmax = vmax*1.0e-6*1.25

!----------------------------------------------------------------------!
! 'a' and 'b' are such that 'jc=a pc  - b'.                            !
!----------------------------------------------------------------------!
a0 = kp*kmq**((t - 25.0)/10.0)/pg
b0 = l*kmq**((t - 25.0)/10.0)/pg

!----------------------------------------------------------------------!
! Reduce equations to a cubic in 'a'.                                  !
!----------------------------------------------------------------------!
a1 =-b0*sc
b1 = a0*sc
c1 =-0.5*po
d1 = tau
e1 =-tau*rd*sc
f1 = tau

a2 = g0*kg*(pa**2.0)
b2 = g0*kg*pa
c2 =-pa*160.0 + g1*rh*kg*pa
d2 = g1*rh*kg

a3 = a1*b2 + a2*b1
b3 = a1*d2 + c2*b1
c3 = c1*b2 + a2*d1
d3 = c1*d2 + c2*d1

f2 =-f1

a4 = f2*c2*d2
b4 = f2*(a2*d2 + b2*c2) + d3*b3 + e1*c2*d2
c4 = f2*a2*b2 + c3*b3 + d3*a3 + e1*(a2*d2 + c2*b2)
d4 = c3*a3 + e1*a2*b2

!----------------------------------------------------------------------!
! Cubic coefficients 'p,q,r'.                                          !
!----------------------------------------------------------------------!
p = b4/a4
q = c4/a4
r = d4/a4

a2t = a2
!----------------------------------------------------------------------!
! Sove the cubic equaiton to find 'a'.                                 !
!----------------------------------------------------------------------!
Qt = (p**2 - 3.0*q)/9.0
Rt = (2.0*p**3 - 9.0*p*q + 27.0*r)/54.0

r2q3 = Rt**2 - Qt**3

if (r2q3<0.0) then
  arg1 = Rt/Qt**1.5
  if (arg1>1.0)  arg1 = 1.0
  if (arg1<-1.0)  arg1 = -1.0
  th = acos(arg1)
  a1 =-2.0*Qt**0.5*cos(th/3.0) - p/3.0
  a2 =-2.0*Qt**0.5*cos((th + 2.0*pi)/3.0) - p/3.0
  a3 =-2.0*Qt**0.5*cos((th + 4.0*pi)/3.0) - p/3.0
  a = a1
  if (a2>a)  a = a2
  if (a3>a)  a = a3
else
  at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5)**(1.0/3.0)
  if (abs(at)>1e-6) then
    bt = Qt/at
  else
    bt = 0.0
  endif
  a1 = at + bt - p/3.0
!        a2 =-0.5*(at + bt) - p/3.0 + i*3.0**0.5*(at - bt)/2.0
!        a3 =-0.5*(at + bt) - p/3.0 - i*3.0**0.5*(at - bt)/2.0
  a = a1
endif

gs = (g0 + g1*a*rh/pa)*kg
pc = pa - a*160.0/gs
w = a0*pc - b0
acheck = (w*(1.0 - 0.5*po/(tau*pc)) - rd)*sc

a2 = a2t

!----------------------------------------------------------------------!
! Compute light dependent assimilation rate.                           !
!----------------------------------------------------------------------!

wq = absx*qrub*f*qg*1.0
p = tau*c2
q = tau*a2 + wq*0.5*po*sc*d2 - tau*c2*(wq*sc - rd*sc)
r = wq*0.5*po*sc*b2 - tau*a2*(wq*sc - rd*sc)
aq = (-q + (abs(q)**2.0 - 4.0*p*r)**0.5)/(2.0*p)

!----------------------------------------------------------------------!
! Compute 'vmax' dependent assimilation rate.                          !
!----------------------------------------------------------------------!

wv = vmax*vmq**((t - 25.0)/10.0)/((1.0 + exp(0.3* &
 (10.0 - t)))*(0.8 + exp(0.14*(t - 36.0))))

q = tau*a2 + wv*0.5*po*sc*d2 - tau*c2*(wv*sc - rd*sc)
r = wv*0.5*po*sc*b2 - tau*a2*(wv*sc - rd*sc)
av = (-q + (abs(q)**2.0 - 4.0*p*r)**0.5)/(2.0*p)

!----------------------------------------------------------------------!
! Find limiting assimilation value.                                    !
!----------------------------------------------------------------------!
if (aq<a) then
  a = aq
  w = wq
  gs = (g0 + g1*a*rh/pa)*kg
  pc = pa - a*160.0/gs
endif

if (av<a) then
  a = av
  w = wv
  gs = (g0 + g1*a*rh/pa)*kg
  pc = pa - a*160.0/gs
endif

!**********************************************************************!
end subroutine assc42
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   xvjmax :: productivity_methods                     !
!                   ------------------------------                     !
!                                                                      !
! This subroutine calculates the value of 'amax' at the maximum        !
! iradience as opposed to the average irradiance in the main program.  !
!                                                                      !
! SUBROUTINE xvjmax(vmx,jmx,j,xt,xq,sum,nup,oi,dresp)                  !
!                                                                      !
!**********************************************************************!
subroutine xvjmax(vmx,jmx,j,xt,xq,sum,nup,oi,dresp)
!**********************************************************************!
real(dp) :: vmx,jmx,j,xt,xq,sum,nup,oi,tk,shapex,maxx,conv,cc,ea,num,&
 denom,up,kc,ko,tau,can,upt,am,vm,jm,irc,q,t,aa,bb,dresp
!----------------------------------------------------------------------!

t = 2.64 + 1.14*xt
tk = 273.0 + t
q = 2.0*xq

nup = 1.5*nup

shapex = 40.8 + 0.01*t - t**2*0.002
maxx = 0.738 - 0.002*t
if (nup<=0.0) then
  write(*,*) 'Problem in XVJMAX nup=',nup
  stop
endif
conv = 97.4115 - 2.5039*log(nup)

cc = 36.0*exp(-t*0.14) + 20.0
ea = 81000.0*exp(-t*0.14) + 50000.0
num = exp(shapex - conv/(0.00831*tk))
denom = 1.0 + exp((maxx*tk - 205.9)/(0.00831*tk))
up = num/denom

kc = exp(35.8 - 80.5/(0.00831*tk))
ko = exp(9.6 - 14.51/(0.00831*tk))*1000.0
tau = exp(-3.949 + 28.99/(0.00831*tk))
aa = 24.5/(24.5 + kc + kc*oi/ko)
bb = 0.5/(24.5 + kc + kc*oi/ko)


can = 1.0/sum
upt = up*can
am = upt*190.0/(360.0 + upt)

dresp = exp(cc-(ea/(8.3144*tk)))*upt/50.0

vm =-1.0*(am + 0.82)/(-1.0*aa + bb*oi/tau)
jm = 29.1 + 1.64*vm
vm = vm/1000000.0
jm = jm/1000000.0

if (t>6.34365) then
  jmx = jm*(1.0 + 0.0409*(t-25.0)- 1.54e-3*(t - 25.0)**2 - &
 9.42e-5*(t - 25.0)**3)
else
  jmx = jm*0.312633
endif

if (t>9.51718) then
  vmx = vm*(1.0 + 0.0505*(t - 25.0)- 0.248e-3*(t - 25.0)**2 - &
 8.09e-5*(t-25.0)**3)
else
  vmx = vm*0.458928
endif

! Mark we have to discuss about this calculation for irc here...
irc = q*exp(-0.5)
j = 0.24*irc/(1.0 + (0.24**2)*(irc**2)/(jmx**2))**0.5

!**********************************************************************!
end subroutine xvjmax
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   jcalc :: productivity_methods                      !
!                   -----------------------------                      !
!                                                                      !
! This subroutine calculates the electron transport rate               !
! extracted from main nppcalc subroutine APW 30 sept 2013              !
!                                                                      !
! SUBROUTINE jcalc(q,jmx,fw1984,j)                                     !
!                                                                      !
!**********************************************************************!
subroutine jcalc(q,jmx,fw1984,j)
!**********************************************************************!
real(dp) :: q,jmx,j,I2,qa,qb,qc
logical :: fw1984
!----------------------------------------------------------------------!

if (.not.fw1984) then
! use Harley 1992 J to Jmax realtionship
! j=0.24*q/(1.0+(0.24**2)*(q**2)/(jmx**2))**0.5
! in the function above 0.24 refers to the apparent quantum efficiency of light capture 
! and is the value used by Harley,  
! however Goudriaan's Law returns absorbed light so already accounts for reflection and transmission
! so the intrinsic quantum efficiency should be used, the 1.15 multiplier in the below function accounts for this
! j=1.15*0.24*q/(1.0+(0.24**2)*(q**2)/(jmx**2))**0.5
    j=0.24*q/(1.0 + (0.24**2)*(q**2)/(jmx**2))**0.5
else
! use Farquhar & Wong 1984 J to Jmax realtionship
  I2 = q*0.36    
  qa = 0.7
  qb = I2 + jmx
  qc = I2*jmx
  j  = (qb - (qb**2.0 - 4.0*qa*qc)**0.5)/(2.0*qa)
endif

!**********************************************************************!
end subroutine jcalc
!**********************************************************************!





!**********************************************************************!
! Brent solver from wikipedia
!**********************************************************************!
function brent_solver(func,i1,i2,oi,ca,vmx,j,rh,kg,rd,ga,t,p,gs_func, &
 g0,g1,dv,i,ci,light,ttype,ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,tmonth,jfv) 

real(dp) :: brent_solver,i1,i2,errortol
real(dp) :: a,b,c,d,s,fa,fb,fc,fs,tmp,tmp2,fa00
real(dp) :: fassv,fassj
real(dp) :: ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,tmonth
real(dp) :: oi,ca,rh,rd,vmx,j,p,ga,dv,t,cs,kg,g0,g1,ci,light
real(dp) :: kc,ko,tau,km,gstar,vt,jt,alp
integer :: func,gs_func,i,n,farq_pars_func,ttype
logical :: mflag,done,jfv

! Error tolerance umol m-2 s-1
errortol = 1.0e-3
done     = .false. 
fs       = -1.0
farq_pars_func = 1

! set parameter values
CALL FARQ_PARS(t,kc,ko,tau,alp,farq_pars_func)
km    = kc*(1 + oi/ko)
gstar = 0.5*oi/tau 
     
! initial guess  
a  = i1
b  = i2
! if (func.eq.0) then
! Error tolerance mol m-2 s-1
errortol = 1.0e-7
! vcmax25 & jmax25 to leaf t scalar 
vt = T_SCALAR(tmonth,'v',ttype,ftToptV,ftHaV,ftHdV,tmonth,jfv)
jt = T_SCALAR(tmonth,'j',ttype,ftToptJ,ftHaJ,ftHdJ,tmonth,jfv)

fa = VCMAX_MAIRE(a,vt,jt,ci,km,gstar,alp,light,farq_pars_func,i)
fb = VCMAX_MAIRE(b,vt,jt,ci,km,gstar,alp,light,farq_pars_func,i)
!      elseif (func.eq.1) then
!        fa = a -FASSV(a,gstar,km,ca,vmx,kg,rd,ga,t,p,gs_func,g0,g1,dv,i)
!        fb = b -FASSV(b,gstar,km,ca,vmx,kg,rd,ga,t,p,gs_func,g0,g1,dv,0)
!        if(fa.ge.0.d0) then
!          ! in this case rd is greater than gross a 
!          ! therefore assume a = 0 and anet = rd i
!          b = 0.d0
!          done = .true.
!        endif
!      elseif(func.eq.2) then
!        fa = a - FASSJ(a,gstar,ca,kg,rd,ga,t,p,j,gs_func,g0,g1,dv,i)
!        fb = b - FASSJ(b,gstar,ca,kg,rd,ga,t,p,j,gs_func,g0,g1,dv,0)
!        if(fa.ge.0.d0) then
!          ! in this case rd is greater than gross a 
!          ! therefore assume a = 0 and anet = rd i
!          b = 0.d0
!          done = .true.
!        endif
!      else 
!        print*, 'incorrect solver function number specified'
!      endif
      
fa00 = fa 
! if(i.eq.1) print*, i,fa,fb
! if f(a) f(b) >= 0 then error-exit
if ((.not.done).and.(fa*fb.ge.0).and.(func.gt.0)) then
  write(*,*) 'f(a) and f(b) do not span 0'
  write(*,*) 'a =',a,'fa =',fa,'b =',b,'fb =',fb
  write(*,*) 'solver function:', func 
  stop
elseif ((.not.done).and.(abs(fa).lt.abs(fb))) then
! if |f(a)| < |f(b)| then swap (a,b) end if
  tmp = a
  a   = b
  b   = tmp
  tmp = fa
  fa  = fb
  fb  = tmp
endif
 
c     = a
fc    = fa
mflag = .TRUE.
 
n = 0 
do 
  if(done.or.(fb.eq.0.0).or.(fs.eq.0.0).or.(abs(b-a).lt.errortol)) exit 

  if((fa.ne.fc).and.(fb.ne.fc)) then
! Inverse quadratic interpolation
    s = a*fb*fc/(fa - fb)/(fa - fc) + b*fa*fc/ &
 (fb - fa)/(fb - fc) + c*fa*fb/(fc - fa)/(fc - fb)
  else
! Secant Rule
    s = b - fb*(b - a)/(fb - fa)
  endif

  tmp2 = (3.0*a + b)/4.0
  if((.not.(((s.gt.tmp2).and.(s.lt.b)).or.((s.lt.tmp2).and.(s.gt.b)))).or. &
 (mflag.and.(abs(s - b).ge.(abs(b - c)/2.0))).or. &
 (.not.mflag.and.(abs(s - b).ge.(abs(c - d)/2.0))).or. &
 (mflag.and.(abs(b - c).lt.errorTol)).or. &
 (.not.mflag.and.(abs(c - d).lt.errorTol))) then
    s = (a + b)/2.0
    mflag = .TRUE.
  else
    mflag = .FALSE.
  endif 
    
! calculates a new value based on functions is the new value 
! fs = function(s)
    
!  if (func.eq.0) then
  fs = VCMAX_MAIRE(s,vt,jt,ci,km,gstar,alp,light,farq_pars_func,i)
!  elseif (func.eq.1) then
!    fs=s-FASSV(s,gstar,km,ca,vmx,kg,rd,ga,t,p,gs_func,g0,g1,dv,0)
!  elseif(func.eq.2) then
!    fs = s - FASSJ(s,gstar,ca,kg,rd,ga,t,p,j,gs_func,g0,g1,dv,0)
!  endif  

  d  = c
  c  = b
  fc = fb
  if ((fa * fs).lt.0) then
    b  = s
    fb = fs
  else
    a  = s
    fa = fs
  endif 

! if |f(a)| < |f(b)| then swap (a,b)
  if (abs(fa).lt.abs(fb)) then
    tmp = a
    a   = b
    b   = tmp
    tmp = fa
    fa  = fb
    fb  = tmp
  endif   

  n = n + 1
enddo
     
brent_solver = b

end function brent_solver






!----------------------------------------------------------------------*
!                                                                      *
!                          FUNCTION MAIRE VCMAX                        *
!                          ********************                        *
!     Gives vcmax AT LEAF TEMPERATURE  in mol m-2 s-1                  *  
!     Calculate vcmax that satisfies the condition that Wc = Wj        *
!     under the past months environmental conditions                   *
!     see Maire et al 2012 PLoS ONE                                    *
!                                                                      *
!----------------------------------------------------------------------*
function VCMAX_MAIRE(vm,vt,jt,ci,km,gstar,alpha,light,farq_pars_func,i)
  
real(dp) :: vcmax_maire,t,ci,oi,light,vm
real(dp)  :: vt,jt,km,gstar,alpha
integer :: farq_pars_func,i,ttype

vcmax_maire = vm*(1.0+(alpha*light/(jt*(exp(1.0)*(vm/vt)**0.89)))**2.0)**0.5 &
 *(4*ci+8*gstar) - alpha*light*(ci+km) 

! vcmax_maire = alpha*light / 
! &(1.d0+(alpha*light/(jt*(exp(1.d0)*(vm/vt)**0.89d0)))**2.d0)**0.5d0
! &*(ci+km)/(4*ci+8*gstar)  
! if(i.eq.1) print*, 'maire vcmax calc'
! if(i.eq.1) print'(4f8.4)', vt,jt,km,gstar
! if(i.eq.1) print'(3f14.8)', vcmax_maire,ci,light

end function vcmax_maire





!----------------------------------------------------------------------*
!                                                                      *
!                          SUBROUTINE FARQ_PARS                        *
!                          *****************!                          *
!     Gives M-M parameters from Farquhar model of photosynthesis       *  
!     Parameters Kc and Ko given in Pa                                 *
!                                                                      *
!----------------------------------------------------------------------*
subroutine FARQ_PARS(t,kc,ko,tau,alpha,farq_pars_func)

real(dp) :: t,tk,kc,ko,tau,alpha,c_p,oi
integer :: farq_pars_func

tk = t + 273.15  
oi = 21000.0 

! see explanation for 1.15 multiplier in subroutine JCALC       
! alpha = 1.150*0.240
alpha = 0.24

! below are the kinetic parameters from Farquhar etal 1980
kc  = exp(35.80 - 80.50/(0.008310*tk))
ko  = exp(9.60 - 14.510/(0.008310*tk))*1000.00
tau = exp(-3.9490 + 28.990/(0.008310*tk))

! the above are now deprecated for the in vivo parameters from Bernacchi etal 2001
!kc  = 40.490 * exp((79430.00 / (8.310 * 
!&298.150)) * (1.00 - (298.150) / (273.150 + t)))
!ko  = 27840.0 * exp((36380.00 / (8.310 *
!&298.150)) * (1.00 - (298.150) / (273.150 + t)))
!c_p = 4.2750 * exp((37830.00 / (8.310 * 
!&298.150)) * (1.00 - (298.150) / (273.150 + t)))
!tau = 0.50*oi/c_p

end subroutine farq_pars



end module productivity_methods
