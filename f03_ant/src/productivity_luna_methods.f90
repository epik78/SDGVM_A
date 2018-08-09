module productivity_luna_methods

use dims
use real_precision
use system_state
use site_parameters
use pft_parameters

implicit none

contains

subroutine assimilation_calc(c3,t,rn,soil2g,wtwp,fshade,fsunlit,qshade,qsunlit, &
 jshade,jsunlit,canga,tk,p, &
 rlai,rh,oi,ca,vmx,kg,ko,kc,tau,rd,ga,t154,xdvy,dv,px,oitau1,oitau2, &
 gsmin,kcoiko,upt,a,ci,gs,gs_type,ii)
integer :: c3,gs_type,ii
real(dp) :: t,rn,soil2g,wtwp,fshade,fsunlit,qshade,qsunlit,canga,tk,p,rlai,rh
real(dp) :: oi,ca,vmx,kg,ko,kc,tau,rd,ga,t154,xdvy,dv,px,oitau1,oitau2,gsmin,kcoiko
real(dp) :: a,ci,gs

real(dp) :: tpaj,tppcj,tpgsj,tppcv,tpgsv,tpgsc4,tpav,tpac4,tppcc4
real(dp) :: rht,ashade,asunlit,pcshade,pcsunlit,jsunlit,jshade,gsshade,gssunlit
real(dp) :: upt,xvmax,cs

save tpav,tppcv,tpgsv

if ((t>0.0).and.(rn>0.0).and.(soil2g>wtwp).and.(fshade+fsunlit>0.1)) then

!----------------------------------------------------------------------!
! Only recompute assimilation every few days to speed things up, but   !
! make sure it is calculated on the first day of any growing season.   !
!----------------------------------------------------------------------!
  rht = rh

  if (c3==1) then
! Calculates assimilation rate,internal pressure and stomatal
! conductance if Wc goes in Eq.1
    call ASSVMAX3(tpav,tppcv,cs,tpgsv,oi,ca,vmx,rh,kg,ko,kc,&
 tau,rd,ga,t,p,t154,xdvy,dv,px,oitau1,oitau2,kcoiko,gsmin,gs_type,ii)

!     make the sum of assimilation, internal concentration and stomatal
!     conductance.
    a = 0.0
    ci = 0.0
    gs = 0.0

    if (fshade>0.01) then
! Calculates assimilation rate,internal pressure and stomatal
! conductance if J goes in Eq.1
      call ASSJ3(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,&
 jshade,t154,xdvy,dv,px,oitau1,oitau2,gsmin,gs_type,ii)

!print*,'tpaj ',tpaj,gsmin,kg
! Checks which assimilation process is smaller and keeps
! those values of assimilation,partial pressure and stomatal conductance 
     call LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,tpgsj,&
 ashade,pcshade,gsshade)
      a = a + fshade*ashade
      ci = ci + fshade*pcshade
      gs = gs + fshade*gsshade
    endif

    if (fsunlit>0.01) then
      call ASSJ4(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,&
 jsunlit,t154,xdvy,dv,px,oitau1,oitau2,gsmin,gs_type,ii)
!print*,'tpaj ',tpaj
      call LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,tpgsj,&
 asunlit,pcsunlit,gssunlit)
      a = a + fsunlit*asunlit
      ci = ci + fsunlit*pcsunlit
      gs = gs + fsunlit*gssunlit
    endif
!    print*,tpav,tpaj,ii,fsunlit

  else

    tpac4 = 0.0
    tppcc4 = 0.0
    tpgsc4 = 0.0

    if (fshade>0.01) then
      call ASSC4(ashade,pcshade,gsshade,ca,xvmax,rht,t,qshade,p,upt)
      tpac4 = tpac4 + fshade*ashade
      tppcc4 = tppcc4 + fshade*pcshade
      tpgsc4 = tpgsc4 + fshade*gsshade
    endif
    if (fsunlit>0.01) then
      call ASSC4(asunlit,pcsunlit,gssunlit,ca,xvmax,rht,t,qsunlit,p,upt)
      tpac4 = tpac4 + fsunlit*asunlit
      tppcc4 = tppcc4 + fsunlit*pcsunlit
      tpgsc4 = tpgsc4 + fsunlit*gssunlit
    endif

    a = tpac4
    gs = tpgsc4
    ci = tppcc4

    if (gs<0.0)  gs = 0.0
    gs = gs/1000.0

  endif
!----------------------------------------------------------------------!
! End of assimilation calculation 'IF' statement.                      !
!----------------------------------------------------------------------!
else
  a = 0.0
  gs = 0.005
  ci = ca
endif

end subroutine assimilation_calc





!**********************************************************************!
!                                                                      !
!                      assj3 :: productivity_methods                   !
!                      -----------------------------                   !
!                                                                      !
! ASSVMAX computes the assimilation rate by forming a cubic equation   !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assj3(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,px,   !
! oitau1,oitau2,gsmin)                                                                     !
!                                                                      !
!**********************************************************************!
subroutine assj3(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,dv,px,&
 oitau1,oitau2,gsmin,gs_type,ii)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step,p,ga,x,y,&
 t,cs,fa,kg,ax,bx,cx,gsmin,t154,xdvy,dv,px,oitau1,oitau2
integer :: gs_type,ii
!----------------------------------------------------------------------!

if (ssv(ssp%cohort)%assj(ssp%lai,1,ii) > 0.0) then
  a = ssv(ssp%cohort)%assj(ssp%lai,1,ii)
  cs = ca - a*px/ga
  gs = gs_calc(pft(ssp%cohort)%g0,pft(ssp%cohort)%g1,t,kg,a,cs,dv, &
 t154,gs_type)
  ci = cs - a*px/gs
  gam = 1.0 - oitau2/ci
  a = (gam*j*ci*0.25/(ci + oitau1) - rd)*1.0e6
else
  call ASSJ(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,gs_type)
endif
ssv(ssp%cohort)%assj(ssp%lai,1,ii) = a

!**********************************************************************!
end subroutine assj3
!**********************************************************************!






!----------------------------------------------------------------------*
!                                                                      *
!                          FUNCTION T_SCALAR                           *
!                          *****************                           *
!     Gives vcmax or jmax temperature scalar                           *  
!                                                                      *
!----------------------------------------------------------------------*
FUNCTION T_SCALAR(t,jv,ttype,ftTopt,ftkHa,ftkHd,tmonth,jfv)

REAL(dp)    :: t,qt,t_scalar,ftTopt,ftkHa,ftkHd,ftHa,ftHd,tmonth
REAL(dp)    :: Tsk,Trk,R,deltaS,dS
INTEGER, intent(in) :: ttype 
LOGICAL, intent(in) :: jfv 
CHARACTER :: jv

if(ttype.eq.0) then
! SDGVM original
  qt = 2.3
  if (t.gt.30.0) t = 30.0
  t_scalar = qt**(t/10.0d0)/qt**(2.5)
elseif (ttype.ge.1) then
! modified Arrhenius
        
  R = 8.31446        

! convert from kJ to J
! parameter values are read in as kJ... to save space
  ftHa = ftkHa*1.0e3
  ftHd = ftkHd*1.0e3

  if(ttype.eq.2) then
! employ Kattge&Knorr scaling based on mean temp of previous month 
          
    ftHd = 2.0e5

    if(jv.eq.'v') then
      ftHa   = 71513.0
      deltaS = 668.390 - 1.070*tmonth
    else
      ftHa   = 49884.0
      deltaS = 659.700 - 0.750*tmonth
    endif
  else 
! use fixed delta S based on input PFT parameters
    deltaS = ftHd/(ftTopt+273.15) + (R*log(ftHa/(ftHd-ftHa)))
  endif

  Tsk = t + 273.15
  Trk = 298.15
  
! see Medlyn etal 2002 or Kattge&Knorr 2007  
  t_scalar = exp(ftHa*(Tsk-Trk)/(R*Tsk*Trk))*( &
 (1 + exp((Trk*deltaS-ftHd)/(Trk*R)))/(1 + exp((Tsk*deltaS-ftHd)/(Tsk*R))))  

! adjust jmax based on kattge&knorr Jmax to Vcmax ratio relationship to temp
  if((jv.eq.'j').and.(ttype.ge.2).and.jfv)  t_scalar = t_scalar*(2.59-0.035*tmonth)/1.715

else
  t_scalar = 0.0
  print*, 'temperature scalar type',ttype,'undefined'
  stop
endif


END function t_scalar





!**********************************************************************!
!                                                                      !
!                   assvmax3 :: productivity_methods                   !
!                   --------------------------------                   !
! ! ASSVMAX computes the assimilation rate by forming a cubic equation !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assvmax3(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p,  !
! t154,xdvy,px,oitau1,oitau2,gsmin,kcoiko)                             !
!                                                                      !
!**********************************************************************!
subroutine assvmax3(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p,&
 t154,xdvy,dv,px,oitau1,oitau2,kcoiko,gsmin,gs_type,ii)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,vmx,rh,ko,kc,tau,rd,fa0,faf,gam,step, &
 p,ga,x,y,t,a0,cs,af,fa,kg,ax,bx,cx,gsmin,t154,xdvy,dv,px,oitau1,&
 oitau2,kcoiko
integer :: gs_type,ii
!----------------------------------------------------------------------!

if (ssv(ssp%cohort)%assv(ssp%lai,ii) > 0.0) then
  a = ssv(ssp%cohort)%assv(ssp%lai,ii)
  cs = ca - a*px/ga
  gs = gs_calc(pft(ssp%cohort)%g0,pft(ssp%cohort)%g1,t,kg,a,cs,dv, &
 t154,gs_type)
  ci = cs - a*px/gs
  gam = 1.0 - oitau2/ci
  a = (gam*vmx*ci/(ci + kcoiko) - rd)*1.0e6
else
  call ASSVMAX(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p,gsmin, &
 gs_type)
endif
ssv(ssp%cohort)%assv(ssp%lai,ii) = a

!**********************************************************************!
end subroutine assvmax3
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   gs_calc :: productivity_methods                    !
!                   -------------------------------                    !
! SUBROUTINE gs_calc(g0,g1,t,kg,a,cs,dv,gs_type)   !
!                                                                      !
!**********************************************************************!
real(dp) function gs_calc(g0,g1,t,kg,a,cs,dv,t154,gs_type)
!**********************************************************************!
real(dp) :: g0,g1,t,kg,a,cs,dv,t154
integer :: gs_type
real(dp) :: x
!----------------------------------------------------------------------!


if (gs_type==0) then
!----------------------------------------------------------------------!
! 070607                                                               !
!----------------------------------------------------------------------!
  x = 3.0 + 0.2*t
  if (x<0.25) x = 0.25
  gs_calc = g0 + ((x*a/(1.0 + dv/1.5)/(cs*10.0 - t154)))*kg

elseif (gs_type==1) then
!----------------------------------------------------------------------!
! GS_WOOD, the default gs model is a form of the Leuning 1995 model    !
! assuming gstar = 1.54t, this is inconsistent with the gstar in the   ! 
! assimilations calculations                                           !
!----------------------------------------------------------------------!
  x = 3.0 + 0.2*t
  if (x<0.25) x = 0.25
  gs_calc = (g0 + (x*a/(1.0 + dv/1.5)/(cs*10.0 - t154)))*kg

elseif (gs_type==2) then
!----------------------------------------------------------------------!
! From Ant's model, was named GS_MED, wasn't commented.                !
!----------------------------------------------------------------------!
  gs_calc = (g0 + (g1/(1.0 + dv**0.5))*(a/(cs/0.1013)))*kg

else
  write(*,*) 'Error: gs_type does not exist. gs_type = ',gs_type
  stop
endif



!**********************************************************************!
end function gs_calc
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assvmax :: productivity_methods                    !
!                   -------------------------------                    !
! ASSVMAX computes the assimilation rate by forming a cubic equation   !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assvmax(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p)   !
!                                                                      !
!**********************************************************************!
subroutine assvmax(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p, &
 gsmin,gs_type)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,vmx,rh,ko,kc,tau,rd,fa0,faf,gam,step, &
 p,ga,x,y,dv,t,a0,cs,af,fa,kg,ax,bx,cx,gsmin
integer :: i,gs_type
!----------------------------------------------------------------------!

gsmin = 0.005

dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)

step = 1.0
a = 0.0
cs = ca - a*1.0e-6*p/ga
gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)

ci = cs - a*1.0e-6*p/gs
gam = 1.0 - 0.5*oi/tau/ci
fa0 = a - (gam*vmx*ci/(ci + kc*(1.0 + oi/ko)) - rd)*1.0e6

a0 = a

50    continue
  a = a + step
  cs = ca - a*1.0e-6*p/ga
  gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  faf = a - (gam*vmx*ci/(ci + kc*(1.0 + oi/ko)) - rd)*1.0e6
if (faf<0.0) then
  fa0 = faf
  a0 = a
  goto 50
endif
af = a

if (af>100) then
  a = 0.0
  cs = ca - a*1.0e-6*p/ga
  gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
else
  a = (a0 + af)/2.0
  cs = ca - a*1.0e-6*p/ga
  gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  fa = a - (gam*vmx*ci/(ci + kc*(1.0 + oi/ko)) - rd)*1.0e6

  bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
  ax = (faf-fa-bx*(af-a))/(af**2-a**2)
  cx = fa-bx*a-ax*a**2

  if (abs(ax)>0.0) then
    if (bx**2-4.0*ax*cx>0.0) then
      a = (-bx+(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
  if (a>af)  a = (-bx-(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
    else
      a = 0.0
    endif
  else
    if (abs(bx)>0.0) then
      a =-cx/bx
    else
      a = 0.0
    endif
  endif

endif

cs = ca - a*1.0e-6*p/ga
gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
ci = cs - a*1.0e-6*p/gs
gam = 1.0 - 0.5*oi/tau/ci

!print*,'gs ',a,gs

!**********************************************************************!
end subroutine assvmax
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                     assj4 :: productivity_methods                    !
!                     -----------------------------                    !
!                                                                      !
! ASSVMAX computes the assimilation rate by forming a cubic equation   !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assj4(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,px,&  !
! oitau1,oitau2,gsmin)                                                 !
!                                                                      !
!**********************************************************************!
subroutine assj4(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,dv,px,&
 oitau1,oitau2,gsmin,gs_type,ii)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step,p,ga,x,y, &
 t,cs,fa,kg,ax,bx,cx,gsmin,t154,xdvy,dv,px,oitau1,oitau2
integer :: gs_type,ii
!----------------------------------------------------------------------!

if (ssv(ssp%cohort)%assj(ssp%lai,2,ii) > 0.0) then
  a = ssv(ssp%cohort)%assj(ssp%lai,2,ii)
  cs = ca - a*px/ga
  gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
  ci = cs - a*px/gs
  gam = 1.0 - oitau2/ci
  a = (gam*j*ci*0.25/(ci + oitau1) - rd)*1.0e6
else
  call ASSJ(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,gs_type)
endif
ssv(ssp%cohort)%assj(ssp%lai,2,ii) = a

!**********************************************************************!
end subroutine assj4
!**********************************************************************!




!**********************************************************************!
!                                                                      !
!                   assj :: productivity_methods                       !
!                   ----------------------------                       !
!                                                                      !
!! ASSVMAX computes the assimilation rate by forming a cubic equation  !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assj(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)                 !
!                                                                      !
!**********************************************************************!
subroutine assj(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,gs_type)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step,p,ga,x,y,&
 dv,t,cs,fa,kg,ax,bx,cx,gsmin
integer :: gs_type
!----------------------------------------------------------------------!

gsmin = 0.005

x = 3.0 + 0.2*t
if (x<0.25) x = 0.25
y = 1.5

dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)

step = 1.0
a = 0.0
cs = ca - a*1.0e-6*p/ga
gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)

ci = cs - a*1.0e-6*p/gs
gam = 1.0 - 0.5*oi/tau/ci
fa0 = a - (gam*j*ci/4.0/(ci + oi/tau) - rd)*1.0e6

a0 = a

50    continue
  a = a + step
  cs = ca - a*1.0e-6*p/ga
  gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  faf = a - (gam*j*ci/4.0/(ci + oi/tau) - rd)*1.0e6
if (faf<0.0) then
  fa0 = faf
  a0 = a
  goto 50
endif
af = a

if (af>100) then
  a = 0.0
  cs = ca - a*1.0e-6*p/ga
  gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
else

  a = (a0 + af)/2.0
  cs = ca - a*1.0e-6*p/ga
  gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  fa = a - (gam*j*ci/4.0/(ci + oi/tau) - rd)*1.0e6

  bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
  ax = (faf-fa-bx*(af-a))/(af**2-a**2)
  cx = fa-bx*a-ax*a**2

  if (abs(ax)>0.0) then
    if (bx**2-4.0*ax*cx>0.0) then
      a = (-bx+(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
  if (a>af)  a = (-bx-(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
    else
      a = 0.0
    endif
  else
    if (abs(bx)>0.0) then
      a =-cx/bx
    else
      a = 0.0
    endif
  endif

endif

cs = ca - a*1.0e-6*p/ga
gs = gs_calc(gsmin,0.0_dp,t,kg,a,cs,dv,t*1.54,gs_type)
ci = cs - a*1.0e-6*p/gs
gam = 1.0 - 0.5*oi/tau/ci


!print*,'gs ',a,gs
!stop

!**********************************************************************!
end subroutine assj
!**********************************************************************!






!**********************************************************************!
!                                                                      !
!              limitationprocess :: productivity_methods               !
!              -----------------------------------------               !
!                                                                      !
!     determin which process light, or Rubisco is the limiting         !
!     factor and return the assimilation, partial presure and          !
!     stomatal conductance corresponding to the limiting process       !
!                                                                      !
! SUBROUTINE limitationprocess(av,pcv,gsv,aj,pcj,gsj,a,pc,gs)          !
!                                                                      !
!**********************************************************************!
subroutine limitationprocess(av,pcv,gsv,aj,pcj,gsj,a,pc,gs)
!**********************************************************************!
real(dp) :: av,pcv,gsv,aj,pcj,gsj,a,pc,gs
!----------------------------------------------------------------------!

if (av<aj) then
   a = av
   gs = gsv
   pc = pcv
else
   a = aj
   gs = gsj
   pc = pcj
endif

!**********************************************************************!
end subroutine limitationprocess
!**********************************************************************!









!**********************************************************************!
!                                                                      !
!                   assc4 :: productivity_methods                      !
!                   -----------------------------                      !
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
! SUBROUTINE ASSC4(a,pc,gs,pa,vmax,rh,t,qg,p,up)                       !
!                                                                      !
!**********************************************************************!
subroutine ASSC4(a,pc,gs,pa,vmax,rh,t,qg,p,up)
!**********************************************************************!
real(dp) :: a,pc,gs,pa,vmax,rh,t,qg,up,p,vmq,absx,qrub,f,jqq,dv,amxt,&
 ac4t
!----------------------------------------------------------------------!

absx = 0.8
qrub = 0.11
f = 0.6
vmq = 2.0

dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)/101.325
! kilo pascals to mol/mol = divide by 101.325

pc = pa*(0.5 - (1.6*dv*p/2500.0/pa)**0.5)

vmax = up*360.0*vmq**((t - 25.0)/10.0)/(360.0 + up)
vmax = vmax/((1 + exp(0.3*(10.0 - t)))*(0.8 + exp(0.14*(t - 36.0))))
vmax = vmax*1.0e-6*1.25

jqq = absx*qrub*f*qg*1.0e+6
vmq = 2.0

amxt = up*190.0/(360.0 + up)
ac4t = amxt*vmq**((t - 25.0)/10.0)/(1.0 + &
 exp(0.3*(10.0 - t)))/(0.8 + exp(0.14*(t - 36.0)))
ac4t = up/2.0*vmq**((t - 25.0)/10.0)/(1.0 +  &
 exp(0.3*(15.0 - t)))/(0.8 + exp(0.14*(t - 33.0)))


if (ac4t<jqq) then
  a = ac4t
else
  a = jqq
endif

gs = a*1.6*p/(pa - pc)*1.0e-3

!**********************************************************************!
end subroutine ASSC4
!**********************************************************************!






!-------------------------------------------------------------------------------------------------------------------------------------------------       
function max_daily_pchg(tleaf10) ! - define this as a function to return max_daily_pchg  
! !LOCAL VARIABLES:
! local pointers to implicit in variables
real(dp)               :: max_daily_pchg                  ! maximum daily percentrage change  for nitrogen allocation
real(dp), intent (in)  :: tleaf10                         ! 10-day running mean of leaf temperature (oC)
real(dp), parameter    :: Q10Enz = 2.0                 ! Q10 value for enzyme decay rate
real(dp), parameter    :: Enzyme_turnover_daily = 0.1  ! the daily turnover rate for photosynthetic enzyme at 25oC in view of ~7 days of half-life time for Rubisco (Suzuki et al. 2001)
real(dp)               :: EnzTurnoverTFactor              ! temperature adjust factor for enzyme decay
!-------------------------------------------------------------------------------------------------------------------------------------------------       
!calculate the enzyme ternover rate
EnzTurnoverTFactor = Q10Enz**(0.1*(min(40.0, tleaf10)- 25.))
max_daily_pchg     = EnzTurnoverTFactor * Enzyme_turnover_daily

end function max_daily_pchg
!-------------------------------------------------------------------------------------------------------------------------------------------------       



end module productivity_luna_methods




