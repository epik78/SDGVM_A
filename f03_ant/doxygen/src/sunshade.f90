module sunshade

use real_precision
implicit none

contains

subroutine SET_GOUD_PARAMS()
real(dp) :: soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo
common /GOUD_PARAMS/soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo

soilalbedo = 0.15
leafalbedo = 0.15
kbeam = 0.5
kdiff = 0.66
m = sqrt(1.0-leafalbedo)
kbeamstar = m*kbeam
canopyalbedo=2.0*kbeam/(kbeam+kdiff)*(1.0-m)/(1.0+m)

end subroutine SET_GOUD_PARAMS





subroutine GOUDRIAANSLAW2(lyr,lai,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade, &
 soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo,albedobeam,albedodiff)
real(dp) :: lyr,lai,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade
real(dp) :: soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo,albedobeam,albedodiff

! soilalbedo: value from Wang
! leafalbedo: leaf transmittance and reflectance are equal. Then albedo is the double
! kbeam: extinction coefficent of the direct beam
! kdiff: extinction coefficient of diffuse radiations
! I found somewhere an approximate relation: kbeam=0.75*kdiff
!... hum, I don't know if albedodiff is REAL*8ly different or not... so equal

qshade=(1.0-albedodiff)*diffrad*kdiff*exp(-kdiff*lyr) &
 +(1.0-albedobeam)*beamrad*kbeam*exp(-kbeam*lyr) &
 -(1.0-leafalbedo)*beamrad*kbeamstar*exp(-kbeamstar*lyr)
if (qshade<0.0) qshade=0.0

qsunlit=(1.0-leafalbedo)*kbeamstar*beamrad+qshade

fsunlit=exp(-kbeam*lyr)
fshade=1-fsunlit

end subroutine GOUDRIAANSLAW2





subroutine goudriaanslaw_ant(lyr,lai,beamrad,diffrad,fsunlit,qsunlit, &
 fshade,qshade,can_clump,cos_zen,s070607)
      
real(dp) :: lyr,lai,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade
real(dp) :: rhosoil,leafscattering,rhoh,rhobeam,rhodiff,rhobeam_can,rhodiff_can
real(dp) :: m,can_clump,cos_zen,Gbeam,canopyalbedo,kbeam,kdiff,kbeamprime,kdiffprime
logical :: s070607

if (.not.(s070607)) then

! Goudriaan's Law as described in Walker ... from Spitters 1986 and 
! Wang 2003 Func. Plant Biol. 
! all below values are for visible wavelengths 

! extinction coefficents of direct & diffuse radiation 
! these account for both canopy clumping and solar zenith angle, also 
! inform by Bodin & Franklin 2012 GMD  
! transmittance equals reflectance
  leafscattering = 2.0*0.075
  m          = sqrt(1.0 - leafscattering)
  Gbeam      = 0.5*can_clump
  kbeam      = Gbeam/cos_zen
  kbeamprime = m*kbeam
! kdiff is Gbeam * can_clump / cos_zen integrated over zenith angle 0 to pi/2  
  kdiff      = 0.8*can_clump 
  kdiffprime = m*kdiff

! calculate albedos
  rhosoil     = 0.15
  rhoh        = (1.0 - m)/(1.0 + m)
  rhobeam_can = rhoh*2.0*kbeam/(kbeam+kdiff)
  rhodiff_can = 4.0*Gbeam*rhoh*(Gbeam*(log(Gbeam)-log(Gbeam+kdiff))/kdiff**2 + 1.0/kdiff)
  rhobeam     = rhobeam_can + (rhosoil-rhobeam_can)*exp(-2.0*kbeamprime*lai)
  rhodiff     = rhodiff_can + (rhosoil-rhodiff_can)*exp(-2.0*kdiffprime*lai)

! calculate absorbed direct & diffuse radiation  
  qshade = (1.0-rhodiff)*diffrad*kdiffprime*exp(-kdiffprime*lyr) &
 + (1.0-rhobeam)*beamrad*kbeamprime*exp(-kbeamprime*lyr) &
 - (1.0-leafscattering)*beamrad*kbeam*exp(-kbeam*lyr)

  if (qshade.LT.0.0) qshade = 0
  qsunlit = (1.0-leafscattering)*kbeam*beamrad + qshade

! calculate fraction sunlit vs shaded leaves
  fsunlit = exp(-kbeam*lyr)
  fshade  = 1.0 - fsunlit

else 

! value from Wang
  rhosoil = 0.15
! leaf transmittance and reflectance are equal. Then albedo is the double
  leafscattering = 2.0*0.075

! extinction coefficent of the direct beam
  kbeam = 0.5*can_clump/cos_zen
! extinction coefficient of diffuse radiations
  kdiffprime = 0.66 
! I found somewhere an approximate relation: kbeam=0.75*kdiff

  m = sqrt(1.0 - leafscattering)
  kbeamprime = m*kbeam

  rhobeam_can = 2.0*kbeam/(kbeam + kdiffprime)*(1.0 - m)/(1.0 + m)
  rhobeam = rhobeam_can + (rhosoil - rhobeam_can)*exp(-kbeam*lai)
! ... hum, I don't know if albedodiff is really different or not... so equal
  rhodiff = rhobeam

  qshade = (1.0 - rhodiff)*diffrad*kdiffprime*exp(-kdiffprime*lyr) &
 +(1.0 - rhobeam)*beamrad*kbeam*exp(-kbeam*lyr) &
 -(1.0 - leafscattering)*beamrad*kbeamprime*exp(-kbeamprime*lyr)
  if (qshade.LT.0.0) qshade = 0.0

  qsunlit = (1.0 - leafscattering)*kbeamprime*beamrad + qshade
  fsunlit = exp(-kbeam*lyr)
  fshade = 1 - fsunlit
endif

end subroutine goudriaanslaw_ant





subroutine GOUDRIAANSLAW(lyr,lai,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade)

real(dp) :: lyr,lai,beamrad,diffrad
real(dp) :: fsunlit,qsunlit,fshade,qshade

real(dp) :: soilalbedo,leafalbedo,canopyalbedo
real(dp) :: albedobeam,albedodiff
real(dp) :: m
real(dp) :: kbeam,kdiff,kbeamstar

!value from Wang
soilalbedo = 0.15
!leaf transmittance and reflectance are equal. Then albedo is the double
leafalbedo = 2.0*0.075

!extinction coefficent of the direct beam
kbeam = 0.5
!extinction coefficient of diffuse radiations
kdiff = 0.66
!I found somewhere an approximate relation: kbeam=0.75*kdiff

m = sqrt(1.0 - leafalbedo)
kbeamstar = m*kbeam

canopyalbedo = 2.0*kbeam/(kbeam + kdiff)*(1.0 - m)/(1.0 + m)
albedobeam=canopyalbedo + (soilalbedo - canopyalbedo)*exp(-kbeam*lai)
!... hum, I don't know if albedodiff is REAL*8ly different or not... so equal
albedodiff=albedobeam

qshade = (1.0 - albedodiff)*diffrad*kdiff*exp(-kdiff*lyr) &
 + (1.0 - albedobeam)*beamrad*kbeam*exp(-kbeam*lyr) &
 - (1.0 - leafalbedo)*beamrad*kbeamstar*exp(-kbeamstar*lyr)
if (qshade<0.0) qshade = 0.0

qsunlit = (1.0 - leafalbedo)*kbeamstar*beamrad + qshade


fsunlit = exp(-kbeam*lyr)
fshade = 1.0 - fsunlit

end subroutine GOUDRIAANSLAW





subroutine BEERSLAW(lyr,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade)
implicit none
real(dp) :: lyr,beamrad,diffrad
real(dp) :: fsunlit,qsunlit,fshade,qshade
real(dp) :: kbeam,kdiff

!extinction coefficent of the direct beam
kbeam=0.5
!extinction coefficient of diffuse radiations
kdiff=0.66     

fshade=1.0
fsunlit=0

qsunlit=0.0
qshade=beamrad*kbeam*exp(-kbeam*lyr)+diffrad*kdiff*exp(-kdiff*lyr)

end subroutine BEERSLAW





subroutine SUNFLECT(lyr,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade)
real(dp) :: lyr,beamrad,diffrad
real(dp) :: fsunlit,qsunlit,fshade,qshade
real(dp) :: kbeam,kdiff

!extinction coefficent of the direct beam
kbeam=0.5
!extinction coefficient of diffuse radiations
kdiff=0.66     

fsunlit=0.5
fshade=0.5

qsunlit=beamrad*kbeam*exp(-kbeam*lyr)+diffrad*kdiff*exp(-kdiff*lyr)
qshade=qsunlit

end subroutine SUNFLECT


end module sunshade
