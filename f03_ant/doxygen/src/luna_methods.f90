module luna_methods
!********************************************************************************************************************************************************************** 
!
! !DESCRIPTION:
! Calculates Nitrogen fractionation within the leaf, based on optimum calculated fractions in rubisco, cholorophyll, 
! Respiration and Storage. Based on Xu et al. 2012 and Ali et al 2015.In Review 
!
! REVISION HISTORY:
! version 1.0, by Chonggang Xu, Ashehad Ali and Rosie Fisher. July 14  2015.
! version 0.1, by Chonggang Xu, Ashehad Ali and Rosie Fisher. October 30 2014. 
!********************************************************************************************************************************************************************** 
     





use real_precision
use productivity_luna_methods

implicit none

contains
     
!********************************************************************************************************************************************************************** 
! this subroutine updates the photosynthetic capacity as determined by Vcmax25 and Jmax25
   
!subroutine Update_Photosynthesis_Capacity(bounds, fn, filterp, 
!&dayl_factor, atm2lnd_inst, temperature_inst, canopystate_inst, 
!&photosyns_inst, 
!&surfalb_inst, solarabs_inst, waterstate_inst, frictionvel_inst)
     
!-------------------------------------------------------------------------------------------------------------------------------------------------       
!Use the LUNA model to calculate the Nitrogen partioning 
!subroutine NitrogenAllocation(FNCa,forc_pbot10, relh10, CO2a10,O2a10, PARi10,PARimx10,rb10, hourpd, tair10, tleafd10, tleafn10, &
!     Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp,&
!     PNstoreold, PNlcold, PNetold, PNrespold, PNcbold, &
!     PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
subroutine LUNA(z,vcmx25_z,jmx25_z,PNlc_z,enzs_z,rb10,sla,lnc,par240d,par240x, &
 CO2a10,o2a10,hourpd,relh10_in,tair10,gs_func,ftg0,ftg1,max_daily_pchg,ttype, &
 ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ) 
     
implicit none
     
!real(dp), intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
real(dp), intent (in) :: sla                        !Specific leaf area m2 g-1 DW
real(dp), intent (in) :: lnc                        !Area based leaf nitrogen content (g N/m2 leaf)
real(dp), intent (in) :: relh10_in                  !10-day mean relative humidity (unitless)
real(dp), intent (in) :: CO2a10                     !10-day meanCO2 concentration in the air (Pa)
real(dp), intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
real(dp), intent (in) :: par240d                    !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp), intent (in) :: par240x                    !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp), intent (in) :: rb10                       !10-day mean boundary layer resistance
real(dp), intent (in) :: hourpd                     !hours of light in a the day (hrs)
real(dp), intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
real(dp), intent (in) :: max_daily_pchg             ! maximum daily proportional change in vcmax or Jmax
! the below thre variables are commented out as the t pars are assumed the same as tair and pressure is assumed constant in SDGVM
!real(dp), intent (in) :: tleafd10                   !10-day running mean of daytime leaf temperature (oC) 
!real(dp), intent (in) :: tleafn10                   !10-day running mean of nighttime leaf temperature (oC) 
!real(dp), intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
! the below variables are commented out because they are declared as parameters below
!real(dp), intent (in) :: Jmaxb0                     !baseline proportion of nitrogen allocated for electron transport rate (unitless)
!real(dp), intent (in) :: Jmaxb1                     !coefficient determining the response of electron transport rate to light availability (unitless) 
!real(dp), intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
!real(dp), intent (in) :: relhExp                    !specifies the impact of relative humidity on electron transport rate (unitless)
! the below variables are commented as they declared as local variables below 
!real(dp), intent (in) :: PNstoreold                 !old value of the proportion of nitrogen allocated to storage (unitless)
!real(dp), intent (in) :: PNlcold                    !old value of the proportion of nitrogen allocated to light capture (unitless)
!real(dp), intent (in) :: PNetold                    !old value of the proportion of nitrogen allocated to electron transport (unitless)
!real(dp), intent (in) :: PNrespold                  !old value of the proportion of nitrogen allocated to respiration (unitless)
!real(dp), intent (in) :: PNcbold                    !old value of the proportion of nitrogen allocated to carboxylation (unitless)  
!real(dp), intent (out):: PNstoreopt                 !optimal proportion of nitrogen for storage 
!real(dp), intent (out):: PNlcopt                    !optimal proportion of nitrogen for light capture 
!real(dp), intent (out):: PNetopt                    !optimal proportion of nitrogen for electron transport 
!real(dp), intent (out):: PNrespopt                  !optimal proportion of nitrogen for respiration 
!real(dp), intent (out):: PNcbopt                    !optial proportion of nitrogen for carboxyaltion  
real(dp), intent (in)   :: ftToptV                  !optimum temp for vcmax (oC)
real(dp), intent (in)   :: ftHaV                    !activation energy for vcmax
real(dp), intent (in)   :: ftHdV                    !deactivation energy for vcmax
real(dp), intent (in)   :: ftToptJ                  !optimum temp for jmax (oC)
real(dp), intent (in)   :: ftHaJ                    !activation energy for jmax
real(dp), intent (in)   :: ftHdJ                    !deactivation energy for jmax
integer, intent (in)  :: z                        !canopy layer
integer, intent (in)  :: ttype                    !temperature scaling method to use
real(dp), intent (inout):: vcmx25_z                 !vcmax25
real(dp), intent (inout):: jmx25_z                  !jmx25
real(dp), intent (inout):: PNlc_z                   !optimal proportion of nitrogen for carboxylation  
real(dp), intent (inout):: enzs_z                   !enzyme decay state

!-------------------------------------------------------------------------------------------------------------------------------
!intermediate variables
real(dp) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
real(dp) :: SNCa                       !Area based structural nitrogen content (g N/m2 leaf)
real(dp) :: par240d_z                  !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp) :: par240x_z                  !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp) :: relh10                     !10-day mean relative humidity (unitless)
real(dp) :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp) :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp) :: PNlcold                    !old value of the proportion of nitrogen allocated to light capture (unitless)
real(dp) :: PNetold                    !old value of the proportion of nitrogen allocated to electron transport (unitless)
real(dp) :: PNrespold                  !old value of the proportion of nitrogen allocated to respiration (unitless)
real(dp) :: PNcbold                    !old value of the proportion of nitrogen allocated to carboxylation (unitless)  
real(dp) :: PNstoreopt                 !optimal proportion of nitrogen for storage 
real(dp) :: PNlcopt                    !optimal proportion of nitrogen for light capture 
real(dp) :: PNetopt                    !optimal proportion of nitrogen for electron transport 
real(dp) :: PNrespopt                  !optimal proportion of nitrogen for respiration 
real(dp) :: PNcbopt                    !optial proportion of nitrogen for carboxyaltion  
real(dp) :: Carboncost1                             !absolute amount of carbon cost associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
real(dp) :: Carboncost2                             !absolute amount of carbon cost associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
real(dp) :: Carbongain1                             !absolute amount of carbon gain associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
real(dp) :: Carbongain2                             !absolute amount of carbon gain associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
real(dp) :: Fc                                      !the temperature adjustment factor for Vcmax 
real(dp) :: Fj                                      !the temperature adjustment factor for Jmax 
real(dp) :: PNlc                                    !the current nitrogen allocation proportion for light capture
real(dp) :: Jmax                                    !the maximum electron transport rate (umol/m2/s) 
real(dp) :: JmaxCoef                                !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
real(dp) :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
real(dp) :: JmaxL                                   !the electron transport rate with maximum daily radiation (umol/m2/s)  
real(dp) :: JmeanL                                  !the electron transport rate with mean radiation (umol/m2/s) 
real(dp) :: Nstore                                  !absolute amount of nitrogen allocated to storage (gN/m2 leaf)
real(dp) :: Nresp                                   !absolute amount of nitrogen allocated to respiration (gN/m2 leaf) 
real(dp) :: Nlc                                     !absolute amount of nitrogen allocated to light capture (gN/m2 leaf) 
real(dp) :: Net                                     !absolute amount of nitrogen allocated to electron transport (gN/m2 leaf) 
real(dp) :: Ncb                                     !absolute amount of nitrogen allocated to carboxylation (gN/m2 leaf) 
real(dp) :: Nresp1                                  !absolute amount of nitrogen allocated to respiration due to increase in light capture nitrogen(gN/m2 leaf)  
real(dp) :: Nlc1                                    !absolute amount of nitrogen allocated to light capture due to increase in light capture nitrogen(gN/m2 leaf) 
real(dp) :: Net1                                    !absolute amount of nitrogen allocated to electron transport due to increase in light capture nitrogen(gN/m2 leaf)
real(dp) :: Ncb1                                    !absolute amount of nitrogen allocated to carboyxlation due to increase in light capture nitrogen(gN/m2 leaf) 
real(dp) :: Nresp2                                  !absolute amount of nitrogen allocated to respiration due to decrease in light capture nitrogen(gN/m2 leaf) 
real(dp) :: Nlc2                                    !absolute amount of nitrogen allocated to light capture due to decrease in light capture nitrogen(gN/m2 leaf) 
real(dp) :: Net2                                    !absolute amount of nitrogen allocated to electron transport due to decrease in light capture nitrogen(gN/m2 leaf) 
real(dp) :: Ncb2                                    !absolute amount of nitrogen allocated to carboxylation due to increase in light capture nitrogen(gN/m2 leaf) 
real(dp) :: PSN                                     !g carbon photosynthesized per day per unit(m2) of leaf
real(dp) :: RESP                                    !g carbon respired per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
real(dp) :: PSN1                                    !g carbon photosynthesized per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
real(dp) :: RESP1                                   !g carbon respired per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
real(dp) :: PSN2                                    !g carbon photosynthesized per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
real(dp) :: RESP2                                   !g carbon respired per day per unit(m2) of leaf
real(dp) :: Npsntarget                              !absolute amount of target nitrogen for photosynthesis(gN/m2 leaf) 
real(dp) :: Npsntarget1                             !absolute amount of target nitrogen for photosynthesis due to increase in light capture nitrogen(gN/m2 leaf) 
real(dp) :: Npsntarget2                             !absolute amount of target nitrogen for photosynthesis due to decrease in light capture nitrogen(gN/m2 leaf) 
real(dp) :: NUEj                                    !nitrogen use efficiency for electron transport under current environmental conditions 
real(dp) :: NUEc                                    !nitrogen use efficiency for carboxylation under current environmental conditions  
real(dp) :: NUEjref                                 !nitrogen use efficiency for electron transport under reference environmental conditions (25oC and 385ppm Co2) 
real(dp) :: NUEcref                                 !nitrogen use efficiency for carboxylation under reference environmental conditions (25oC and 385ppm Co2) 
real(dp) :: NUEr                                    !nitrogen use efficiency for respiration 
real(dp) :: PARi10c                                 !10-day mean constrained photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp) :: PARimx10c                               !10-day mean constrained 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp) :: Kj2Kcref                                !the ratio of rubisco-limited photosynthetic rate (Wc) to light limited photosynthetic rate (Wj)
real(dp) :: PNlcoldi                                !old value of the proportion of nitrogen allocated to light capture (unitless) 
real(dp) :: Kj2Kc                                   !the ratio of Wc to Wj under changed conditions 
real(dp) :: Kc                                      !conversion factors for Vc,max to Wc 
real(dp) :: Kj                                      !conversion factor for electron transport rate to Wj 
real(dp) :: theta                                   !efficiency of light energy conversion (unitless) 
real(dp) :: chg                                     !the nitrogen change per interation
real(dp) :: chg_constrn                             !the nitrogen change per interation
real(dp) :: chg_per_step                            !the nitrogen change per interation
real(dp) :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
real(dp) :: ci                                      !inter-cellular CO2 concentration (Pa)
!real(dp) :: theta_cj                                !interpolation coefficient
real(dp) :: tleafd10                   !10-day running mean of daytime leaf temperature (oC) 
real(dp) :: tleafn10                   !10-day running mean of nighttime leaf temperature (oC) 
real(dp) :: forc_pbot10                !10-day mean air pressure (Pa)
real(dp) :: tleafd10c                               !10-day mean daytime leaf temperature, contrained for physiological range (oC)
real(dp) :: tleafn10c                               !10-day mean leaf temperature for night, constrained for physiological range (oC)
real(dp) :: Vcmax                                   !the maximum carboxyaltion rate (umol/m2/s) 
real(dp) :: vcmx25_opt   
real(dp) :: jmx25_opt   
real(dp) :: rabsorb
real(dp) :: radmax2mean 
integer :: KcKjFlag                                !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
integer :: jj                                      !index record fo the number of iterations
integer :: increase_flag                           !whether to increase or decrease

!-------------------------------------------------------------------------------------------------------------------------------
!pass thru variables
real(dp), intent (in) :: ftg0                       !the intercept of the stomatal conductance equation
real(dp), intent (in) :: ftg1                       !the slope of the stomatal conductance equation 
integer, intent (in):: gs_func                    !flag to indicate which stomatal conductance function to use 

!-------------------------------------------------------------------------------------------------------------------------------
! functions
!real(dp) :: RespTBernacchi
!real(dp) :: t_scalar                                !temperature scaling factor for Vcmax or Jmax 

!------------------------------------------------------------------------------ 
!Constants  
real(dp), parameter :: rhol = 0.075                     ! leaf reflectance: 1=vis, 2=nir  
real(dp), parameter :: taul = 0.075                     ! leaf transmittance: 1=vis, 2=nir 
real(dp), parameter :: Cb   = 1.78                    ! nitrogen use effiency for choloraphyll for light capture, see Evans 1989  
real(dp), parameter :: Cv   = 1.2d-5 * 3600.          ! conversion factor from umol CO2 to g carbon
real(dp), parameter :: Kc25 = 40.49                   ! Mechalis constant of CO2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
real(dp), parameter :: Ko25 = 27840                   ! Mechalis constant of O2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
real(dp), parameter :: Cp25 = 4.275                   ! CO2 compensation point at 25C (Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
real(dp), parameter :: Fc25 = 294.2                   ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
real(dp), parameter :: Fj25 = 1257.0                  ! Fj25 = 8.06*156 # #see COSTE 2005 and Xu et al 2012
real(dp), parameter :: NUEr25 = 33.69                 ! nitrogen use efficiency for respiration, see Xu et al 2012
real(dp), parameter :: O2ref = 209460.0               ! ppm of O2 in the air
real(dp), parameter :: CO2ref = 380.0                 ! reference CO2 concentration for calculation of reference NUE. 
real(dp), parameter :: forc_pbot_ref = 101325.0       ! reference air pressure for calculation of reference NUE
real(dp), parameter :: Jmaxb0 = 0.0311                ! the baseline proportion of nitrogen allocated for electron transport (J)     
real(dp), parameter :: Jmaxb1 = 0.1745                ! the baseline proportion of nitrogen allocated for electron transport (J)    
real(dp), parameter :: Wc2Wjb0 = 0.8054               ! the baseline ratio of rubisco limited rate vs light limited photosynthetic rate (Wc:Wj) 
real(dp), parameter :: relhExp = 6.0999               ! electron transport parameters related to relative humidity
real(dp), parameter :: NMCp25 = 0.715                 ! estimated by assuming 80% maintenance respiration is used for photosynthesis enzyme maintenance
real(dp), parameter :: Trange1 = 5.0                  ! lower temperature limit (oC) for nitrogen optimization  
real(dp), parameter :: Trange2 = 42.0                 ! upper temperature limit (oC) for nitrogen optimization
real(dp), parameter :: SNC = 0.004                    ! structural nitrogen concentration (g N g-1 C)
real(dp), parameter :: mp = 9.0                       ! slope of stomatal conductance; this is used to estimate model parameter, but may need to be updated from the physiology file, 
!SDGVM PAR units are in molm-2s-1
real(dp), parameter :: PARLowLim = 200.0d-6             ! minimum photosynthetically active radiation for nitrogen optimization
real(dp), parameter :: minrelh = 0.25                 ! minimum relative humdity for nitrogen optimization
     
!-------------------------------------------------------------------------------------------------------------------------------------------------       
   
!print*, z,vcmx25_z,jmx25_z,PNlc_z,enzs_z,lnc,par240d
   
! SDGVM asssumes daily average temperatures for the leaf & constant pressure 
tleafd10    = tair10
tleafn10    = tair10
forc_pbot10 = 101325.
 
!SDGVM PAR units are in molm-2s-1
!par240d_z   = par240d * 1d6
!par240x_z   = par240x * 1d6
par240d_z   = par240d
par240x_z   = par240x

! LUNA expects RH in proportion i.e. 0-1, SDGVM in %      
relh10      = relh10_in * 1.d-2
 
! SLA is fixed through the canopy in SDGVM and lnc is pre-determined according to Beer's Law scaling or similar 
! - therefore relCLNCa, PARTop are not needed
! PAR is already passed to this routine in umol m-2 s-1 (check not day)
! - therefore the intermediary steps to convert from wm-2 are not necessary

!------------------------------------------------------------------
!SNCa     =  1.0/slatop(ft) * SNC !(Eq A3)
SNCa = 1.0/sla*0.48 * SNC     !(Eq A3)
FNCa = lnc - SNCa                 !(Eq A4)

!------------------------------------------------------------------
! for SDGVM the distinction between sunlit and non-sunlit leaves is not made in the LUNA model
! - this is because during the day leaves experience both sunlit and non-sunlit states so  
! par240d/x_z is 10-day mean/maximum PAR for leaves in a canopy layer
! rabsorb is ratio of absorbed light to incident light (rhol and taul are reflectance and transmittance)
rabsorb     = 1.0-rhol-taul
radmax2mean = par240x_z / par240d_z
PARi10      = par240d_z / rabsorb
PARimx10    = PARi10 * radmax2mean

!------------------------------------------------------------------
!nitrogen allocation model-start          
!PNlcold     = PNlc_z
PNlcold     = 0.2 / lnc 
PNetold     = 0.0
PNrespold   = 0.0
PNcbold     = 0.0                                     
!end brought in from above subroutine 

! this call to N allocation has been replaced by the N allocation code 
!call NitrogenAllocation(FNCa,forc_pbot10(p), relh10, CO2a10, O2a10, PARi10, PARimx10, rb10, hourpd, &
!     tair10, tleafd10, tleafn10, &
!     Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp,  PNstoreold, PNlcold, PNetold, PNrespold, &
!     PNcbold, PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)

! set referfence NUE parameters 
! - using sinlge common NUE function fed with default values as arguments
! call NUEref(NUEjref, NUEcref, Kj2Kcref)
!print*, 'LUNA start'
call NUE(O2ref*0.1013, 0.7*CO2ref*0.1013, 25.0_dp, 25.0_dp, &
 NUEjref, NUEcref, Kj2Kcref,ttype,ftToptV,ftHaV,ftHdV, &
 ftToptJ,ftHaJ,ftHdJ)

! appears theta_cj is not necessary for this subroutine
!theta_cj = 0.95
! it's not currently clear to me where PNstoreold is passed from in the original code
!Nstore   = PNstoreold * FNCa                         !FNCa * proportion of storage nitrogen in functional nitrogen
Nlc      = PNlcold    * FNCa                         !FNCa * proportion of light capturing nitrogen in functional nitrogen
Net      = PNetold    * FNCa                         !FNCa * proportion of light harvesting (electron transport) nitrogen in functional nitrogen
Nresp    = PNrespold  * FNCa                         !FNCa * proportion of respiration nitrogen in functional nitrogen
Ncb      = PNcbold    * FNCa                         !FNCa * proportion of carboxylation nitrogen in functional nitrogen
if (Nlc > FNCa * 0.5) Nlc = 0.5 * FNCa
PNlc     = PNlcold
PNlcoldi = PNlcold  - 0.001

! constrain the physiological range of PAR and t
PARi10c   = max(PARLowLim, PARi10)
PARimx10c = max(PARLowLim, PARimx10)
tleafd10c = min(max(tleafd10, Trange1), Trange2)  
tleafn10c = min(max(tleafn10, Trange1), Trange2) 

! initialse solver parameters 
ci            = 0.7 * CO2a10 
JmaxCoef      = Jmaxb1 * ((hourpd / 12.0)**2.0) * &
 ( 1.0 - exp( -relhExp * max(relh10 - minrelh, 0.0) / &
 (1.0 - minrelh) ) )
chg_per_step  = 0.02* FNCa
increase_flag = 0

!------------------------------------------------------------------
jj = 1
do while ( (PNlcoldi .NE. PNlc) .and. (jj .lt. 100) )      
     
 ! Fc is the scaling factor to go from leaf N invested in RuBisCO to Vcmax  
 ! Fj is the scaling factor to go from leaf N invested in elec trans to Jmax  
 !Fc   = VcmxTKattge(tair10, tleafd10c) * Fc25
 !Fj   = JmxTKattge(tair10, tleafd10c)  * Fj25
 Fc = T_SCALAR(tleafd10c,'v',ttype,ftToptV,ftHaV,ftHdV,tair10,.FALSE.)*Fc25
 Fj = T_SCALAR(tleafd10c,'j',ttype,ftToptJ,ftHaJ,ftHdJ,tair10,.FALSE.)*Fj25
 
 NUEr = Cv * NUEr25 * ( RespTBernacchi(tleafd10c) * hourpd + &
  RespTBernacchi(tleafn10c) * (24.0 - hourpd) ) !nitrogen use efficiency for respiration (g biomass/m2/day/g N)

!print*, 'LUNA do while loop'
 
  call NUE(O2a10, ci, tair10, tleafd10c, NUEj, NUEc, Kj2Kc,ttype, &
 ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ)
 
 KcKjFlag = 0
  call LUNA_Ninvestments (KcKjFlag,FNCa, Nlc, forc_pbot10,relh10, &
 CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, tair10, &
 tleafd10c,tleafn10c, &
 Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr, &
 Kc, Kj, ci, Vcmax, Jmax,JmeanL,JmaxL, Net, Ncb, Nresp, PSN, RESP, &
 gs_func,ftg0,ftg1)      

! target nitrogen allocated to photosynthesis, which may be lower or higher than Npsn_avail
  Npsntarget = Nlc + Ncb + Net       
  PNlcoldi   = Nlc / FNCa
  Nstore     = FNCa - Npsntarget - Nresp
 
!test the increase of light capture nitrogen
  if ( ((Nstore > 0.0) .and.(increase_flag .eq. 1)) .or. (jj .eq. 1) ) then
    Nlc2 = Nlc + chg_per_step
    if (Nlc2 / FNCa > 0.95) Nlc2 = 0.95 * FNCa
  
    KcKjFlag = 1
    call LUNA_Ninvestments (KcKjFlag,FNCa, Nlc2, forc_pbot10, &
 relh10, CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, & 
 tair10, tleafd10c,tleafn10c, &
 Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr, &
 Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL, Net2, Ncb2, Nresp2, PSN2, RESP2, &
 gs_func,ftg0,ftg1)
  
    Npsntarget2 = Nlc2 + Ncb2 + Net2
  
!update the nitrogen change
    Carboncost2 = (Npsntarget2 - Npsntarget) * NMCp25 * Cv * &
 ( RespTBernacchi(tleafd10c) * hourpd + RespTBernacchi(tleafn10c) * &
 (24.0  - hourpd) )
    Carbongain2 =  PSN2 - PSN
  
    if( (Carbongain2 > Carboncost2).and.(Npsntarget2 + Nresp2 < 0.95 * FNCa) ) then
      Nlc = Nlc2
      Net = Net2
      Ncb = Ncb2
      Nstore = FNCa - Npsntarget2 - Nresp2 
      if (jj == 1) increase_flag = 1
    endif
  endif
 
 !------------------------------------------------------------------------------------
 !test the decrease of light capture nitrogen
  if (increase_flag == 0) then  
 
    if (Nstore < 0.0) then
      Nlc1 = Nlc * 0.8 !bigger step of decrease if it is negative            
    else
      Nlc1 = Nlc - chg_per_step
    endif
 
    if (Nlc1 < 0.05) Nlc1 = 0.05
    KcKjFlag = 1
    call LUNA_Ninvestments(KcKjFlag,FNCa, Nlc1,forc_pbot10,relh10, &
 CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, &
 tair10, tleafd10c,tleafn10c, &
 Kj2Kc, Wc2Wjb0,JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr, &
 Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL,Net1, Ncb1, Nresp1, PSN1, RESP1, &
 gs_func,ftg0,ftg1)
  
    Npsntarget1 = Nlc1 + Ncb1 + Net1
    Carboncost1 = (Npsntarget - Npsntarget1) * NMCp25 * Cv * &
 ( RespTBernacchi(tleafd10c) * hourpd + RespTBernacchi(tleafn10c) * &
 (24.0  - hourpd) )
  Carbongain1 =  PSN - PSN1
  
    if( ((Carbongain1 < Carboncost1) .and. (Nlc1 > 0.05)) .or.&
 (Npsntarget + Nresp > 0.95 * FNCa) ) then
      Nlc    = Nlc1 
      Net    = Net1   
      Ncb    = Ncb1
      Nstore = FNCa - Npsntarget1 - Nresp1  
    endif
  
  endif
  PNlc = Nlc / FNCa
 
  jj = jj + 1  
enddo
!------------------------------------------------------------------
! record optimum proportions      
PNlcopt    = Nlc    / FNCa
PNstoreopt = Nstore / FNCa
PNcbopt    = Ncb    / FNCa
PNetopt    = Net    / FNCa
PNrespopt  = Nresp  / FNCa 

!brought in from parent subroutine
PNlc_z      = PNlcopt

! determine change in vcmax and jmax 
vcmx25_opt  = PNcbopt * FNCa * Fc25
jmx25_opt   = PNetopt * FNCa * Fj25

chg         = vcmx25_opt-vcmx25_z
chg_constrn = min(abs(chg),vcmx25_z*max_daily_pchg)
vcmx25_z    = vcmx25_z+sign(1.0_dp,chg)*chg_constrn
!if(z.eq.1) print*, vcmx25_opt,chg,vcmx25_z      
!print*, vcmx25_opt,chg,vcmx25_z      
 
chg         = jmx25_opt-jmx25_z
chg_constrn = min(abs(chg),jmx25_z*max_daily_pchg)
jmx25_z     = jmx25_z+sign(1.0_dp,chg)*chg_constrn 
!print*, jmx25_opt,chg,jmx25_z      

if(enzs_z<1.0) enzs_z = enzs_z * (1.0 + max_daily_pchg)

!nitrogen allocation subroutine end  
!------------------------------------------------------------------

if( isnan(vcmx25_z) .or. (vcmx25_z>1000.) .or. (vcmx25_z<0.) ) then
  !write(iulog, *) 'Error: Vc,mx25 become unrealistic (NaN,>1000,
  write(*, *) 'Error: Vc,mx25 become unrealistic (NaN,>1000, or negative) for z=', z
  write(*, *) 'Error: Vcx25:', vcmx25_z
  write(*, *) 'Error: Jmx25:', jmx25_z
  !write(iulog, *) 'LUNA env:',FNCa,forc_pbot10, relh10, CO2a10, 
  write(*, *) 'LUNA env:',FNCa,forc_pbot10, relh10, CO2a10, &
 O2a10, PARi10, PARimx10, rb10, hourpd, tair10, tleafd10, tleafn10
  !call endrun(msg=errmsg(__FILE__, __LINE__))
  stop
endif
if( isnan(jmx25_z) .or. (jmx25_z>1000.) .or. (jmx25_z<0.) ) then
 !write(iulog, *) 'Error: Jmx25 become unrealistic (NaN,>1000, 
  write(*, *) 'Error: Jmx25 become unrealistic (NaN,>1000,or negative)for z=', z
  write(*, *) 'Error: Jmx25:', jmx25_z
  write(*, *) 'Error: Vcx25:', vcmx25_z
  !write(iulog, *) 'LUNA env:', FNCa,forc_pbot10, relh10, CO2a10,
  write(*, *) 'LUNA env:', FNCa,forc_pbot10, relh10, CO2a10, &
 O2a10, PARi10, PARimx10, rb10,hourpd, tair10, tleafd10, tleafn10
 !call endrun(msg=errmsg(__FILE__, __LINE__))
  stop
endif
!end brought in from above subroutine 

end subroutine LUNA
!-------------------------------------------------------------------------------------------------------------------------------------------------       



!-------------------------------------------------------------------------------------------------------------------------------------------------       
subroutine LUNA_Ninvestments (KcKjFlag, FNCa, Nlc, forc_pbot10, &
 relh10,CO2a10,O2a10,PARi10,PARimx10,rb10,hourpd,tair10,tleafd10, &
 tleafn10,Kj2Kc, Wc2Wjb0, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, &
 NUEjref,NUEr,Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL,Net,Ncb,Nresp,PSN, &
 RESP,gs_func,ftg0,ftg1)
!calculate the nitrogen investment for electron transport, carb10oxylation, respiration given a specified value 
!of nitrogen allocation in light capture [Nlc]. This equation are based on Ali et al 2015b.

implicit none

integer,intent (in) :: KcKjFlag                   !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
real(dp), intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
real(dp), intent (in) :: Nlc                        !nitrogen content for light capture(g N/m2 leaf)
real(dp), intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
real(dp), intent (in) :: relh10                     !10-day mean relative humidity (unitless)
real(dp), intent (in) :: CO2a10                     !10-day mean CO2 concentration in the air (Pa)
real(dp), intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
real(dp), intent (in) :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp), intent (in) :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
real(dp), intent (in) :: rb10                       !10-day mean boundary layer resistance (s/m)
real(dp), intent (in) :: hourpd                     !hours of light in a the day (hrs)
real(dp), intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
real(dp), intent (in) :: tleafd10                   !10-day mean daytime leaf temperature (oC) 
real(dp), intent (in) :: tleafn10                   !10-day mean nighttime leaf temperature (oC) 
real(dp), intent (in) :: Kj2Kc                      !ratio:  Kj / Kc
real(dp), intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
real(dp), intent (in) :: JmaxCoef                   !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
real(dp), intent (in) :: Fc                         !the temperature adjusted ratio of Vcmax to N invested in Vcmax  
real(dp), intent (in) :: Fj                         !the temperature adjusted ratio of  Jmax to N invested in electron trans 
real(dp), intent (in) :: NUEc                       !nitrogen use efficiency for carboxylation 
real(dp), intent (in) :: NUEj                       !nitrogen use efficiency for electron transport
real(dp), intent (in) :: NUEcref                    !nitrogen use efficiency for carboxylation under reference climates
real(dp), intent (in) :: NUEjref                    !nitrogen use efficiency for electron transport under reference climates
real(dp), intent (in) :: NUEr                       !nitrogen use efficiency for respiration
real(dp), intent (inout) :: Kc                      !conversion factors from Vc,max to Wc 
real(dp), intent (inout) :: Kj                      !conversion factor from electron transport rate to Wj 
real(dp), intent (inout) :: ci                      !inter-cellular CO2 concentration (Pa) 
real(dp), intent (out)   :: Vcmax                   !the maximum carboxyaltion rate (umol/m2/s) 
real(dp), intent (out)   :: Jmax                    !the maximum electron transport rate (umol/m2/s) 
real(dp), intent (out)   :: JmaxL                   !the electron transport rate with maximum daily radiation (umol/m2/s)  
real(dp), intent (out)   :: JmeanL                  !the electron transport rate with mean radiation (umol/m2/s) 
real(dp), intent (out)   :: Net                     !nitrogen content for electron transport(g N/m2 leaf)
real(dp), intent (out)   :: Ncb                     !nitrogen content for carboxylation(g N/m2 leaf)
real(dp), intent (out)   :: Nresp                   !nitrogen content for respiration(g N/m2 leaf)
real(dp), intent (out)   :: PSN                     !daily photosynthetic rate(g C/day/m2 leaf)
real(dp), intent (out)   :: RESP                    !daily respiration rate(g C/day/m2 leaf)

!-------------------------------------------------------------------------------------------------------------------------------
!intermediate variables
real(dp) :: A                                       !Gross photosynthetic rate (umol CO2/m2/s)
real(dp) :: Wc2Wj                                   !ratio: Wc/Wj  
real(dp) :: ELTRNabsorb                             !absorbed electron rate, umol electron/m2 leaf /s
real(dp) :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
real(dp) :: theta_cj                                !interpolation coefficient
real(dp) :: theta                                   !light absorption rate (0-1)
real(dp) :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
real(dp) :: Wc                                      !rubisco-limited photosynthetic rate (umol/m2/s)
real(dp) :: Wj                                      !light limited photosynthetic rate (umol/m2/s)
real(dp) :: k_c      
real(dp) :: k_o       
real(dp) :: tau        
real(dp) :: c_p        
real(dp) :: awc        
real(dp) :: gs        
real(dp) :: NUECHG                                  !the nitrogen use efficiency change under current conidtions compared to reference climate conditions (25oC and 385 ppm )

!-------------------------------------------------------------------------------------------------------------------------------
!pass thru variables
real(dp), intent (in) :: ftg0                       !the intercept of the stomatal conductance equation
real(dp), intent (in) :: ftg1                       !the slope of the stomatal conductance equation 
integer,intent (in) :: gs_func                    !flag to indicate which stomatal conductance function to use 

!-------------------------------------------------------------------------------------------------------------------------------
!parameters
real(dp), parameter :: Cb = 1.78                      ! nitrogen use effiency for choloraphyll for light capture, see Evans 1989  
real(dp), parameter :: Cv = 1.2d-5 * 3600.            ! conversion factor from umol CO2 to g carbon
real(dp), parameter :: Jmaxb0 = 0.0311                ! the baseline proportion of nitrogen allocated for electron transport (J)     
real(dp), parameter :: Kc25 = 40.49                   ! Mechalis constant of CO2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
real(dp), parameter :: Ko25 = 27840.                  ! Mechalis constant of O2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
real(dp), parameter :: Cp25 = 4.275                   ! CO2 compensation point at 25C (Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
!-------------------------------------------------------------------------------------------------------------------------------------------------       

theta_cj    = 0.95
theta       = 0.292 / (1.0 + 0.076 / (Nlc * Cb))
ELTRNabsorb = theta  * PARi10
Jmaxb0act   = Jmaxb0 * FNCa * Fj
! adjusted for SDGVM units molm-2s-1
Jmaxb0act   = Jmaxb0act * 1d-6
Jmax        = Jmaxb0act + JmaxCoef * ELTRNabsorb
JmaxL       = theta  * PARimx10 / ( sqrt(1.0 + (theta * PARimx10 / Jmax)**2.0) )        
NUEchg      = (NUEc / NUEcref) * (NUEjref / NUEj)
Wc2Wj       = Wc2Wjb0 * (NUEchg**0.5)
Wc2Wj       = min(1.0, Wc2Wj)
Vcmax       = Wc2Wj * JmaxL * Kj2Kc
     
!added the two below lines for stability when not calculating LUNA everyday
! - in the lowest canopy layers the solver sometimes returns negative values 
! - and with larger intervals between LUNA calcultions (~10 days) the max proportional change can equal 1
! - allowing returned Vcmax and Jmax to be below 0
Vcmax       = max(1d-7,Vcmax)
Jmax        = max(1d-7,Jmax)
     
JmeanL      = theta * PARi10 / ( sqrt(1.0 + (ELTRNabsorb / Jmax)**2.0) )
 
if(KcKjFlag.eq.0) then      !update the Kc,Kj, anc ci information
 
! From SDGVM - for consistency
! in LUNA ko, kc, and tau are at the 10-day mean leaf temp
  k_c = Kc25 * exp((79430.0 / (8.31 * &
 298.15)) * (1.0 - (298.15) / (273.15 + tleafd10)))
  k_o = Ko25 * exp((36380.0 / (8.31 * &
 298.15)) * (1.0 - (298.15) / (273.15 + tleafd10)))
  c_p = Cp25 * exp((37830.0 / (8.31 * &
 298.15)) * (1.0 - (298.15) / (273.15 + tleafd10)))
  !print*, 'LUNA  orig.', k_c,k_o,c_p,tleafd10
 
  !k_c  = exp(35.8   - 80.5 / (0.00831*(tleafd10+273.15)) )
  !k_o  = exp(9.6    - 14.51/ (0.00831*(tleafd10+273.15)) )
  !k_o  = k_o * 1.0d3
  !tau  = exp(-3.949 + 28.99/ (0.00831*(tleafd10+273.15)) )
  !c_p  = 0.5*O2a10/tau
  !print*, 'SDGVM orig.', k_c,k_o,c_p,tleafd10
 
! this function is over-parameterised due to an old attempt to use the Brent solver that is currently not implemented
! Vcmax and J are also assumed to be in molm-2s-1 in SDGVM
! ga and gs must be in molm-2s-1 - gs is a local variable calculated during the call to assimilation_calc
! rd needs to be unpacked, ko, kc, and tau are also needed here 
! z is the LAI layer that is currently under calculation - in SDGVM this only triggers print to screen for LAI layer 1        
! call Photosynthesis_luna(forc_pbot10, tleafd10, relh10, CO2a10, O2a10,rb10, Vcmax, JmeanL, ci, Kc, Kj, A) 
! CALL ASSIMILATION_CALC(fshade,fsunlit,
!&vmx,xvmax,rd,jshade,jsunlit,upt,qshade,
!&qsunlit,t,rn,soil2g,wtwp,ga,rh,C3,kg,ko,kc,tau,p,oi,ca,
!&a,gs,ci,day,mnth,oday,omnth,.FALSE.,gs_func,ftg0,ftg1,i)
 
! should rd be zero below? seems like maybe it should according to the code but by not including rd in the a,ci,gs solution this will underestimate ci
! the below function calls the SDGVM calculation of assimilation, simultaneously solving A, gs, & ci
! for the purposes of the LUNA model A (umolm-2s-1)  and ci (Pa) are returned and used in further calculations
! inputs to this model are in Pa units and molm-2s-1 
!  CALL ASSIMILATION_CALC(0.0,1.0, &
! Vcmax,0.0,0.0,0.0,JmeanL,0.0,0.0, &
! 0.0,tleafd10,1.0,1.0,0.0,rb10*1d2,relh10,1,1.0, &
! k_o,k_c,tau,forc_pbot10,O2a10,CO2a10,A,gs,ci,1,2,1,1,.FALSE., &
! gs_func,ftg0,ftg1,2)

  awc = k_c * (1.0 + O2a10 / k_o)
  Kj  = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p)
  Kc  = max(ci - c_p, 0.0) / (ci + awc)
 
! below else changed to endif
! different from original LUNA code as A in LUNA is smoothed
! between Wc and Wj, and parameters have been calibrated to this
! smoothing. Somewhat decouples LUNA from the SDGVM photosynthesis
! method which assumes the straight minimum of wc or wj but is
! probably a reasonable compromise 
!else
endif
 
! as described above, satisfy the LUNA assumption of co-limitation
! between wc and wj
Wc = Kc * Vcmax
Wj = Kj * JmeanL
!A  = (1.0 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj) 
! Vcmax and J are in molm-2s-1 because that's what SDGVM expects, so the above calculation gives Kc and Kj in molm-2s-1
! the below smoothing function is used for numercial stability,
! at the solution Wj and Wc ought to be pretty similar due to the
! encoded co-ordination hypothesis, therefore this smoothing
! should not have a big effect on A calculated at the optimum
A  = ( (1.0 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj) ) * 1e6
 
!endif

! Cv converts from umolm-2s-1 to gC m-2 hour-1
PSN        = Cv * A * hourpd
!      Vcmaxnight = VcmxTKattge(tair10, tleafn10) / 
!     &VcmxTKattge(tair10, tleafd10) * Vcmax
! SDGVM assumes mean 24 hr air temp 
Vcmaxnight = Vcmax
! again these calculations are adjusted to suit SDGVM molm-2s-1 units for Vcmax and Jmax 
RESP       = Cv * 0.015 * ( Vcmax * hourpd + Vcmaxnight * (24. - hourpd) ) * 1e6                
Net        = Jmax  * 1e6 / Fj
Ncb        = Vcmax * 1e6 / Fc
Nresp      = RESP  / NUEr
     
!print*, A,RESP,Vcmax,Jmax
 
end subroutine LUNA_Ninvestments
!-------------------------------------------------------------------------------------------------------------------------------------------------       



!-------------------------------------------------------------------------------------------------------------------------------------------------       
!Calculate the Nitrogen use effieciency dependence on CO2 and leaf temperature
subroutine NUE(o2a,ci,tgrow,tleaf,NUEj,NUEc,Kj2Kc, &
 ttype,ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ)

! this function can be used to calculate NUEref aswell, it just needs to be called with both temps at 25oC and a ci that maintains a ci:ca of 0.7

! uses module constants of Fc25, Fj25, tfrz, rgas

implicit none
real(dp), intent (in) :: o2a                        !air O2 partial presuure (Pa)
real(dp), intent (in) :: ci                         !leaf inter-cellular [CO2] (Pa) - originally incorrectly labelled as PPM
real(dp), intent (in) :: tgrow                      !10 day growth temperature (oC), 24 hour mean temperature
real(dp), intent (in) :: tleaf                      !leaf temperature (oC)
real(dp), intent (in) :: ftToptV                    !optimum temp for vcmax (oC)
real(dp), intent (in) :: ftHaV                      !activation energy for vcmax
real(dp), intent (in) :: ftHdV                      !deactivation energy for vcmax
real(dp), intent (in) :: ftToptJ                    !optimum temp for jmax (oC)
real(dp), intent (in) :: ftHaJ                      !activation energy for jmax
real(dp), intent (in) :: ftHdJ                      !deactivation energy for jmax
!integer, intent (in):: ttype                      !temperature scaling method to use
integer:: ttype                      !temperature scaling method to use
real(dp), intent (out):: NUEj                       !nitrogen use efficiency for electron transport under refernce environmental conditions (25oC and 385 ppm co2)
real(dp), intent (out):: NUEc                       !nitrogen use efficiency for carboxylation under reference environmental conditions  (25oC and 385 ppm co2)
real(dp), intent (out):: Kj2Kc                      !the ratio of Kj to Kc 
!------------------------------------------------
!intermediate variables
real(dp) :: Fj                                      !the temperature adjusted factor for Jmax 
real(dp) :: Fc                                      !the temperature adjusted factor for Vcmax 
real(dp) :: Kc                                      !conversion factor from Vcmax to Wc 
real(dp) :: Kj                                      !conversion factor from J to W 
real(dp) :: k_o                                     !Rubsico O2 specifity
real(dp) :: k_c                                     !Rubsico CO2 specifity
real(dp) :: awc                                     !second deminator term for rubsico limited carboxylation rate based on Farquhar model
real(dp) :: c_p                                     !CO2 compenstation point (Pa)
real(dp) :: tau
!real(dp) :: t_scalar                                !temperature scaling factor for Vcmax or Jmax 

real(dp), parameter :: Fc25 = 294.2               ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
real(dp), parameter :: Fj25 = 1257.0              ! Fj25 = 8.06*156  #see COSTE 2005 and Xu et al 2012
real(dp), parameter :: Kc25 = 40.49               ! Mechanis constant of CO2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
real(dp), parameter :: Ko25 = 27840.              ! Mechanis constant of O2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
real(dp), parameter :: Cp25 = 4.275               ! CO2 compensation point at 25C (Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
   
!-------------------------------------------------------------------------------------------------------------------------------------------------       
!print*, 'LUNA NUE, ttype:', ttype
!Fc  = VcmxTKattge(tgrow, tleaf) * Fc25
!Fj  = JmxTKattge(tgrow, tleaf)  * Fj25
Fc  = T_SCALAR(tleaf,'v',ttype,ftToptV,ftHaV,ftHdV,tgrow,.FALSE.)*Fc25
Fj  = T_SCALAR(tleaf,'j',ttype,ftToptJ,ftHaJ,ftHdJ,tgrow,.FALSE.)*Fj25
 
k_c = Kc25 * exp((79430.0 / (8.31 * &
 298.15)) * (1.0 - (298.15) / (273.15 + tleaf)))
k_o = Ko25 * exp((36380.0 / (8.31 * &
 298.15)) * (1.0 - (298.15) / (273.15 + tleaf)))
c_p = Cp25 * exp((37830.0 / (8.31 * &
 298.15)) * (1.0 - (298.15) / (273.15 + tleaf)))
!print*, 'LUNA  orig.', k_c,k_o,c_p,tleaf
 
! From SDGVM - for consistency
! in LUNA ko, kc, and tau are at the 10-day mean leaf temp
! in SDGVM these are all in Pa
!k_c = exp(35.8   - 80.5 /(8.31d-3*(tleaf+273.15)))
!k_o = exp(9.6    - 14.51/(8.31d-3*(tleaf+273.15)))
!k_o = k_o*1.d3
!tau = exp(-3.949 + 28.99/(8.31d-3*(tleaf+273.15)))
!c_p = 0.5*o2a/tau
!print*, 'SDGVM orig.', k_c,k_o,c_p
     
awc = k_c * ( 1.0 + o2a/k_o )
Kj  = max( ci-c_p,0.0 ) / ( 4.0*ci + 8.0*c_p )
Kc  = max( ci-c_p,0.0 ) / ( ci+awc )
     
NUEj  = Kj * Fj
NUEc  = Kc * Fc  
Kj2Kc = Kj / Kc
     
end subroutine NUE
!-------------------------------------------------------------------------------------------------------------------------------------------------       



!-------------------------------------------------------------------------------------------------------------------------------------------------       
!Calculate the temperature response for respiration, following Bernacchi PCE 2001
function RespTBernacchi(tleaf)

implicit none

real(dp) RespTBernacchi
real(dp), intent(in):: tleaf  !leaf temperature (oC)

!RespTBernacchi= exp(18.72-46.39/(rgas*1.e-6 *
!&(tleaf+tfrz)))
RespTBernacchi= exp( 18.72-46.39/(8.31d-3 * &
 (tleaf+273.15)) )

end function RespTBernacchi
!-------------------------------------------------------------------------------------------------------------------------------------------------       

     

     

!-------------------------------------------------------------------------------------------------------------------------------------------------       
subroutine LUNA_nogrowth(max_daily_pchg,vcmx25_z,jmx25_z,enzs_z) 

real(dp), intent (in)    :: max_daily_pchg      ! maximum daily percentrage change for nitrogen allocation
real(dp), intent (inout) :: vcmx25_z            ! leaf Vc,max25 (umol/m2 leaf/s) for canopy layer 
real(dp), intent (inout) :: jmx25_z             ! leaf Jmax25 (umol electron/m**2/s) for canopy layer
real(dp), intent (inout) :: enzs_z              ! enzyme decay status 1.0-fully active; 0-all decayed during stress
real(dp)                 :: max_daily_decay     ! maximum daily percentrage change for nitrogen allocation
!INTEGER                :: nrad                ! number of canopy layers
!INTEGER                :: z                   ! canopy layer index

! need to determine if these really need to have each layer or are passed back to each layer in npp calc
!-------------------------------------------------------------------------------------------------------------------------------------------------       
!nrad            = 12

!assume enzyme turnover under maintenance is 10 times lower than enzyme change under growth
max_daily_decay = min(0.5, 0.1 * max_daily_pchg)

!print*, 'LUNA nogrowth:', max_daily_decay
!print*, vcmx25_z(1)

!do z = 1 , nrad
 !decay is set at only 50% of original enzyme in view that plant will need to maintain their basic functionality
 if(enzs_z>0.5) then 
  enzs_z   = enzs_z   * (1.0 - max_daily_decay)
  jmx25_z  = jmx25_z  * (1.0 - max_daily_decay) 
  vcmx25_z = vcmx25_z * (1.0 - max_daily_decay) 
 endif
!end do              

!print*, vcmx25_z(1)

end subroutine LUNA_nogrowth 
!-------------------------------------------------------------------------------------------------------------------------------------------------       

end module luna_methods


