! ==================================================================== !
! Name        : sdgvm.f90
! Author      : Mark Lomas
! Version     : 0.1 Comments
! Copyright   : Your copyright notice
! Description : sdgvm
! ==================================================================== !

program sdgvm

use real_precision,  only: dp
use dims,            only: max_cohorts, max_years, max_outputs, str_len
use system_state,    only: ssv
use pft_parameters
use state_methods
use output_methods
use input_methods
use phenology_methods
use daily_step
use hydrology_methods
use light_methods
use soil_methods
use misc_values
use crops
use file_class
use file_object
use fnames_class
use fnames_object
use input_file

implicit none

character(len=str_len), parameter :: stver = 'VERSION current'
real(dp) :: defaulttopsl
parameter(defaulttopsl = 5.0)
logical :: check_closure
parameter(check_closure = .false.)

real(dp), dimension(max_cohorts) :: lai,evt,sresp,rof,gpp,ftprop, &
 nppstoreold,trn,lch,bioo,ht,soilc,soiln,minn,ftcov,covo, &
 leafnpp,stemnpp,rootnpp,bioleaf

real(dp) :: ca,resp,soilt,grassrc,tmp(12,31),prc(12,31), &
 hum(12,31),cld(12),latdel,londel,leafper,stemper,rootper,avnpp,avgpp, &
 avlai,co20,co2f,avrof,infix,avnppst,sum1,oscale,yield,co2(max_years), &
 sumcov,maxcov,maxbio,barerc,avtrn,firec,pet2,f2,f3,avevt,sumbio,kd, &
 kx,stembio,rootbio,sum,lutab(255,100),awl(4),nci(4),gsn,eemm,etmm,rn, &
 xlatf,xlatres,xlon0,xlonres,lat_lon(max_sites,2),nleaf,canga,s1in,  &
 dp2,h2o,adp(4),sfc(4),sw(4),sswc(4),nupc,swc,swf,ssm,xfprob, &
 xno_fires,nfix,daygpp,evap,tran,roff,pet, &
 ans(12,31),interc,evbs,soil_chr(10),cluse(max_cohorts,max_years), &
 cirr(max_cohorts,max_years),cfert(max_cohorts,max_years,3), &
 soil_chr2(10),fprob,q,qdirect,qdiff,topsl,hrs,fpr,resp_s,resp_r, &
 lmor_sc(3600,max_cohorts),leaflitter,laiinc,daynpp,resp_l,resp_m, &
 pm_tmp,pm_prc,pm_hum,pm_evbs,pm_swc,pm_swf,pm_ssm,cos_zen, &
 daily_out(max_outputs,max_cohorts,12,31)

real(dp), dimension(30,12,max_cohorts) :: ce_light,ce_ci,ce_ga, &
 ce_maxlight
real(dp), dimension(30,max_cohorts) :: ce_t,ce_rh

integer :: sites,yr0,yrf,snp_no,snpshts(max_years),snp_year,day,d, &
 isite,du,otagsn(max_outputs),otagsnft(max_outputs),nft,nat_map(8), &
 ilanduse,siteno,iofn, iofnft,iofngft,nomdos,i,j,k,ft,site,year, &
 covind,bioind,xyearf,mnth,fireres,xyear0,omav(max_outputs),year_out, &
 luse(max_years),fno,iyear,oymd,oymdft,iargc,yearind(max_years),idum, &
 outyears,thty_dys,xlatresn,xlonresn,day_mnth,yearv(max_years),nyears, &
 narg,seed1,seed2,seed3,spinl,xseed1,site_dat,site_out, &
 country_id,outyears1,outyears2,budo(max_cohorts),seno(max_cohorts), &
 sit_grd,co,nn1,fid,kode,imap,par_loops,gs_type,nco

real(dp), dimension(500,12,31) :: xtmpv,xprcv,xhumv,xswrv
real(dp), dimension(500,12) :: xcldv
real(dp), dimension(12,31) :: swr

real(dp) :: soilCtoN,soilp,tleaf_n,tleaf_p


integer :: daily_co2,read_par,hw_j,phen_cor,soilcn_map,t,total_t
integer :: fid_output_state
logical :: fire_ant(max_years)
logical :: harvest_ant(max_years)
real(dp) :: swrv(500,12,31)

character(len=str_len) :: st1
character(len=str_len), dimension(max_outputs) :: &
 ofmt_daily,ofmt_monthly,ofmt_yearly,otags
character(len=str_len) :: country_name
character(len=80) :: buff1

character(len=20) :: sttemp

logical :: speedc,crand,xspeedc,withcloudcover,l_clim,l_lu,l_soil(20)
logical :: l_stats,l_regional,out_cov,out_bio,out_bud,out_sen,l_b_and_c

logical :: met_seq,goudriaan_old
!----------------------------------------------------------------------!

!call fnms%set_names()

!----------------------------------------------------------------------!
! Some temporary settings from ants work.                              !
!----------------------------------------------------------------------!
daily_co2 = -1
read_par = -1
par_loops = 10
gs_type = 1

ce_t           = 15.0
ce_rh          = 75.0 
ce_ci          = 28.0
ce_ga          = 0.10
ce_light       = 5.0e-4 
ce_maxlight    = 1.0e-3 


!----------------------------------------------------------------------!
! in b_input_methods.f90
! Get input filename from the command line.                            !
! buff1 is the path to the input file name                             !
!----------------------------------------------------------------------!
call get_input_filename(buff1)

!----------------------------------------------------------------------!
! Reads the input file.inp is a variable of Type Input declared in     !
! input_file.f90.variable inp holds other Derived Data Type Variables  !
! with the values read from the input file organized with them.        !
!----------------------------------------------------------------------!
call inp%read_input_file(trim(buff1))

!----------------------------------------------------------------------!
! in b_input_methods.f90                                               !
! Check the number of command line arguments are consistent with the   !
! input file.You can pass some variables from the command line or the  !
! input file e.g. output directory.This checks that they are all there !
!----------------------------------------------------------------------!
call command_line_argument_check()

call read_param(stver)
!----------------------------------------------------------------------!
! ! in b_input_methods.f90                                             !
! Process the input file data.                                         !
!----------------------------------------------------------------------!
call process_input_file(buff1,xlatf,xlon0,xlatres,xlonres, &
 speedc,xspeedc,crand,outyears,nyears,yr0,yrf,yearind,idum,yearv, &
 nomdos,otags,omav,ofmt_daily,ofmt_monthly,ofmt_yearly,outyears1, &
 outyears2,oymd,otagsn,otagsnft,snpshts,snp_no,out_cov, &
 out_bio,out_bud,out_sen,lutab,grassrc,barerc,fireres,luse,l_b_and_c, &
 soil_chr,topsl,defaulttopsl,sites,latdel,londel,lat_lon,day_mnth, &
 thty_dys,xlatresn,xlonresn,ilanduse,nft,xyear0,xyearf,lmor_sc, &
 oymdft,iofnft,sit_grd,du,narg,fire_ant,harvest_ant,met_seq,par_loops)
!----------------------------------------------------------------------!
! in b_input_methods.f90                                               ! 
! Read internal parameters from "param.dat" file, and io               !
! parameters from "misc_params.dat".                                   !
!----------------------------------------------------------------------!
!call read_param(stver)

!----------------------------------------------------------------------!
! in b_output_methods.f90
! set up mapping for output files for and place in inp%output%*%map.   !
!----------------------------------------------------------------------!
call output_mapping(nomdos,otags)

!----------------------------------------------------------------------!
! in open_files.f90                                                    !
! Opens default output files.                                          !
!----------------------------------------------------------------------!
call open_default()

!----------------------------------------------------------------------!
! in open_files.f90                                                    !
! Opens site_info.dat file                                             !
!----------------------------------------------------------------------!
call open_site_info()

if (inp%run%read_in_state) then
  call fun%openu_old(trim(inp%dirs%input)//'/state.dat',fid)
endif

if (inp%run%output_state) then
  call fun%openu(trim(inp%dirs%output)//'/state.dat',fid_output_state)
endif

call crop_outputs(nft,0)

if (outyears2>0) call open_optional(out_cov,out_bio,out_bud,out_sen)

!call open_state()

!----------------------------------------------------------------------!
! in open_files.f90                                                    !
! Opens state input files                                              !
!----------------------------------------------------------------------! 
call open_tile_pft()

if(inp%run%gs_func==3) goudriaan_old = .TRUE.

if (inp%run%s070607) then
  inp%run%s140129 = .false.
  inp%run%read_daily_co2 = .false.
  inp%run%read_par = .false.
  inp%run%subdaily = .false.
  inp%run%read_clump = .false.
  inp%run%calc_zen = .false.
  inp%run%no_soil_water_limitation = .false.

  inp%run%cstype     = 0  
  inp%run%ncalc_type = 0
  inp%run%ttype      = 0
  inp%run%vcmax_type = 0
  inp%run%soilp_map  = 0
  inp%run%gs_func    = 0

  soilcn_map = 0
  phen_cor   = 0 
  hw_j       = 0
endif

if (inp%run%s140129) then
  inp%run%s070607 = .false.
  inp%run%read_daily_co2 = .false.
  inp%run%read_par = .true.
  inp%run%subdaily = .true.
  inp%run%read_clump = .false.
  inp%run%calc_zen = .true.
  inp%run%no_soil_water_limitation = .false.

  inp%run%cstype     = 0
  inp%run%ncalc_type = 1
  inp%run%ttype      = 0
  inp%run%vcmax_type = 1
  inp%run%soilp_map  = 0
  inp%run%gs_func    = 1

  soilcn_map = 0
  phen_cor   = 1
  hw_j       = 0
endif

if(goudriaan_old) hw_j = 3

if (inp%run%s070607) then
  tgp%p_et = 0.7
  tgp%p_pet = 0.7
endif
if (inp%run%s140129) then
  tgp%p_et = 1.0
  tgp%p_pet = 1.0
endif

if (.not.inp%run%read_clump) then
  pft_tab%can_clump = 1.0
endif

!----------------------------------------------------------------------!
! in data.f90
! Read co2 file.                                                      !
!----------------------------------------------------------------------!
call readco2(yr0,yrf,co2,daily_co2,nyears)
!----------------------------------------------------------------------!

call set_goud_params()

site_dat = 0

!======================================================================!
!                         SITE LOOP                                    !
!======================================================================!
do site=1,sites
  
  ssp%jday = 5000
  speedc = xspeedc
  ssp%lat=lat_lon(site,1)
  ssp%lon=lat_lon(site,2)
  seed1 = inp%run%random_seed
  seed2 = 2*seed1
  seed3 = 3*seed1

!----------------------------------------------------------------------!
! Write lat/lon in crop output file                                    !
!----------------------------------------------------------------------!
  call crop_outputs(nft,2)

!----------------------------------------------------------------------!
! Read in climate.                                                     !
!----------------------------------------------------------------------!
  call read_climate(ssp%lat,ssp%lon,xlatf, &
 xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv, &
 xcldv,xswrv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats, &
 siteno,day_mnth,thty_dys,sit_grd,withcloudcover)

!----------------------------------------------------------------------!
! Read in soil parameters.                                             !
!----------------------------------------------------------------------!
  call read_soil(ssp%lat,ssp%lon,soil_chr,soil_chr2,du,l_soil)
  soilCtoN = soil_chr2(9)
  soilp = soil_chr2(10)

!----------------------------------------------------------------------!
! Read in canopy clumping index from a map                             !
!----------------------------------------------------------------------!
!      IF(read_clump==0)THEN
!      !if no clumping set canopy clumping index to 1
!        ftcan_clump(:) = 1
!      ELSEIF(read_clump==2)THEN
!        CALL EX_CLUMP(lat,lon,map_clump,du)
!        ftcan_clump(:) = map_clump
!      ENDIF

!----------------------------------------------------------------------!
! in data.f90                                                          !
! Read land use.cluse(nft,year) the fraction of each ft per year       !
!----------------------------------------------------------------------!
  call read_landuse(ilanduse,yr0,yrf,du,nft,ssp%lat,ssp%lon,lutab,luse,&
 cluse,l_lu)
  
!----------------------------------------------------------------------!
! in crops.f90                                                         !
! Read fertilizers.cfert(nft,year,NPK).Units are kg/ha,scale is 0.1    !
!----------------------------------------------------------------------!
  call read_fertilizers(du,yr0,yrf,ssp%lat,ssp%lon,nft,cfert)

!----------------------------------------------------------------------!
! in crops.f90                                                         !
! Read irrigation.cirr(nft,year).Units are %
!----------------------------------------------------------------------!
  call read_irrigation(du,yr0,yrf,ssp%lat,ssp%lon,nft,cirr)
   
!======================================================================!
!                         DRIVING DATA IF                              !
!======================================================================!
  if ((l_clim).and.(l_stats).and.(l_soil(1)).and.(l_soil(3)).and. &
 (l_soil(5)).and.(l_soil(8)).and.(l_lu)) then

    site_dat = site_dat + 1
       
!----------------------------------------------------------------------!
! in b_output_methods
! Write lat & lon for output files.                                    !
!----------------------------------------------------------------------!
    call write_lat_lon(ssp%lat,ssp%lon)

!----------------------------------------------------------------------!
! Computation of hydrological parameters.                              !
!----------------------------------------------------------------------!
    call wsparam(l_b_and_c,nupc,awl,kd,kx,nci,infix,adp,topsl,sfc,sw,sswc)

!----------------------------------------------------------------------!
! in data.f90
! Extract initial and final C02 values.                                !
!----------------------------------------------------------------------!
    call co2_0_f(co20,co2f,yearv,yr0,co2,nyears)
    
!----------------------------------------------------------------------!
! in data.f90
! Extract the country/state corresponding to the site.                 !
!----------------------------------------------------------------------!
    call country(ssp%lat,ssp%lon,country_name,country_id,l_regional)

!----------------------------------------------------------------------!
! Write site info to 'site_info.dat'.                                  !
!----------------------------------------------------------------------!
    fid = fun%get_id('site_info.dat')
    write(fid,'(a15,i6,f9.3,f9.3,1x,2f6.1,1x,2f7.1,f7.3,1x,3f7.3,f7.1,1x,200(2f6.1,1x))') &
 country_name,country_id,ssp%lat,ssp%lon,co20,co2f,ssp%sand,ssp%silt,ssp%bulk,ssp%wilt, &
 ssp%field,ssp%sat,ssp%soil_depth, &
 (cluse(ft,yearv(1)-yr0+1),cluse(ft,yearv(nyears)-yr0+1),ft=1,nft)

!----------------------------------------------------------------------!
! Initialise the system state.                                         !
!----------------------------------------------------------------------!
    call initialise_state(nft,cluse,xtmpv,soilt)

    if (mod(site_dat,max(site_out,1))==min(1,site_out)-1) &
 write(*,'( '' Site no. '',i0,'', Lat ='',f6.2,'', Lon ='',f7.2,'' Cohorts = '',i0)') site_dat,ssp%lat,ssp%lon,ssp%cohorts

!======================================================================!
!                         YEARLY LOOP                                  !
!======================================================================!
    do iyear=1,nyears

      year = yearv(iyear)

      !The indexed year 
      ssp%iyear = iyear
      !The calendar year
      ssp%year = year
           
      nfix = infix

!----------------------------------------------------------------------!
! in data.f90
! Set CO2 value 'ca' from 'co2' or 'co2const'.                         !
!----------------------------------------------------------------------!
      call set_co2(ca,iyear,speedc,co2,year,yr0)

      if (mod(iyear,max(sop%year_out,1))==min(1,sop%year_out)-1) then
        write(*,'('' Year no. '',i3,'' '',i4,'', ca = '',f6.2,'', cohorts = '', i0,''.'')') &
 iyear,year,ca,ssp%cohorts
      endif

!----------------------------------------------------------------------!
! in data.f90                                                          !
! Set 'tmp' 'hum' 'prc' 'swr' and 'cld', and calc monthly and yearly av!
! in crops.f90                                                         !
! Set sow day pft_tab(3:nft)%sowday and find an estimate of required   !
! GDD for the grid cell and crop pft_tab(ft)%cropgdd(1,:)              !
!----------------------------------------------------------------------!
      do nn1=1,0,-1
        call set_climate(xtmpv,xprcv,xhumv,xcldv,xswrv,withcloudcover, &
 yearv,iyear,tmp,prc,hum,cld,swr,thty_dys,yr0,year,nyears,nn1)
        call seasonality(tmp,prc,cld,thty_dys,nft,year,nn1)
      enddo
      
      do ft=1,nft
        !Irrigation in fraction of gridcell that is irrigated per crop
        pft_tab(ft)%irrig(3)=0.01*cirr(ft,year-yr0+1)
        !Nitrogen in kg/ha
        pft_tab(ft)%fert(1)=10*cfert(ft,year-yr0+1,1)
        !Phosphorus in kg/ha
        pft_tab(ft)%fert(2)=10*cfert(ft,year-yr0+1,2)
        !Potassium in kg/ha
        pft_tab(ft)%fert(3)=10*cfert(ft,year-yr0+1,3)
      enddo
      
!----------------------------------------------------------------------!
! in crops.f90                                                         !
! Figures the optimal crop lai pft_tab(ft)%optlai based on the applied !
! fertilizer                                                           !
!----------------------------------------------------------------------!
      call fert_crops(nft)

!      call read_opt_par()

!----------------------------------------------------------------------!
! in data.f90                                                          ! 
! Set land use through ftprop.ftprop(nft) assumes the values of        !
! cluse(nft,year) with the % for each class                            !
!----------------------------------------------------------------------!
      call set_landuse(ftprop,tmp,prc,nat_map,nft,cluse,year,yr0)
    
!----------------------------------------------------------------------!
! in veg_dynamics.f90                                                  !
! Calculates ftprop(nft) which now becomes the fraction of cover that  !
! must be added for each pft                                           !
!----------------------------------------------------------------------!      
      call cover(nft,tmp,prc,firec,fireres,fprob,ftprop,check_closure)

      call initialise_new_cohorts(nft,ftprop,check_closure)
 
!----------------------------------------------------------------------!
! in soil_methods.f90                                                  !
!----------------------------------------------------------------------!
      call mkdlit()

!----------------------------------------------------------------------!
! in state_methods.f90                                                 !
!----------------------------------------------------------------------!
      call restrict_cohort_numbers()
      
!----------------------------------------------------------------------!
! Initialisations that were in doly at the beginning of the year       !
!----------------------------------------------------------------------!
      do ft=1,ssp%cohorts
        leafnpp(ft) = 0.0
        stemnpp(ft) = 0.0
        rootnpp(ft) = 0.0

        ssv(ft)%evp = 0.0

        budo(ft) = 0
        seno(ft) = 0
        
!----------------------------------------------------------------------!
! Height (m) ccn  ht(ft) = 0.807*(laimax(ft)**2.13655)                 !
!----------------------------------------------------------------------!
        ht(ft) = 0.807*(5.0**2.13655)
        if (ht(ft)>50.0)  ht(ft)=50.0

      enddo
      
!======================================================================!
!                         MONTHLY LOOP                                 !
!======================================================================!
      do mnth=1,12
        
        ssp%mnth = mnth
!----------------------------------------------------------------------!
! in soil_methods.f90                                                  !
!----------------------------------------------------------------------!
        call sum_soilcn(soilc,soiln,minn)

!======================================================================!
!                         DAILY LOOP                                   !
!======================================================================!
        do day=1,no_days(year,mnth,thty_dys)

          !Day of the month
          ssp%day = day
          !Counter that starts at 5000 for each site
          ssp%jday = ssp%jday + 1
          
          fpr=0.0

!----------------------------------------------------------------------!
! Updata daily temperature memory.                                     !
!----------------------------------------------------------------------!
          do i=1,299
            ssp%tmem(301-i) = ssp%tmem(300-i)
          enddo
          ssp%tmem(1) = tmp(mnth,day)
          
!----------------------------------------------------------------------!
! Radiation calculation.                                               !
!----------------------------------------------------------------------!
          hrs = dayl(ssp%lat,no_day(year,mnth,day,thty_dys))
!          call PFD(lat,no_day(year,mnth,day,thty_dys),hrs,cld(mnth), &
!     qdirect,qdiff,q)

!----------------------------------------------------------------------!
! in light_methods.f90                                                 !
! Calculates direct(qdirect),diffused(qdiff) and total(q) PAR radiation!
! for day of the month in mol/m2/sec.                                  !
!----------------------------------------------------------------------!
          call pfd_ant(ssp%lat,no_day(year,mnth,day,thty_dys),hrs, &
 cld(mnth),qdirect,qdiff,q,swr(mnth,day),inp%run%read_par,.false., &
 t,total_t,inp%run%calc_zen,cos_zen)
 
!----------------------------------------------------------------------!
! Mix water resources.                                                 !
!----------------------------------------------------------------------!
!          call MIX_WATER(ftcov,nft)

!======================================================================!
!                         COHORT LOOP                                  !
!======================================================================!
          do ft=1,ssp%cohorts
            fpr=0.0
            ssp%cohort = ft

            if (ssv(ft)%cov>0.0) then
              
              !Provisional nfert value linked to N fertilizer
              if(pft(ft)%fert(1)<=0.0) then
                nfix=0.5
              else
                nfix=0.1*pft(ft)%fert(1)          
              endIF
              
!----------------------------------------------------------------------!
! in state_methods.f90                                                 !
! Calculates leaf molecular weight msv%mv_leafmol based on sla         !
! and msv%mv_respref which is a variable used in respiration           !
! calculations.                                                        !
!----------------------------------------------------------------------!
              call set_misc_values(pft(ft)%sla,tmp(mnth,day))
              
!----------------------------------------------------------------------!
! nppstore leafnpp stemnpp rootnpp leaflit stemlit rootlit in mols
!----------------------------------------------------------------------!
              soilt = 0.97*soilt + 0.03*tmp(mnth,day) !unused
!----------------------------------------------------------------------!
! in crops.f90                                                         !
!----------------------------------------------------------------------!            
              call irrigate(ssp%cohort,sfc,sw) 

!----------------------------------------------------------------------!
! in daily_step.f90                                                    !
! Outputs daygpp and resp_l which are daily gpp and leaf respiration.  !
! Also canga and gsn which are the boundary and canopy conductance     !
! which are used in evapotranspiration sub below.                      !
!----------------------------------------------------------------------!              
              call dailyStep(tmp(mnth,day),prc(mnth,day),hum(mnth,day), &
 cld(mnth),ca,soilc(ft),soiln(ft),minn(ft),kd,kx,daygpp,resp_l,lai(ft), &
 ht(ft),ft,nleaf,hrs,q,qdirect,qdiff,fpr,tleaf_n,tleaf_p,canga,gsn,rn, &
 ce_light(:,:,ft),ce_ci(:,:,ft),ce_t,ce_maxlight(:,:,ft),ce_ga(:,:,ft), &
 ce_rh,check_closure,par_loops,ssp%lat,year,mnth,day,thty_dys, &
 inp%run%gs_func,swr(mnth,day))

!----------------------------------------------------------------------!
! in daily_step.f90                                                    !
! Calculates evapotransiration etmm and evaporation eemm.              !
!----------------------------------------------------------------------!
              call evapotranspiration(tmp(mnth,day),hum(mnth,day),rn,canga,gsn,hrs,eemm,etmm)
              pet = eemm
              pet2 = pet

!----------------------------------------------------------------------!
! in hydrology_methods.f90                                             !
!----------------------------------------------------------------------!            
              call hydrology(adp,sfc,sw,sswc,awl,kd,kx,eemm,etmm,pet2,prc(mnth,day), &
     s1in,tmp(mnth,day),ssv(ft)%lai%tot(1),evap,tran,roff,interc,evbs,f2,f3,ft)

!----------------------------------------------------------------------!
! in phenology_methods.f90                                             !
!----------------------------------------------------------------------!            
              call phenology(yield,laiinc)
    
!----------------------------------------------------------------------!
! in phenology_methods.f90                                             !
!----------------------------------------------------------------------!            
              call allocation(laiinc,daygpp,resp_l,lmor_sc(:,pft(ft)%itag),resp, &
     leaflitter,stemnpp(ft),rootnpp(ft),resp_s,resp_r,resp_m,check_closure)

              ssv(ft)%slc = ssv(ft)%slc + leaflitter*ssv(ft)%cov

              call soil_dynamics2(pet,prc(mnth,day),tmp(mnth,day),f2/10.0,f3/10.0,nfix, &
     nci,sresp(ft),lch(ft),ca,site,year,yr0,yrf,speedc,ft,check_closure)
    
!----------------------------------------------------------------------!
              swc = ssv(ft)%soil_h2o(1)+ssv(ft)%soil_h2o(2)+ssv(ft)%soil_h2o(3)+ssv(ft)%soil_h2o(4)
              swf = (swc-sw(1)-sw(2)-sw(3)-sw(4))/(sfc(1)+sfc(2)+sfc(3)+sfc(4)-sw(1)-sw(2)-sw(3)-sw(4))
              ssm = ssv(ft)%soil_h2o(1)/10.0/topsl
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Compute fire inputs and call fire routine on the last day of the     !
! month.                                                               !
!----------------------------------------------------------------------!
!       CALL BFIRE_INPUTS(pm_tmp,pm_prc,pm_hum,pm_evbs,pm_swc,pm_swf,pm_ssm, &
! tmp(mnth,day),prc(mnth,day),hum(mnth,day),evbs,swc,swf,ssm,day)
!       IF ((day == 30).and.(ft == ssp%cohorts)) CALL FIRE_B(xfprob,pm_tmp,pm_prc,pm_hum, &
! pm_evbs,pm_swc,pm_swf,pm_ssm,lat,mnth,ssp%cohort)

!----------------------------------------------------------------------!
! Set daily memories for output.                                       !
!----------------------------------------------------------------------!
              daily_out(1,ft,mnth,day) = ssv(ft)%lai%tot(1)
              daily_out(2,ft,mnth,day) = roff
              daily_out(3,ft,mnth,day) = evap+tran
              daily_out(4,ft,mnth,day) = tran
              daily_out(5,ft,mnth,day) = ssv(ft)%npp
              daily_out(6,ft,mnth,day) = daygpp
              daily_out(7,ft,mnth,day) = sresp(ft)
              daily_out(8,ft,mnth,day) = daily_out(5,ft,mnth,day) - &
 daily_out(7,ft,mnth,day)
              daily_out(9,ft,mnth,day) = tmp(mnth,day)
              daily_out(10,ft,mnth,day) = prc(mnth,day)
              daily_out(11,ft,mnth,day) = hum(mnth,day)
              daily_out(12,ft,mnth,day) = ssv(ft)%nppstore(1)
              daily_out(13,ft,mnth,day) = swf
              daily_out(14,ft,mnth,day) = pet
              daily_out(15,ft,mnth,day) = interc
              daily_out(16,ft,mnth,day) = evbs
              daily_out(17,ft,mnth,day) = &
 min(1.0,ssv(ft)%soil_h2o(1)/10.0/topsl)
              daily_out(18,ft,mnth,day) = swc
              daily_out(19,ft,mnth,day) = 1.0
!        daily_out(19,ft,mnth,day) = resp
              daily_out(20,ft,mnth,day) = qdirect*hrs*3600.0
              daily_out(21,ft,mnth,day) = qdiff*hrs*3600.0
              daily_out(22,ft,mnth,day) = 1.0
!        daily_out(22,ft,mnth,day) = nleaf
!        daily_out(23,ft,mnth,day) = leaflit(ft) - lflitold
              daily_out(24,ft,mnth,day) = cld(mnth)
              daily_out(25,ft,mnth,day) = fpr
              daily_out(26,ft,mnth,day) = 0.0
              if (pft(ft)%sla > 0.0)  daily_out(26,ft,mnth,day) = &
 ssv(ft)%lai%tot(1)*12.0/pft(ft)%sla/25.0
              daily_out(27,ft,mnth,day) = &
 ssv(ft)%bio(1) + ssv(ft)%stem%tot(1) + ssv(ft)%nppstore(1)
              daily_out(28,ft,mnth,day) = ssv(ft)%bio(2) + ssv(ft)%root%tot(1)
              daily_out(29,ft,mnth,day) = soilc(ft)
              daily_out(30,ft,mnth,day) = soiln(ft)
              daily_out(31,ft,mnth,day) = lch(ft)
              daily_out(32,ft,mnth,day) = firec/360.0
              daily_out(33,ft,mnth,day) = yield

!----------------------------------------------------------------------!
              if ((ssv(ft)%bb==day+(mnth-1)*30).and.(budo(ft)==0)) &
 budo(ft) = ssv(ft)%bb
              if ((ssv(ft)%ss==day+(mnth-1)*30).and.(seno(ft)==0)) &
 seno(ft) = ssv(ft)%ss

            endif
!======================================================================!
!                         END OF COVER IF                              !
!======================================================================!

          enddo
!======================================================================!
!                         END OF COHORT LOOP                           !
!======================================================================!

        enddo
!======================================================================!
!                         END OF DAILY                                 !
!======================================================================!

      enddo
!======================================================================!
!                         END OF MONTHLY LOOP                          !
!======================================================================!

      call growth(nft,lai,stembio,rootbio,check_closure)
  
!----------------------------------------------------------------------!
! Average outputs by cover proportions.                                !
!----------------------------------------------------------------------!
      avlai   = 0.0
      avnpp   = 0.0
      avnppst = 0.0
      avrof   = 0.0
      avtrn   = 0.0
      avevt   = 0.0
      avgpp   = 0.0
      do ft=1,ssp%cohorts
        avnppst = avnppst + ssv(ft)%cov*ssv(ft)%nppstore(1)
        avrof   = avrof   + ssv(ft)%cov*rof(ft)
        avtrn   = avtrn   + ssv(ft)%cov*trn(ft)
        avevt   = avevt   + ssv(ft)%cov*evt(ft)
        avgpp   = avgpp   + ssv(ft)%cov*gpp(ft)
      enddo

!----------------------------------------------------------------------!
      sumbio = 0.0
      bioind = 0
      covind = 0
      maxbio = 0.0
      maxcov = 0.0
      leafper = 0.0
      stemper = 0.0
      rootper = 0.0

      do ft=1,nft
        bioo(ft) = 0.0
        covo(ft) = 0.0
        if (ssp%co2ftmap(ft,1) > 0) then
          budo(ft) = ssv(ssp%co2ftmap(ft,2))%bb
          seno(ft) = ssv(ssp%co2ftmap(ft,2))%ss
        else
          budo(ft) = 0
          seno(ft) = 0
        endif
        do j=1,ssp%co2ftmap(ft,1)
          co = ssp%co2ftmap(ft,j+1)
          do mnth=1,12
            do day=1,no_days(year,mnth,thty_dys)
              bioo(ft) = bioo(ft) + (daily_out(26,co,mnth,day) + &
 daily_out(27,co,mnth,day) + daily_out(28,co,mnth,day))*ssv(co)%cov
            enddo
          enddo
          covo(ft) = covo(ft) + ssv(co)%cov
          leafper = leafper + ssv(co)%npl*ssv(co)%cov
          stemper = stemper + ssv(co)%nps*ssv(co)%cov
          rootper = rootper + ssv(co)%npr*ssv(co)%cov
        enddo
        bioo(ft) = bioo(ft)/360.0
        sumbio = sumbio + bioo(ft)
        if (covo(ft)>maxcov) then
          maxcov = covo(ft)
          covind = ft
        endif
        if (bioo(ft)>maxbio) then
          maxbio = bioo(ft)
          bioind = ft
        endif
      enddo

!----------------------------------------------------------------------!
! Write outputs.                                                       !
!----------------------------------------------------------------------!
! General.                                                             !
!----------------------------------------------------------------------!
      if (iyear>=nyears-outyears+1) then
write(fun%get_id('lai.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,1,'Max')
write(fun%get_id('npp.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,5,'Add')
write(fun%get_id('scn.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,29,'Average')  !soil carbon
write(fun%get_id('snn.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,30,'Average')  !soil nitrogen
write(fun%get_id('nep.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,5,'Add') - outputs(daily_out,7,'Add') !nep
write(fun%get_id('swc.dat'),'('' '',f8.1)',advance='NO') min(swc,9999.0)
write(fun%get_id('biomass.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,26,'Average') + &
 outputs(daily_out,27,'Average') + outputs(daily_out,28,'Average') ! biomass
write(fun%get_id('bioind.dat'),'('' '',i2)',advance='NO')   bioind
write(fun%get_id('covind.dat'),'('' '',i2)',advance='NO')   covind
write(fun%get_id('rof.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,2,'Add')
write(fun%get_id('fcn.dat'),'('' '',f8.2)',advance='NO') outputs(daily_out,32,'Add')
write(fun%get_id('nppstore.dat'),'('' '',f8.1)',advance='NO') avnppst
write(fun%get_id('stembio.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,27,'Average')  !stem biomass
write(fun%get_id('rootbio.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,28,'Average')  !root biomass
write(fun%get_id('leafper.dat'),'('' '',f8.1)',advance='NO') leafper
write(fun%get_id('stemper.dat'),'('' '',f8.1)',advance='NO') stemper
write(fun%get_id('rootper.dat'),'('' '',f8.1)',advance='NO') rootper
write(fun%get_id('sresp.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,7,'Add')       !soil respiration
write(fun%get_id('evt.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,3,'Add')
write(fun%get_id('gpp.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,6,'Add')       !gpp
write(fun%get_id('lch.dat'),'('' '',f8.3)',advance='NO') outputs(daily_out,31,'Add')      !soil leaching
write(fun%get_id('prc.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,10,'Add') !precipitation
write(fun%get_id('nbp.dat'),'('' '',f9.3)',advance='NO') outputs(daily_out,5,'Add') - outputs(daily_out,7,'Add') - &
 outputs(daily_out,32,'Add') - outputs(daily_out,31,'Add') - outputs(daily_out,33,'Add') !nbp avnpp-avsresp-firec-avlch-avyield
write(fun%get_id('trn.dat'),'('' '',f8.1)',advance='NO') outputs(daily_out,4,'Add')
write(fun%get_id('fab.dat'),'('' '',f8.5)',advance='NO') fprob
write(fun%get_id('tmp.dat'),'('' '',f8.2)',advance='NO') outputs(daily_out,9,'Average') !temperature
write(fun%get_id('hum.dat'),'('' '',f8.2)',advance='NO') outputs(daily_out,11,'Average') !humidity
!----------------------------------------------------------------------!
! Write var in crop output file                                        !
!----------------------------------------------------------------------!
    call crop_outputs(nft,3)
      endif

      if (iyear==nyears) then
write(fun%get_id('lai.dat'),*)
write(fun%get_id('npp.dat'),*)
write(fun%get_id('scn.dat'),*)
write(fun%get_id('snn.dat'),*)
write(fun%get_id('nep.dat'),*)
write(fun%get_id('swc.dat'),*)
write(fun%get_id('biomass.dat'),*)
write(fun%get_id('bioind.dat'),*)
write(fun%get_id('covind.dat'),*)
write(fun%get_id('rof.dat'),*)
write(fun%get_id('fcn.dat'),*)
write(fun%get_id('nppstore.dat'),*)
write(fun%get_id('stembio.dat'),*)
write(fun%get_id('rootbio.dat'),*)
write(fun%get_id('leafper.dat'),*)
write(fun%get_id('stemper.dat'),*)
write(fun%get_id('rootper.dat'),*)
write(fun%get_id('sresp.dat'),*)
write(fun%get_id('evt.dat'),*)
write(fun%get_id('gpp.dat'),*)
write(fun%get_id('lch.dat'),*)
write(fun%get_id('prc.dat'),*)
write(fun%get_id('nbp.dat'),*)
write(fun%get_id('trn.dat'),*)
write(fun%get_id('fab.dat'),*)
write(fun%get_id('tmp.dat'),*)
write(fun%get_id('hum.dat'),*)
      endif


!----------------------------------------------------------------------!
! Write optional cov bio bud sen.                                      !
!----------------------------------------------------------------------!
      iofn = iofngft
      if (iyear>=nyears-outyears2+1) then
        do ft=1,nft
          if (out_cov) then
write(fun%get_id('cov_'//trim(pft_tab(ft)%tag)//'.dat'),'('' '',f8.6)',advance='NO') covo(ft)
          endif
          if (out_bio) then
write(fun%get_id('bio_'//trim(pft_tab(ft)%tag)//'.dat'),'('' '',f8.6)',advance='NO') bioo(ft)
          endif
      if (out_bud) then
write(fun%get_id('bud_'//trim(pft_tab(ft)%tag)//'.dat'),'('' '',f8.6)',advance='NO') budo(ft)
      endif
      if (out_sen) then
write(fun%get_id('sen_'//trim(pft_tab(ft)%tag)//'.dat'),'('' '',f8.6)',advance='NO') seno(ft)
      endif
    enddo
  endif

!----------------------------------------------------------------------!
! Write selected tile and pft outputs.                                 !
!----------------------------------------------------------------------!
! Yearly tile
  if (iyear>=nyears-inp%output%tile_yearly+1) then
      do i=1,inp%output%tile_vars_yearly%n
          imap = inp%output%tile_vars_yearly%map(i)
          do mnth=1,12
            ans(mnth,1) = 0.0
          enddo
          do ft=1,nft
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                if (omav(imap)==1) then
                  oscale = 1.0/real(no_days(year,mnth,thty_dys))
                else
                  oscale = 1.0
                endif
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,1) = ans(mnth,1) + &
 daily_out(imap,co,mnth,day)*ssv(co)%cov*oscale
                enddo
              enddo
            enddo
          enddo
          iofn = iofn + 1
          sum1 = 0.0
          if (omav(imap)==1) then
            oscale = 1.0/12.0
          else
            oscale = 1.0
          endif
          do mnth=1,12
            sum1 = sum1 + ans(mnth,1)*oscale
          enddo
          write(fun%get_id('yearly_'//trim(otags(imap))//'.dat'),ofmt_yearly(imap),advance='no') sum1
      enddo
  endif


! Yearly pft
  if (iyear>=nyears-inp%output%pft_yearly+1) then
      do i=1,inp%output%pft_vars_yearly%n
          imap = inp%output%pft_vars_yearly%map(i)
          do ft=1,nft
            sumcov = 0.0
            do mnth=1,12
              ans(mnth,1) = 0.0
            enddo
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                if (omav(imap)==1) then
                  oscale = 1.0/real(no_days(year,mnth,thty_dys))
                else
                  oscale = 1.0
                endif
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,1) = ans(mnth,1) + daily_out(imap,co,mnth,day)*oscale*ssv(co)%cov
                enddo
              enddo
              sumcov = sumcov + ssv(co)%cov
            enddo
            iofn = iofn + 1
            if (sumcov > 0.0) then
              do mnth=1,12
                ans(mnth,1) = ans(mnth,1)/sumcov
              enddo
            endif
          sum1 = 0.0
          if (omav(imap)==1) then
            oscale = 1.0/12.0
          else
            oscale = 1.0
          endif
          do mnth=1,12
            sum1 = sum1 + ans(mnth,1)*oscale
          enddo
            write(fun%get_id('yearly_'//trim(otags(imap))//'_'//trim(pft_tab(ft)%tag)//'.dat'), &
  ofmt_yearly(imap)) sum1
          enddo
      enddo
  endif

! Monthly tile
  if (iyear>=nyears-inp%output%tile_monthly+1) then
      do i=1,inp%output%tile_vars_monthly%n
          imap = inp%output%tile_vars_monthly%map(i)
          do mnth=1,12
            ans(mnth,1) = 0.0
          enddo
          do ft=1,nft
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                if (omav(imap)==1) then
                  oscale = 1.0/real(no_days(year,mnth,thty_dys))
                else
                  oscale = 1.0
                endif
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,1) = ans(mnth,1) + &
 daily_out(imap,co,mnth,day)*ssv(co)%cov*oscale
                enddo
              enddo
            enddo
          enddo
          iofn = iofn + 1
          sum1 = 0.0
          if (omav(imap)==1) then
            oscale = 1.0/12.0
          else
            oscale = 1.0
          endif
          do mnth=1,12
            sum1 = sum1 + ans(mnth,1)*oscale
          enddo
write(fun%get_id('monthly_'//trim(otags(imap))//'.dat'),ofmt_monthly(imap)) year,(ans(mnth,1),mnth=1,12),sum1
      enddo
  endif

! Monthly pft
  if (iyear>=nyears-inp%output%pft_monthly+1) then
      do i=1,inp%output%pft_vars_monthly%n
          imap = inp%output%pft_vars_monthly%map(i)
          do ft=1,nft
            sumcov = 0.0
            do mnth=1,12
              ans(mnth,1) = 0.0
            enddo
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                if (omav(imap)==1) then
                  oscale = 1.0/real(no_days(year,mnth,thty_dys))
                else
                  oscale = 1.0
                endif
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,1) = ans(mnth,1) + daily_out(imap,co,mnth,day)*oscale*ssv(co)%cov
                enddo
              enddo
              sumcov = sumcov + ssv(co)%cov
            enddo
            iofn = iofn + 1
            if (sumcov > 0.0) then
              do mnth=1,12
                ans(mnth,1) = ans(mnth,1)/sumcov
              enddo
            endif
          sum1 = 0.0
          if (omav(imap)==1) then
            oscale = 1.0/12.0
          else
            oscale = 1.0
          endif
          do mnth=1,12
            sum1 = sum1 + ans(mnth,1)*oscale
          enddo
            write(fun%get_id('monthly_'//trim(otags(imap))//'_'//trim(pft_tab(ft)%tag)//'.dat'), &
  ofmt_monthly(imap)) year,(ans(mnth,1),mnth=1,12),sum1
          enddo
      enddo
  endif

! Daily tile
  if (iyear>=nyears-inp%output%tile_daily+1) then
      do i=1,inp%output%tile_vars_daily%n
          imap = inp%output%tile_vars_daily%map(i)
          write(fun%get_id('daily_'//trim(otags(imap))//'.dat'),'(i4)') year
          do mnth=1,12
            do day=1,no_days(year,mnth,thty_dys)
              ans(mnth,day) = 0.0
            enddo
          enddo
          do ft=1,nft
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,day) = ans(mnth,day) + daily_out(imap,co,mnth,day)*ssv(co)%cov
                enddo
              enddo
            enddo
          enddo
          do mnth=1,12
write(fun%get_id('daily_'//trim(otags(imap))//'.dat'),ofmt_daily(imap)) (ans(mnth,day),day=1,no_days(year,mnth,thty_dys))
          enddo
      enddo
  endif

! Daily pft
  if (iyear>=nyears-inp%output%pft_daily+1) then
      do i=1,inp%output%pft_vars_daily%n
          imap = inp%output%pft_vars_monthly%map(i)
          do ft=1,nft
            write(fun%get_id('daily_'//trim(otags(imap))//'_'//trim(pft_tab(ft)%tag)//'.dat'),'(i4)') year
            do mnth=1,12
              do day=1,no_days(year,mnth,thty_dys)
                ans(mnth,day) = 0.0
              enddo
            enddo
            sumcov = 0.0
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,day) = ans(mnth,day) + &
 daily_out(imap,co,mnth,day)*ssv(co)%cov
                enddo
              enddo
              sumcov = sumcov + ssv(co)%cov
            enddo
            if (sumcov > 0.0) then
              do mnth=1,12
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,day) = ans(mnth,day)/sumcov
                enddo
              enddo
            endif
            do mnth=1,12
              write(fun%get_id('daily_'//trim(otags(imap))//'_'//trim(pft_tab(ft)%tag)//'.dat'), &
  ofmt_daily(imap)) (ans(mnth,day),day=1,no_days(year,mnth,thty_dys))
            enddo
          enddo
      enddo
  endif

!----------------------------------------------------------------------!
! Write snapshots
!----------------------------------------------------------------------!
  snp_year = 1
  if (snp_no>0) then
    if (snp_year<=snp_no) then
      if (year==snpshts(snp_year)) then
        fno = 100 + (snp_year-1)*4
        write(fno+1,'(F7.3,F9.3,4500F9.1)')  &
 ssp%lat,ssp%lon,(((ssv(ft)%bio(j),j=1,2),ft=1,pft(k)%mort),k=1,nft)
        write(fno+2,'(F7.3,F9.3,1500F12.9)') &
 ssp%lat,ssp%lon,((ssv(ft)%cov,ft=1,pft(j)%mort),j=1,nft)
        write(fno+3,'(F7.3,F9.3,1500F12.7)') &
 ssp%lat,ssp%lon,((ssv(ft)%ppm,ft=1,pft(j)%mort),j=1,nft)
        write(fno+4,'(F7.3,F9.3,1500F8.3)')  &
 ssp%lat,ssp%lon,((ssv(ft)%hgt,ft=1,pft(j)%mort),j=1,nft)
        snp_year = snp_year + 1
      endif
    endif
  endif

    enddo
!======================================================================!
!                         END OF YEARLY LOOP                           !
!======================================================================!

!----------------------------------------------------------------------!
!                      End record for default output files             !
!----------------------------------------------------------------------!
!  do i=1,fun%n
!    if (fun%opened(i)) then
!      write(fun%get_id_n(i),*)
!      call fun%close_n(i)
!    endif
!  enddo


!----------------------------------------------------------------------!
!                      End record for default output files             !
!----------------------------------------------------------------------!

  iofn = iofngft
  if (outyears2>0) then
    do ft=1,nft
      if (out_cov) then
write(fun%get_id('cov_'//trim(pft_tab(ft)%tag)//'.dat'),*)
      endif
      if (out_bio) then
write(fun%get_id('bio_'//trim(pft_tab(ft)%tag)//'.dat'),*)
      endif
      if (out_bud) then
write(fun%get_id('bud_'//trim(pft_tab(ft)%tag)//'.dat'),*)
      endif
      if (out_sen) then
write(fun%get_id('sen_'//trim(pft_tab(ft)%tag)//'.dat'),*)
      endif
    enddo
  endif

!----------------------------------------------------------------------!
! Open file to output system state.                                    !
!----------------------------------------------------------------------!
if (inp%run%output_state) then
  nco = ssp%cohorts
  nft = ssp%nft

  fid = fid_output_state
  write(fid) ssp%cohorts
  write(fid) ssp%nft
  write(fid) ssv(1:nco)
  write(fid) pft(1:nco)
  write(fid) ssp%day,ssp%mnth,ssp%year,ssp%iyear,ssp%lat,ssp%lon &
 ,ssp%latres,ssp%lonres,ssp%nft,ssp%lai,ssp%cohort,ssp%jday &
 ,ssp%soil_depth,ssp%sand,ssp%silt,ssp%clay,ssp%bulk,ssp%orgc,ssp%wilt &
 ,ssp%field,ssp%sat,ssp%tmem,ssp%new_soil_h2o,ssp%new_snow &
 ,ssp%new_l_snow,ssp%new_c,ssp%new_n,ssp%new_minn,ssp%new_slc &
 ,ssp%new_rlc,ssp%new_sln,ssp%new_rln,ssp%new_cov &
 ,ssp%xnew_soil_h2o(1:nft,:),ssp%xnew_snow(1:nft),ssp%xnew_l_snow(1:nft) &
 ,ssp%xnew_c(1:nft,:),ssp%xnew_n(1:nft,:),ssp%xnew_minn(1:nft,:) &
 ,ssp%xnew_slc(1:nft),ssp%xnew_rlc(1:nft),ssp%xnew_sln(1:nft) &
 ,ssp%xnew_rln(1:nft),ssp%xnew_cov(1:nft),ssp%mnthtmp,ssp%mnthprc &
 ,ssp%mnthhum,ssp%emnthtmp,ssp%emnthprc,ssp%emnthhum,ssp%iseas &
 ,ssp%cohorts,ssp%co2ftmap(1:nco,1:nco),ssp%ftcov(1:nft)
    endif

  else
  write(fun%get_id('diag.dat'),*) &
 '                  clm stt ssc blk wfs dep lus'
  write(fun%get_id('diag.dat'),'(f7.3,f9.3,1x,20L4)') &
 ssp%lat,ssp%lon,l_clim,l_stats,l_soil(1),l_soil(3),l_soil(5),l_soil(8),l_lu

  endif
!======================================================================!
!                         END OF DRIVING DATA EXISTS IF                !
!======================================================================!

!----------------------------------------------------------------------!
! Skip line in crop output files.                                      !
!----------------------------------------------------------------------!
call crop_outputs(nft,4)

enddo
!======================================================================!
!                         END OF SITE LOOP                             !
!======================================================================!


!----------------------------------------------------------------------!
! Open file to record version number, command line, input file and     !
! parameter file.                                                      !
!----------------------------------------------------------------------!
call fun%open(trim(inp%dirs%output)//'/simulation.dat',fid)

call getarg(0,buff1)
write(fid,'(A)') stver(1:28)
write(fid,*)

write(fid,'(''************************************************'')')
write(fid,'(''* Command line arguments                       !'')')
write(fid,'(''************************************************'')')

do i=0,iargc()
  call getarg(i,buff1)
  write(fid,'(A,'' '')',advance='no') trim(buff1)
enddo
write(fid,*)

write(fid,'(''************************************************'')')
write(fid,'(''* Parameter file (param.dat)                   !'')')
write(fid,'(''************************************************'')')

call fun%open_old('inc/param.dat',i)
do
  read(i,'(A)',iostat=kode) st1
  if (kode/=0) exit
  write(fid,'(A)') trim(st1)
enddo
call fun%close(i)
write(fid,*)

write(fid,'(''************************************************'')')
write(fid,'(''* Input file                                   !'')')
write(fid,'(''************************************************'')')

call getarg(1,buff1)
call fun%open_old(buff1,i)
do
  read(i,'(A)',iostat=kode) st1
  if (kode/=0) exit
  write(fid,'(A)') trim(st1)
enddo
call fun%close(i)

write(fid,'(''************************************************'')')
!----------------------------------------------------------------------!

call fun%close(fid)
call fun%close(fun%get_id('site_info.dat'))
call fun%close(fun%get_id('diag.dat'))

iofn = iofngft
if (outyears2>0) then
  do ft=1,nft
    if (out_cov) then
call fun%close(fun%get_id('cov_'//trim(pft_tab(ft)%tag)//'.dat'))
    endif
    if (out_bio) then
call fun%close(fun%get_id('bio_'//trim(pft_tab(ft)%tag)//'.dat'))
    endif
    if (out_bud) then
call fun%close(fun%get_id('bud_'//trim(pft_tab(ft)%tag)//'.dat'))
    endif
    if (out_sen) then
call fun%close(fun%get_id('sen_'//trim(pft_tab(ft)%tag)//'.dat'))
    endif
  enddo
endif

do i=70,89
  close(i)
enddo
do i=91,93
  close(i)
enddo
do i=98,99
  close(i)
enddo

if (snp_no>0) then
  do i=1,snp_no
    fno = 100+(i-1)*4
    close(fno+1)
    close(fno+1)
    close(fno+3)
    close(fno+4)
  enddo
endif

!----------------------------------------------------------------------!
! Close crop output files.                                             !
!----------------------------------------------------------------------!
call crop_outputs(nft,1)

!call fun%print()

end program sdgvm



