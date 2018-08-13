module read_input

use real_precision
use dims
use pft_parameters
use site_parameters
use sdgvm1
use input_methods
use func
use misc_parameters
use tuning_parameters
use open_files
use file_class
use file_object

implicit none

contains

!**********************************************************************!
!                                                                      !
!                    read_input_file :: read_input                     !
!                    -----------------------------                     !
!                                                                      !
! SUBROUTINE read_input_file(buff1,stco2,xlatf,xlon0,xlatres,xlonres,  !
! speedc,xspeedc,xseed1,spinl,                             !
! crand,yr0p,yrfp,outyears,nyears,yr0,yrf,yearind,idum,yearv,nomdos,   !
! otags,omav,ofmt,outyears1,outyears2,oymd,otagsn,otagsnft,snpshts,    !
! snp_no,out_cov,out_bio,out_bud,out_sen,lutab,grassrc,barerc,fireres, !
! luse,l_b_and_c,soil_chr,topsl,defaulttopsl,sites,latdel,londel,      !
! lat_lon,day_mnth,thty_dys,xlatresn,xlonresn,ilanduse,nft,xyear0,     !
! xyearf,lmor_sc,oymdft,iofnft,sit_grd,du,narg)                        !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine process_input_file(buff1,xlatf,xlon0,xlatres,xlonres,&
 speedc,xspeedc,crand,&
 outyears,nyears,yr0,yrf,yearind,idum,yearv,nomdos,otags,omav, &
 ofmt_daily,ofmt_monthly,ofmt_yearly, &
 outyears1,outyears2,oymd,otagsn,otagsnft,snpshts,snp_no,out_cov,&
 out_bio,out_bud,out_sen,lutab,grassrc,barerc,fireres,luse,l_b_and_c,&
 soil_chr,topsl,defaulttopsl,sites,latdel,londel,lat_lon,day_mnth,&
 thty_dys,xlatresn,xlonresn,ilanduse,nft,xyear0,xyearf,lmor_sc,&
 oymdft,iofnft,sit_grd,du,narg,fire_ant,harvest_ant,met_seq,par_loops)
!**********************************************************************!
integer :: fireres,sites,per,site,ibox,jbox,l,recl1, &
 site0,sitef,luse(max_years),ilanduse,persum,clim_type
real(dp) :: latdel,londel,lat_lon(max_sites,2),iadj,jadj, &
 grassrc,barerc,topsl,defaulttopsl,soil_chr(10),lutab(255,100),&
 lmor_sc(3600,max_cohorts),lutab2(255,100)
character(len=str_len) :: st1,st2,st3,st4,st5
character(len=str_len), dimension(max_outputs) :: otags,ofmt_daily, &
 ofmt_monthly,ofmt_yearly,fttags
character(len=80) :: buff1
integer, dimension(max_years) :: snpshts
integer :: i,du,xyear0,xyearf,xlatresn,xlonresn,kode,ii, &
 xseed1,spinl,yr0s,cycle,yr0p,yrfp,outyears,yr0,yrf,idum, &
 nyears,yearv(max_years),yearind(max_years),nomdos,omav(max_outputs),&
 oymd,otagsn(max_outputs),otagsnft(max_outputs),outyears2,snp_no, &
 sit_grd,day_mnth,thty_dys,narg,j,iofnft,iofn,ft,nft,k,oymdft,outyears1
real(dp) :: xlatf,xlon0,lon0,lonf,lat0,latf,xlatres,xlonres
logical :: speedc,xspeedc,crand
logical :: out_cov,out_bio,out_bud,out_sen,l_b_and_c
logical :: fire_ant(max_years),harvest_ant(max_years),met_seq

integer :: subd_par,read_clump,calc_zen,no_slw_lim,cstype,ncalc_type
integer :: ttype,vcmax_type,coilp_map,s070607,gs_func,soilcn_map
integer :: phen_cor,hw_j,read_par,soil_map,soilp_map,daily_co2,par_loops
logical :: goudriaan_old
real(dp) :: p_pet,p_et
integer :: yr0ms,yrfms,yr0m,yrfm,met_seqv(max_years),met_yearv(max_years)
integer :: fid

logical :: logic

integer :: nlat,nlon
!----------------------------------------------------------------------!

      met_seq = .FALSE.
      if(ii>1) THEN
        st3 = 'seq'
        CALL STRIPBS(st1,st2)
        if (stcmp(st2,st3)==1) then
          met_seq = .TRUE.
        ELSE
          write(*,*) &
 'If it exists, second field of second line in input.dat must read "seq"' 
           stop 
        endIF
      endIF

      open(newunit=fid,FILE=trim(inp%dirs%climate)//'/readme.dat', &
 STATUS='OLD',iostat=kode)
      if (kode/=0) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'Climate data file does not exist:'
        write(*,'('' "'',A,''/readme.dat"'')') trim(inp%dirs%climate)
        stop
      endIF
      READ(fid,'(A)') st1

!----------------------------------------------------------------------!
! Determine whether climate is daily or monthly, from first line of    !
! readme file.                                                         !
!----------------------------------------------------------------------!
st2 = 'DAILY'
st3 = 'MONTHLY'
st4 = 'SITED'
st5 = 'SITEM'
sit_grd = 0
if (stcmp(st1,st2)==1) then
  day_mnth = 1
  thty_dys = 1
  read(fid,*)
  read(fid,*) xlatf,xlon0
  read(fid,*)
  read(fid,*) xlatres,xlonres
  read(fid,*)
  read(fid,*) xlatresn,xlonresn
  read(fid,*)
  read(fid,*) xyear0,xyearf
elseif (stcmp(st1,st3)==1) then
  clim_type = 2
  day_mnth = 0
  thty_dys = 1
  read(fid,*)
  read(fid,*) xlatf,xlon0
  read(fid,*)
  read(fid,*) xlatres,xlonres
  read(fid,*)
  read(fid,*) xlatresn,xlonresn
  read(fid,*)
  read(fid,*) xyear0,xyearf
elseif (stcmp(st1,st4)==1) then
  clim_type = 3
  day_mnth = 1
  thty_dys = 0
  sit_grd = 1
  read(fid,*)
  read(fid,*) xlatf,xlon0
  read(fid,*)
  read(fid,*) xyear0,xyearf
elseif (stcmp(st1,st5)==1) then
  clim_type = 4
  day_mnth = 0
  thty_dys = 1
  sit_grd = 1
  read(fid,*)
  read(fid,*) xlatf,xlon0
  read(fid,*)
  read(fid,*) xyear0,xyearf
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) &
 'First line of climat readme.dat file must read DAILY MONTHLY or SITE.'
  stop
endif
close(fid)

!----------------------------------------------------------------------!
! read input switches                                                  !
!----------------------------------------------------------------------!

      !line 9 in <input.dat> determines SDVGM version
      !0 is current version 1 is old 070607 version
!      READ(fid_input,*)      

      !read next line of input switches. ln 10 in input.dat
!      READ(fid_input,*)
!!      ii = n_fields(st1) 
!      IF (ii==6) THEN 
!        CALL STRIPBN(st1,i)
!        !read in daily co2 
!        IF (i>-1)  daily_co2  = i 
!        CALL STRIPBN(st1,i) 
!        !read PAR data 
!        IF (i>-1)  read_par   = i 
!        CALL STRIPBN(st1,i)
!        !use sub-daily PAR scaling
!        IF (i>-1)  subd_par   = i 
!        CALL STRIPBN(st1,i) 
!        !canopy clumping index: 0-do not use; 1-from pft values; 2-from a map
!        IF (i>-1)  read_clump = i 
!        CALL STRIPBN(st1,i) 
!        !calculate solar zenith angle and use in canopy light interception 
!        IF (i>-1)  calc_zen   = i 
!        CALL STRIPBN(st1,i) 
!        !no soil water limitation
!        IF (i>-1)  no_slw_lim = i 
!      ELSE !
!        WRITE(*,'('' PROGRAM TERMINATED'')') !
!        WRITE(*,*) 'Line 10 must contain 6 fields' !
!        WRITE(*,'('' "'',A,''"'')') st1(1:30) !
!        STOP !       
!      ENDIF !

      if(inp%run%subdaily) then
        write(*,*) ''
        write(*,'('' SDGVM running with '',i3, &
 '' sub-daily photosynthesis time-points'')') par_loops
      else
        write(*,*) ''
        write(*,*) 'SDGVM running with 1 sub-daily photosynthesis time-points'
      endif

      !read next line of input switches. ln 11 in input.dat
!      READ(fid_input,*)
!      ii = n_fields(st1) 
!      IF (ii==5) THEN 
!        CALL STRIPBN(st1,i)
!        !select canopy nitrogen calculation method
!        cstype = 0
!        IF (i>=10) THEN
!          cstype = int(real(i)/10.0)
!          i      = i - 10*cstype
!        ENDIF
!        IF (i>-1)  ncalc_type = i 
!        CALL STRIPBN(st1,i) 
!        !select vcmax parameterisation
!        ttype = 0
!        IF (i>=10) THEN
!          ttype = int(real(i)/10.0)
!          i     = i - 10*ttype
!        ENDIF
!        IF (i>-1)  vcmax_type = i 
!        CALL STRIPBN(st1,i) 
!        !use soil P from a map to calculate Vcmax (Vcmax switch must also be 1 for P to be used in Vcmax calc) 
!        !- even if this variable is 0, it is expected to be read from the soils database 
!        IF (i>-1)  soilp_map = i 
!        CALL STRIPBN(st1,i) 
!        !switch to run with 070607 routines that have been changed in the new version but don't have individual switches (rd, swlim, et multiplier etc)
!        IF (i>-1)  s070607   = i 
!        CALL STRIPBN(st1,i) 
!        !switch stomatal conductance function
!        IF (i>-1)  gs_func = i 
!        !CALL STRIPBN(st1,i) 
!        !switch electron transport function
!        !IF (i>-1)  hw_j = i 
!      ELSE !
!        WRITE(*,'('' PROGRAM TERMINATED'')') !
!        WRITE(*,*) 'Line 11 must contain 5 fields' !
!        WRITE(*,'('' "'',A,''"'')') st1(1:30) !
!        STOP !
!      ENDIF !

      
      !the below switches are hard coded as they are unlikely to need changing
      ! - they are only changed if the old (070607) version of the model is used  

      !use soil C:N ratio from a map 
      !- even if this variable is 0, it is expected to be read from the soils database
      ! - 0 is the default for this switch, use this switch to implement a routine that uses soil C:N read from the soils database 
      soilcn_map = 0  

      !grass lai can only take 50% of stored C and leaf growth subject to growth respiration
      phen_cor   = 1 

      !switch electron transport function, 0 - Harley 1992, 1 - Farquhar
      !& Wong 1984
      hw_j       = 0

      !this is commented out for ease of adding new swithces to the input.dat for model development
      !read next line of input switches. ln 12 in input.dat
!      READ(98,'(1000a)') st1 
!      ii = n_fields(st1) 
!      IF (ii==4) THEN 
!      CALL STRIPBN(st1,i) 
!      !master switch for these switches
!      IF (i>-1)  mswitch = i 
!      CALL STRIPBN(st1,i) 
!      !switch 1 
!      IF (i>-1)  switch1 = i 
!      !switch 2
!      IF (i>-1)  switch2 = i 
!      !switch 3
!      IF (i>-1)  switch3 = i 
!      ELSE !
!        WRITE(*,'('' PROGRAM TERMINATED'')') !
!        WRITE(*,*) 'Line 12 must contain 4 fields' !
!        WRITE(*,'('' "'',A,''"'')') st1(1:30) !
!        STOP !   
!      ENDIF !

      if(gs_func==3) goudriaan_old = .TRUE.

      !set default configurations for standard versions
      if (inp%run%subdaily) then
        daily_co2  = 0 
        read_par   = 1
        subd_par   = 1
        read_clump = 0
        calc_zen   = 1 
        no_slw_lim = 0

        cstype     = 0
        ncalc_type = 1
        ttype      = 0
        vcmax_type = 1
        soilp_map  = 0 
        s070607    = 0
        gs_func    = 0

        soilcn_map = 0
        phen_cor   = 1
        hw_j       = 0
      ELSE
        daily_co2  = 0 
        read_par   = 0
        subd_par   = 0 
        read_clump = 0
        calc_zen   = 0
        no_slw_lim = 0

        cstype     = 0  
        ncalc_type = 0
        ttype      = 0
        vcmax_type = 0
        soilp_map  = 0
        s070607    = 1
        gs_func    = 0

        soilcn_map = 0
        phen_cor   = 0 
        hw_j       = 0
      endIF

      if(goudriaan_old) hw_j = 3

      !parameters for 070607 version
      if(s070607==1) then
        p_et  = 0.7
        p_pet = 0.7
      endif

!read(fid_input,*)

!read(fid_input,*)

st2 = 'ARGUMENT'
if (stcmp(inp%dirs%input,st2)==1) then
  call GETARG(narg,buff1)
  do i=1,80
    inp%dirs%input(i:i) = buff1(i:i)
  enddo
endif

!read(fid_input,*)

st2 = 'ARGUMENT'
if (stcmp(inp%dirs%output,st2)==1) then
  call GETARG(narg,buff1)
  do i=1,80
    inp%dirs%output(i:i) = buff1(i:i)
  enddo
endif
!read(fid_input,*)

!----------------------------------------------------------------------!
! Check if output directory exists.                                    !
!----------------------------------------------------------------------!
inquire(file=trim(inp%dirs%output)//'/.',exist=logic)
if (.not.logic) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'SDGVM output directory does not exist.'
  write(*,'('' "'',A,''"'')') trim(inp%dirs%output)
  stop
endif
!----------------------------------------------------------------------!

!read(fid_input,*)
!ii = n_fields(st1)
!if ((ii==3).or.(ii==4)) then
!  call STRIPBN(st1,i)
!  call STRIPBN(st1,j)
!  call STRIPBN(st1,k)
!  speedc = .true.
!  if (k==0)  speedc = .false.
!else
!  write(*,'('' PROGRAM TERMINATED'')')
!  write(*,*) 'Line 16 must contain 3 or 4 arguments'
!  write(*,'('' "'',A,''"'')') st1(1:30)
!  stop
!endif
xspeedc = speedc

!read(fid_input,*)
!call STRIPB(st1)
!i = n_fields(st1)
!if (i==7) then
!  call STRIPBN(st1,spinl)
!  call STRIPBN(st1,yr0s)
!  call STRIPBN(st1,cycle)
!  call STRIPBN(st1,j)
!  crand = .true.
!  if (j==0)  crand = .false.
!  call STRIPBN(st1,yr0p)
!  call STRIPBN(st1,yrfp)
!  call STRIPBN(st1,outyears)
!elseif (i==6) then
!  call STRIPBN(st1,spinl)
!  call STRIPBN(st1,yr0s)
!  call STRIPBN(st1,cycle)
!  call STRIPBN(st1,j)
!  crand = .true.
!  if (j==0)  crand = .false.
!  call STRIPBN(st1,yr0p)
!  call STRIPBN(st1,yrfp)
!  outyears = yrfp - yr0p + 1
!elseif (i==5) then
!  call STRIPBN(st1,spinl)
!  call STRIPBN(st1,yr0s)
!  call STRIPBN(st1,cycle)
!  call STRIPBN(st1,j)
!  crand = .true.
!  if (j==0)  crand = .false.
!  yr0p = yr0s+1
!  yrfp = yr0s
!  call STRIPBN(st1,outyears)
!elseif (i==4) then
!  call STRIPBN(st1,spinl)
!  call STRIPBN(st1,yr0s)
!  call STRIPBN(st1,cycle)
!  call STRIPBN(st1,j)
!  crand = .true.
!  if (j==0)  crand = .false.
!  yr0p = yr0s+1
!  yrfp = yr0s
!  outyears = inp%run%spinup_cycle_length + 1
!elseif (i==3) then
!  spinl = 0
!  inp%run%spinup_cycle_length = max_years
!  j = 0
!  crand = .false.
!  call STRIPBN(st1,yr0p)
!  call STRIPBN(st1,yrfp)
!  call STRIPBN(st1,outyears)
!  yr0s = yr0p
!elseif (i==2) then
!  spinl = 0
!  inp%run%spinup_cycle_length = max_years
!  j = 0
!  crand = .false.
!  call STRIPBN(st1,yr0p)
!  call STRIPBN(st1,yrfp)
!  yr0s = yr0p
!  outyears = yrfp - yr0p + 1
!else
!  write(*,'('' PROGRAM TERMINATED'')')
!  write(*,*) 'Line 17 must contain 2-7 fields'
!  stop
!endif
nyears = inp%run%spinup_length + inp%run%yearf - inp%run%year0 + 1
outyears = min(nyears,inp%output%nyears)

!yrfp = 0
!yr0p = 0

if (nyears>max_years) then
  write(*,'('' PROGRAM TERMINATED'')')
 write(*,'('' Trying to simulate '',i4,&
 &'' years, maximum allowable is '',i4,''.'')') nyears,max_years
 write(*,*) &
 'Either reduce the length of the simulation, or increase "max_years"'
  write(*,*) &
 'this is set in array_param.txt, you must re-comile after altering'
  write(*,*) 'this file.'
  stop
endif
if ((j/=0).and.(j/=1)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,'('' Fourth field of line 14 in the input file must be &
 &either 0 or 1.'')')
  write(*,'('' Currently set to '',i5,''.'')') j
  stop
endif


!----------------------------------------------------------------------!
! Read in a met sequence if specified                                  !
!----------------------------------------------------------------------!
      yr0ms = yr0p
      yrfms = yrfp 
      if(met_seq) THEN
        open(99,FILE='met_seq.dat',iostat=kode)
        if (kode/=0) then
          write(*,'('' PROGRAM TERMINATED'')')
          write(*,*) &
 'non-sequential met sequence selected but file does not exist.'
          write(*,'('' "'',A,''"'')') 'met_seq.dat'
          stop
        endIF
        do i=1,(yrfp-yr0p+1)
          READ(99,*) met_seqv(i)
          if(i==1) yr0ms = met_seqv(i)
          if(i>1) yr0ms = min(yr0ms,met_seqv(i))
          if(i==1) yrfms = met_seqv(i)
          if(i>1) yrfms = max(yrfms,met_seqv(i))
        endDO
        close(99)
      endIF

!----------------------------------------------------------------------!
! Set 'yr0' and 'yrf' which are the years of actual climate required   !
! for the run. And check that the climate exists in the climate        !
! database                                                             !
!----------------------------------------------------------------------!
      yr0 = min(inp%run%spinup_year0,inp%run%year0)
      yrf = max(inp%run%spinup_year0+min(inp%run%spinup_length,inp%run%spinup_cycle_length)-1,inp%run%yearf)
      yr0m = yr0
      yrfm = yrf
if(met_seq) yr0m = min(yr0s,yr0ms)
if(met_seq) yrfm = max(yr0s+min(inp%run%spinup_length,inp%run%spinup_cycle_length)-1,yrfms)
if ((yr0m<xyear0).or.(yrfm>xyearf)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,'('' Trying to use '',i4,''-'',i4,'' climate.'')') yr0m,yrfm 
  write(*,'('' Climate database runs from '',i4,''-'',i4,''.'')') xyear0,xyearf
  stop
endIF

!----------------------------------------------------------------------!
! Set up year climate sequence.                                        !
!----------------------------------------------------------------------!
do i=1,max_years
  yearind(i) = i
enddo
idum = 1
do i=1,nyears
  if (i<=inp%run%spinup_length) then
    if (crand) then
      if ((mod(i-1,inp%run%spinup_cycle_length)==0).and.(i>2)) &
      call RANDOMV(yearind,1,inp%run%spinup_cycle_length,idum)
      yearv(i) = yearind(mod(i-1,inp%run%spinup_cycle_length)+1) + inp%run%spinup_year0 - 1
    else
      yearv(i) = mod(i-1,inp%run%spinup_cycle_length) + inp%run%spinup_year0
    endif
  else
    yearv(i) = i - inp%run%spinup_length + inp%run%year0 - 1
    if (met_seq) then
      met_yearv(i) = met_seqv(i-inp%run%spinup_length)
    ELSE
      met_yearv(i) = yearv(i)
    endIF
  endif
enddo

!----------------------------------------------------------------------!
! Open monthly/daily PIXEL output files if required.                   !
!----------------------------------------------------------------------!
call OUTPUT_OPTIONS(nomdos,otags,omav,ofmt_daily,ofmt_monthly,ofmt_yearly)

!read(fid_input,*)
!call STRIPB(st1)
!st2 = 'PIXEL'
!if (stcmp(st1,st2)==0) then
!  write(*,'('' PROGRAM TERMINATED'')')
!  write(*,*) &
! 'The first field of line 15 of the input file must read ''PIXEL''.'
!  write(*,*) st1(1:30)
!  write(*,*) 'Output variable options:'
!  write(*,'(1x,20a4)') (otags(i),i=1,15)
!  write(*,'(1x,20a4)') (otags(i),i=16,nomdos)
!  stop
!endif


call SET_PIXEL_OUT(st1,st2,outyears1,nomdos,otagsn,otags,oymd)
outyears1 = min(outyears1,nyears)

st1 = inp%dirs%output


!----------------------------------------------------------------------!
! Determine whether daily or monthly subpixel outputs are required.    !
!----------------------------------------------------------------------!

!read(fid_input,*)

call SET_SUBPIXEL_OUT(st1,st2,outyears2,nomdos,otagsnft,otags,oymdft,&
 out_cov,out_bio,out_bud,out_sen)
outyears2 = min(outyears2,nyears)

!----------------------------------------------------------------------!
! Read in snapshot years.                                              !
!----------------------------------------------------------------------!
!read(fid_input,*)
!call STRIPB(st1)
!st2 = 'SNAPSHOTS'
!if (stcmp(st1,st2)==0) then
!  write(*,'('' PROGRAM TERMINATED'')')
!  write(*,*) 'First field of line 17 must read SNAPSHOTS'
!  stop
!endif
!call STRIPBS(st1,st2)
!call STRIPB(st1)
!call ST2ARR(st1,snpshts,100,snp_no)
!read (fid_input,*)

!----------------------------------------------------------------------!
! Read in compulsory functional types.                                 !
!----------------------------------------------------------------------!
pft_tab(1)%tag = 'BARE'
pft_tab(2)%tag = 'CITY'
pft_tab(1)%itag = 1
pft_tab(2)%itag = 2
!read(fid_input,*)
!call STRIPB(st1)
!st2 = 'BARE'
!if ((stcmp(st2,st1)==0).or.(n_fields(st1)/=2)) then
!  write(*,'('' PROGRAM TERMINATED'')')
!  write(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1).'
!  stop
!endif
!call STRIPBS(st1,st2)
!read(st1,*) pft_tab(1)%mix
!read(fid_input,*)
!call STRIPB(st1)
!st2 = 'CITY'
!if ((stcmp(st2,st1)==0).or.(n_fields(st1)/=2)) then
!  write(*,'('' PROGRAM TERMINATED'')')
!  write(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1).'
!  stop
!endif
!call STRIPBS(st1,st2)
!read(st1,*) pft_tab(2)%mix
!read (fid_input,*)

! Initialise redundant parameterisation
ft = 1
pft_tab(ft)%c3c4 = .true.
pft_tab(ft)%phen = 0
pft_tab(ft)%crop = 0.0
pft_tab(ft)%d2h = 0
pft_tab(ft)%mort = 1
pft_tab(ft)%wden = 0.0
pft_tab(ft)%xyl = 0.0
pft_tab(ft)%pdif = 0.0
pft_tab(ft)%sla = 0.0
pft_tab(ft)%lls = 0
pft_tab(ft)%sls = 0
pft_tab(ft)%rls = 0
pft_tab(ft)%lmor = 0.0
pft_tab(ft)%lrat = 0.0
pft_tab(ft)%bbmem = 0
pft_tab(ft)%bb0 = 0.0
pft_tab(ft)%bbmax = 0.0
pft_tab(ft)%bblim = 0.0
pft_tab(ft)%senm = 0
pft_tab(ft)%sens = 0
pft_tab(ft)%senlim = 0.0
pft_tab(ft)%stemx = 0.0
pft_tab(ft)%gr0 = 0.0
pft_tab(ft)%grf = 0.0
pft_tab(ft)%ppm0 = 0.0
pft_tab(ft)%sowthresh(1)=0.0
pft_tab(ft)%sowthresh(2)=0.0
pft_tab(ft)%lethal(1)=0.0
pft_tab(ft)%lethal(2)=0.0
pft_tab(ft)%cardinal(1)=0.0
pft_tab(ft)%cardinal(2)=0.0
pft_tab(ft)%cardinal(3)=0.0
pft_tab(ft)%cardinal(4)=0.0
pft_tab(ft)%cardinal(5)=0.0
pft_tab(ft)%cardinal(6)=0.0
pft_tab(ft)%cardinal(7)=0.0
pft_tab(ft)%cardinal(8)=0.0
pft_tab(ft)%cardinal(9)=0.0
pft_tab(ft)%croptype(1)=0.0
pft_tab(ft)%croptype(2)=0.0
pft_tab(ft)%photoperiod(1)=0.0
pft_tab(ft)%photoperiod(2)=0.0
pft_tab(ft)%photoperiod(3)=0.0
pft_tab(ft)%photoperiod(4)=0.0
pft_tab(ft)%photoperiod(5)=0.0
pft_tab(ft)%photoperiod(6)=0.0
pft_tab(ft)%croprange(1)=0.0
pft_tab(ft)%croprange(2)=0.0
pft_tab(ft)%croprange(3)=0.0
pft_tab(ft)%croprange(4)=0.0
pft_tab(ft)%cropphen(1)=0.0
pft_tab(ft)%cropphen(2)=0.0
pft_tab(ft)%cropphen(3)=0.0
pft_tab(ft)%cropphen(4)=0.0
pft_tab(ft)%cropphen(5)=0.0
pft_tab(ft)%cropphen(6)=0.0
pft_tab(ft)%irrig(1)=0.0
pft_tab(ft)%irrig(2)=0.0
pft_tab(ft)%irrig(3)=0.0
pft_tab(ft)%sowday(1)=0
pft_tab(ft)%sowday(2)=0
pft_tab(ft)%sowday(3)=0
pft_tab(ft)%cropgdd(1,1)=0
pft_tab(ft)%cropgdd(1,2)=0
pft_tab(ft)%cropgdd(1,3)=0
pft_tab(ft)%cropgdd(2,1)=0
pft_tab(ft)%cropgdd(2,2)=0
pft_tab(ft)%cropgdd(2,3)=0
pft_tab(ft)%fert(1)=0.
pft_tab(ft)%fert(2)=0.
pft_tab(ft)%fert(3)=0.
pft_tab(ft)%fert(4)=0.
pft_tab(ft)%fert(5)=0.
pft_tab(ft)%fert(6)=0.
pft_tab(ft)%optlai=0.
pft_tab(ft)%harvindx=0.
pft_tab(ft)%limdharv=0

ft = 2
pft_tab(ft)%c3c4 = .true.
pft_tab(ft)%phen = 0
pft_tab(ft)%crop = 0.0
pft_tab(ft)%d2h = 0
pft_tab(ft)%mort = 1
pft_tab(ft)%wden = 0.0
pft_tab(ft)%xyl = 0.0
pft_tab(ft)%pdif = 0.0
pft_tab(ft)%sla = 0.0
pft_tab(ft)%lls = 0
pft_tab(ft)%sls = 0
pft_tab(ft)%rls = 0
pft_tab(ft)%lmor = 0.0
pft_tab(ft)%lrat = 0.0
pft_tab(ft)%bbmem = 0
pft_tab(ft)%bb0 = 0.0
pft_tab(ft)%bbmax = 0.0
pft_tab(ft)%bblim = 0.0
pft_tab(ft)%senm = 0
pft_tab(ft)%sens = 0
pft_tab(ft)%senlim = 0.0
pft_tab(ft)%stemx = 0.0
pft_tab(ft)%gr0 = 0.0
pft_tab(ft)%grf = 0.0
pft_tab(ft)%ppm0 = 0.0
pft_tab(ft)%sowthresh(1)=0.0
pft_tab(ft)%sowthresh(2)=0.0
pft_tab(ft)%lethal(1)=0.0
pft_tab(ft)%lethal(2)=0.0
pft_tab(ft)%cardinal(1)=0.0
pft_tab(ft)%cardinal(2)=0.0
pft_tab(ft)%cardinal(3)=0.0
pft_tab(ft)%cardinal(4)=0.0
pft_tab(ft)%cardinal(5)=0.0
pft_tab(ft)%cardinal(6)=0.0
pft_tab(ft)%cardinal(7)=0.0
pft_tab(ft)%cardinal(8)=0.0
pft_tab(ft)%cardinal(9)=0.0
pft_tab(ft)%croptype(1)=0.0
pft_tab(ft)%croptype(2)=0.0
pft_tab(ft)%photoperiod(1)=0.0
pft_tab(ft)%photoperiod(2)=0.0
pft_tab(ft)%photoperiod(3)=0.0
pft_tab(ft)%photoperiod(4)=0.0
pft_tab(ft)%photoperiod(5)=0.0
pft_tab(ft)%photoperiod(6)=0.0
pft_tab(ft)%croprange(1)=0.0
pft_tab(ft)%croprange(2)=0.0
pft_tab(ft)%croprange(3)=0.0
pft_tab(ft)%croprange(4)=0.0
pft_tab(ft)%cropphen(1)=0.0
pft_tab(ft)%cropphen(2)=0.0
pft_tab(ft)%cropphen(3)=0.0
pft_tab(ft)%cropphen(4)=0.0
pft_tab(ft)%cropphen(5)=0.0
pft_tab(ft)%cropphen(6)=0.0
pft_tab(ft)%irrig(1)=0.0
pft_tab(ft)%irrig(2)=0.0
pft_tab(ft)%irrig(3)=0.0
pft_tab(ft)%sowday(1)=0
pft_tab(ft)%sowday(2)=0
pft_tab(ft)%sowday(3)=0
pft_tab(ft)%cropgdd(1,1)=0
pft_tab(ft)%cropgdd(1,2)=0
pft_tab(ft)%cropgdd(1,3)=0
pft_tab(ft)%cropgdd(2,1)=0
pft_tab(ft)%cropgdd(2,2)=0
pft_tab(ft)%cropgdd(2,3)=0
pft_tab(ft)%fert(1)=0.
pft_tab(ft)%fert(2)=0.
pft_tab(ft)%fert(3)=0.
pft_tab(ft)%fert(4)=0.
pft_tab(ft)%fert(5)=0.
pft_tab(ft)%fert(6)=0.
pft_tab(ft)%optlai=0.
pft_tab(ft)%harvindx=0.
pft_tab(ft)%limdharv=0

!----------------------------------------------------------------------!
! Read in functional type parameterisation.                            !
!----------------------------------------------------------------------!
do ft=1,inp%npft

 pft_tab(ft)%itag=ft

 pft_tab(ft)%tag = inp%pft(ft)%tag
 pft_tab(ft)%c3c4 = inp%pft(ft)%c3c4
 pft_tab(ft)%phen = inp%pft(ft)%phen
 pft_tab(ft)%mix = inp%pft(ft)%mix
 pft_tab(ft)%crop = inp%pft(ft)%crop
 pft_tab(ft)%d2h = inp%pft(ft)%d2h
 pft_tab(ft)%mort = inp%pft(ft)%mort
 pft_tab(ft)%wden = inp%pft(ft)%wden
 pft_tab(ft)%xyl = inp%pft(ft)%xyl
 pft_tab(ft)%pdif = inp%pft(ft)%pdif
 pft_tab(ft)%sla = inp%pft(ft)%sla
 pft_tab(ft)%lls = inp%pft(ft)%lls
 pft_tab(ft)%sls = inp%pft(ft)%sls
 pft_tab(ft)%rls = inp%pft(ft)%rls
 pft_tab(ft)%lmor = inp%pft(ft)%lmor
 pft_tab(ft)%lrat = inp%pft(ft)%lrat
 pft_tab(ft)%bbmem = inp%pft(ft)%bbmem
 pft_tab(ft)%bb0 = inp%pft(ft)%bb0
 pft_tab(ft)%bbmax = inp%pft(ft)%bbmax
 pft_tab(ft)%bblim = inp%pft(ft)%bblim
 pft_tab(ft)%senm = inp%pft(ft)%senm
 pft_tab(ft)%sens = inp%pft(ft)%sens
 pft_tab(ft)%senlim = inp%pft(ft)%senlim
 pft_tab(ft)%stemx = inp%pft(ft)%stemx
 pft_tab(ft)%gr0 = inp%pft(ft)%gr0
 pft_tab(ft)%grf = inp%pft(ft)%grf
 pft_tab(ft)%ppm0 = inp%pft(ft)%ppm0

 pft_tab(ft)%can_clump = inp%pft(ft)%can_clump
 pft_tab(ft)%vna = inp%pft(ft)%vna
 pft_tab(ft)%vnb = inp%pft(ft)%vnb
 pft_tab(ft)%jva = inp%pft(ft)%jva
 pft_tab(ft)%jvb = inp%pft(ft)%jvb
 pft_tab(ft)%g0 = inp%pft(ft)%g0
 pft_tab(ft)%g1 = inp%pft(ft)%g1

 pft_tab(ft)%sowthresh(1) = inp%pft(ft)%sowthresh(1)
 pft_tab(ft)%sowthresh(2) = inp%pft(ft)%sowthresh(2)
 pft_tab(ft)%lethal(1) = inp%pft(ft)%lethal(1)
 pft_tab(ft)%lethal(2) = inp%pft(ft)%lethal(2)
 pft_tab(ft)%cardinal(1) = inp%pft(ft)%cardinal(1)
 pft_tab(ft)%cardinal(2) = inp%pft(ft)%cardinal(2)
 pft_tab(ft)%cardinal(3) = inp%pft(ft)%cardinal(3)
 pft_tab(ft)%cardinal(4) = inp%pft(ft)%cardinal(4)
 pft_tab(ft)%cardinal(5) = inp%pft(ft)%cardinal(5)
 pft_tab(ft)%cardinal(6) = inp%pft(ft)%cardinal(6)
 pft_tab(ft)%cardinal(7) = inp%pft(ft)%cardinal(7)
 pft_tab(ft)%cardinal(8) = inp%pft(ft)%cardinal(8)
 pft_tab(ft)%cardinal(9) = inp%pft(ft)%cardinal(9)
 pft_tab(ft)%croptype(1) = inp%pft(ft)%croptype(1)
 pft_tab(ft)%croptype(2) = inp%pft(ft)%croptype(2)
 pft_tab(ft)%photoperiod(1) = inp%pft(ft)%photoperiod(1)
 pft_tab(ft)%photoperiod(2) = inp%pft(ft)%photoperiod(2)
 pft_tab(ft)%photoperiod(3) = inp%pft(ft)%photoperiod(3)
 pft_tab(ft)%photoperiod(4) = inp%pft(ft)%photoperiod(4)
 pft_tab(ft)%photoperiod(5) = inp%pft(ft)%photoperiod(5)
 pft_tab(ft)%photoperiod(6) = inp%pft(ft)%photoperiod(6)
 pft_tab(ft)%croprange(1) = inp%pft(ft)%croprange(1)
 pft_tab(ft)%croprange(2) = inp%pft(ft)%croprange(2)
 pft_tab(ft)%croprange(3) = inp%pft(ft)%croprange(3)
 pft_tab(ft)%croprange(4) = inp%pft(ft)%croprange(4)
 pft_tab(ft)%cropphen(1) = inp%pft(ft)%cropphen(1)
 pft_tab(ft)%cropphen(2) = inp%pft(ft)%cropphen(2)
 pft_tab(ft)%cropphen(3) = inp%pft(ft)%cropphen(3)
 pft_tab(ft)%cropphen(4) = inp%pft(ft)%cropphen(4)
 pft_tab(ft)%cropphen(5) = inp%pft(ft)%cropphen(5)
 pft_tab(ft)%cropphen(6) = inp%pft(ft)%cropphen(6)
 pft_tab(ft)%irrig(1) = inp%pft(ft)%irrig(1)
 pft_tab(ft)%irrig(2) = inp%pft(ft)%irrig(2)
 pft_tab(ft)%irrig(3) = inp%pft(ft)%irrig(3)
 pft_tab(ft)%sowday(1) = inp%pft(ft)%sowday(1)
 pft_tab(ft)%sowday(2) = inp%pft(ft)%sowday(2)
 pft_tab(ft)%sowday(3) = inp%pft(ft)%sowday(3)
 pft_tab(ft)%cropgdd(1,1) = inp%pft(ft)%cropgdd(1,1)
 pft_tab(ft)%cropgdd(1,2) = inp%pft(ft)%cropgdd(1,2)
 pft_tab(ft)%cropgdd(1,3) = inp%pft(ft)%cropgdd(1,3)
 pft_tab(ft)%cropgdd(2,1) = inp%pft(ft)%cropgdd(2,1)
 pft_tab(ft)%cropgdd(2,2) = inp%pft(ft)%cropgdd(2,2)
 pft_tab(ft)%cropgdd(2,3) = inp%pft(ft)%cropgdd(2,3)
 pft_tab(ft)%fert(1) = inp%pft(ft)%fert(1)
 pft_tab(ft)%fert(2) = inp%pft(ft)%fert(2)
 pft_tab(ft)%fert(3) = inp%pft(ft)%fert(3)
 pft_tab(ft)%fert(4) = inp%pft(ft)%fert(4)
 pft_tab(ft)%fert(5) = inp%pft(ft)%fert(5)
 pft_tab(ft)%fert(6) = inp%pft(ft)%fert(6)
 pft_tab(ft)%optlai = inp%pft(ft)%optlai
 pft_tab(ft)%harvindx = inp%pft(ft)%harvindx
 pft_tab(ft)%limdharv = inp%pft(ft)%limdharv

 if (pft_tab(ft)%sla<0.0) then
   pft_tab(ft)%sla = 10.0**(2.35 - 0.39*log10(real(pft_tab(ft)%lls)/30.0))*2.0/10000.0
!      pft_tab(ft)%sla = 10.0**(2.43-&
! 0.46*log10(real(pft_tab(ft)%lls)/30.0))*2.0/10000.0
 endif

 pft_tab(ft)%sla = pft_tab(ft)%sla/tgp%p_sla

 if (pft_tab(ft)%lls<0.0) then
   pft_tab(ft)%lls = int(10.0**((2.35 - log10(pft_tab(ft)%sla*&
 10000.0/2.0))/0.39)*30.0+0.5)
 endif

enddo

nft = inp%npft

!----------------------------------------------------------------------!
! Open optional tile and pft output files.                             !
!----------------------------------------------------------------------!
!call open_tile_pft(nft)

!----------------------------------------------------------------------!
! Change the units of xylem and water potential difference.            !
!----------------------------------------------------------------------!
do ft=1,inp%npft
  pft_tab(ft)%xyl = pft_tab(ft)%xyl*1.0e-9
  pft_tab(ft)%pdif = pft_tab(ft)%pdif*1.0e+3
enddo

!----------------------------------------------------------------------!
! Create leaf mortality scales values.                                 !
!----------------------------------------------------------------------!
lmor_sc(1,1) = 0.0_dp
lmor_sc(2,1) = 0.0_dp
do ft=3,nft
  do i=1,pft_tab(ft)%lls
    lmor_sc(i,ft)=(real(pft_tab(ft)%lls-i)/&
 real(pft_tab(ft)%lls))**pft_tab(ft)%lmor
  enddo
enddo


!----------------------------------------------------------------------!
! Read in landuse index mapping.                                       !
!----------------------------------------------------------------------!
do ft=1,nft
  fttags(ft) = pft_tab(ft)%tag
enddo

lutab = 0.0
do i=1,inp%npft_mapping
  do j=1,inp%pft_mapping(i)%n
    k = 0
    do
      k=k+1
      if (inp%pft(k)%tag == inp%pft_mapping(i)%pft(j)) then
        lutab(inp%pft_mapping(i)%i,k) = inp%pft_mapping(i)%percent(j)
        exit
      endif
      if (k==inp%npft) then
        write(*,*) 'Error can''t match pft mapping.'
        write(*,*) inp%pft_mapping(i)%pft(j)
      endif
    enddo
  enddo
! Set bare ground to any empty space.
  lutab(inp%pft_mapping(i)%i,1) = 100.0
  do j=2,inp%npft
    lutab(inp%pft_mapping(i)%i,1) = lutab(inp%pft_mapping(i)%i,1) - lutab(inp%pft_mapping(i)%i,j)
  enddo
enddo

!----------------------------------------------------------------------!
! Grass reclimation, Bare reclimation and fire resistance (age in years)
grassrc = 0.05
barerc = 0.5
fireres = 1000

!----------------------------------------------------------------------!
! Read in sites.                                                       !
!----------------------------------------------------------------------!
latdel = 0.0
londel = 0.0

if (trim(inp%sites%style)=='list') then
  sites = inp%sites%n
  do i=1,inp%sites%n
    lat_lon(i,1) = inp%sites%list(i,1)
    lat_lon(i,2) = inp%sites%list(i,2)
  enddo
elseif (trim(inp%sites%style)=='box') then
  lat0 = inp%sites%box(1)
  latf = inp%sites%box(2)
  lon0 = inp%sites%box(3)
  lonf = inp%sites%box(4)
  latdel = inp%sites%res(1)
  londel = inp%sites%res(2)
  nlat = int((latf-lat0)/latdel + .5)
  nlon = int((lonf-lon0)/londel + .5)
  sites = 0
  do i=1,nlat
    do j=1,nlon
      sites = sites + 1
      lat_lon(sites,1) = lat0 + (i-0.5)*latdel
      lat_lon(sites,2) = lon0 + (j-0.5)*londel
    enddo
  enddo

elseif (trim(inp%sites%style)=='file') then

  site0 = inp%sites%site0
  sitef = inp%sites%sitef

  recl1 = 17

  open(99,FILE=trim(inp%sites%filename),STATUS='OLD',FORM='FORMATTED', &
 ACCESS='DIRECT',RECL=recl1,iostat=kode)
  if (kode/=0) then
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) 'File containing list of sites does not exist.'
    write(*,'('' "'',A,''"'')') trim(inp%sites%filename)
    stop
  endIF
  do site=site0,sitef
    READ(99,'(F7.3,F9.3)',REC=site,iostat=kode) &
 lat_lon(site-site0+1,1),lat_lon(site-site0+1,2)
    if (kode/=0) then
      exit
    endif
  endDO

  if (kode/=0) then
    sites = site - site0
    sitef = site0 + sites - 1
  else
    sites = site - site0
  endif
  latdel = inp%sites%res(1)
  londel = inp%sites%res(2)

else
  write(*,*) 'Sites Style must be list, box, or file.'
  stop
endif

!----------------------------------------------------------------------!
! Read in soil characteristics. Defaults used when zero                !
!----------------------------------------------------------------------!
soil_chr(1) = inp%soil%sand
soil_chr(2) = inp%soil%silt
soil_chr(3) = inp%soil%bulk
soil_chr(4) = inp%soil%orgc
soil_chr(5) = inp%soil%wilt
soil_chr(6) = inp%soil%field
soil_chr(7) = inp%soil%sat
soil_chr(8) = inp%soil%depth

if (abs(inp%soil%topsl)<1.0e-6) then
  topsl = defaulttopsl
else
  topsl = inp%soil%topsl
endif


!----------------------------------------------------------------------!
! Open diagnostics file.                                               !
!----------------------------------------------------------------------!
call OPEN_DIAG()

!----------------------------------------------------------------------!
! Check sites against land mask and disregard when no land.            !
!----------------------------------------------------------------------!
call LAND_SITE_CHECK(st1,sites,lat_lon,latdel,londel,du)

!----------------------------------------------------------------------!
! Read in type of landuse: 0 = defined by map; 1 = defined explicitly  !
! in the input file; 2 = natural vegetation based on average monthly   !
! temperatures.                                                        !
!----------------------------------------------------------------------!
if (inp%land_use%read_from_landuse_dir) then
  ilanduse = 0
else
  ilanduse = 1
endif

fire_ant(:)    = .FALSE.
harvest_ant(:) = .FALSE.
if (ilanduse==1) then
! Use landuse defined in input file.
  call LANDUSE1(luse,yr0,yrf,fire_ant,harvest_ant)
endif
if ((ilanduse<0).or.(ilanduse>2)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'No landuse defined'
  write(*,*) '0:= Defined from a map.'
  write(*,*) '1:= Defined in the input file.'
  write(*,*) '2:= Natural vegetation.'
  stop
endif
call OPEN_SNAPSHOTS(snp_no,snpshts)

ssp%nft = nft

end subroutine process_input_file






!**********************************************************************!
!                                                                      !
!                    read_param :: read_input                          !
!                    ------------------------                          !
!                                                                      !
! subroutine read_param(l_regional,site_out,year_out,stver)            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details First reads misc_params.dat file in /f03 with input on year
!! step to output on screen and whether run is regional or country.Saved
!! in msp structure defined in misc_parameters.f90.
!! It then reads the param.dat file in /inc with the tuning parameters
!! and saves in structure tgp defined in tuning_parameters.f90.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_param(stver)
!**********************************************************************!
integer :: kode,fid
character(len=str_len) :: st1,st2,st3,stver
logical :: l_regional

!----------------------------------------------------------------------!
! Read 'misc_params.dat'.                                              !
!----------------------------------------------------------------------!
open(newunit=fid,file='misc_params.dat',status='OLD',iostat=kode)

if (kode/=0) then
  write(*,*) ' File does not exist: "misc_params.dat"'
  write(*,*) ' Using screen output options: 1 0. '
  write(*,*) ' Using countries, not regions. '
  msp%site_out = 1
  msp%year_out = 0
  msp%l_regional = .false.
else
  read(fid,*)
  read(fid,*) msp%site_out,msp%year_out
  read(fid,*)
  read(fid,*) msp%l_regional
endif
close(fid)

!----------------------------------------------------------------------!
! Read internal parameters from "param.dat".                           !
!----------------------------------------------------------------------!
open(newunit=fid,file='inc/param.dat',status='OLD',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) ' File does not exist: "param.dat"'
  stop
endif

read(fid,*)
read(fid,*) st1
call STRIPB(st1)
st2 = stver
call STRIPBS(st2,st3)
call STRIPB(st2)
st1=st1(1:blank(st1))
st2=st2(1:blank(st2))
if (stcmp(st1,st2)/=1) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) &
 'Version number mismatch between ''sdgvm0.f'' and ''param.dat''.'
  stop
endif

read(fid,*)
read(fid,*) tgp%p_sla        ! 4
read(fid,*)
read(fid,*) tgp%p_stemfr     ! 6
read(fid,*)
read(fid,*) tgp%p_rootfr     ! 8
read(fid,*)
read(fid,*) tgp%p_opt        !10
read(fid,*)
read(fid,*) tgp%p_stmin      !12
read(fid,*)
read(fid,*) tgp%p_laimem     !14
read(fid,*)
read(fid,*) tgp%p_resp       !16
read(fid,*)
read(fid,*) tgp%p_kscale     !18
read(fid,*)
read(fid,*) tgp%p_nu1,tgp%p_nu2,tgp%p_nu3,tgp%p_nu4
read(fid,*)
read(fid,*) tgp%p_nleaf
read(fid,*)
read(fid,*) tgp%p_dresp
read(fid,*)
read(fid,*) tgp%p_vm
read(fid,*)
read(fid,*) tgp%p_kgw
read(fid,*)
read(fid,*) tgp%p_v1,tgp%p_v2,tgp%p_v3
read(fid,*)
read(fid,*) tgp%p_j1,tgp%p_j2,tgp%p_j3
read(fid,*)
read(fid,*) tgp%p_pet
read(fid,*)
read(fid,*) tgp%p_bs
read(fid,*)
read(fid,*) tgp%p_et
read(fid,*)
read(fid,*) tgp%p_roff
read(fid,*)
read(fid,*) tgp%p_roff2
read(fid,*)
read(fid,*) tgp%p_fprob
read(fid,*)
read(fid,*) tgp%p_city_dep
close(fid)

end subroutine read_param



end module read_input
