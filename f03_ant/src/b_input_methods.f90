module input_methods

use real_precision
use dims
use pft_parameters
use site_parameters
use output_methods
use func
use screen_output_parameters
use tuning_parameters
use open_files
use file_class
use file_object

implicit none

contains

subroutine command_line_argument_check()
logical :: arg_input_not_set,arg_output_not_set,arg_n_range_1_not_set
logical :: arg_n_range_2_not_set,arg_filename_not_set
integer :: additional_arguments,i
character(len=200) :: arg
!----------------------------------------------------------------------!
! Check additional command line arguemnts are consistent with input    !
! and set.                                                             !
!----------------------------------------------------------------------!
arg_input_not_set     = .true.
arg_output_not_set    = .true.
arg_filename_not_set  = .true.
arg_n_range_1_not_set = .true.
arg_n_range_2_not_set = .true.

additional_arguments = 0
if (inp%dirs%input_argument) additional_arguments = &
 additional_arguments + 1
if (inp%dirs%output_argument) additional_arguments = &
 additional_arguments + 1
if (inp%sites%filename_argument) additional_arguments = &
 additional_arguments + 1
if (inp%sites%n_range_argument) additional_arguments = &
 additional_arguments + 2

i = 0
do
call get_command_argument(i,arg)
  if (len_trim(arg) == 0) exit
  i = i + 1
enddo

! check number of additoinal args is as required by the input file.
if (additional_arguments/=i-2) then
  write(*,'(''There are'',i2,'' arguments to sdgvm,'',i2 &
 ,'' were expected.'')') i-1,additional_arguments+1
  write(*,'(''Expected'')')
  i = 1
  write(*,'(i2,''. input file'')') i
  if (inp%dirs%input_argument) then
    i = i + 1
    write(*,'(i2''. input directory'')') i
  endif
  if (inp%dirs%output_argument) then
    i = i + 1
    write(*,'(i2''. output directory'')') i
  endif
  if (inp%sites%filename_argument) then
    i = i + 1
    write(*,'(i2''. land sites file directory'')') i
  endif
  if (inp%sites%n_range_argument) then
    i = i + 1
    write(*,'(i2,''. site0'')') i
  endif
  if (inp%sites%n_range_argument) then
    i = i + 1
    write(*,'(i2''. sitef'')') i
  endif
  stop
endif

i = 0
do
call get_command_argument(i,arg)
  if (len_trim(arg) == 0) exit
  if (i>1) then
    if ((inp%dirs%input_argument).and.(arg_input_not_set)) then
      read(arg,'(A)') inp%dirs%input
      arg_input_not_set = .false.
    elseif ((inp%dirs%output_argument).and.(arg_output_not_set)) then
      read(arg,'(A)') inp%dirs%output
      arg_output_not_set = .false.
    elseif ((inp%sites%filename_argument).and.(arg_filename_not_set)) then
      read(arg,'(A)') inp%sites%filename
      arg_filename_not_set = .false.
    elseif ((inp%sites%n_range_argument).and.(arg_n_range_1_not_set)) then
      read(arg,*) inp%sites%site0
      arg_n_range_1_not_set = .false.
    elseif ((inp%sites%n_range_argument).and.(arg_n_range_2_not_set)) then
      read(arg,*) inp%sites%sitef
      arg_n_range_2_not_set = .false.
    endif
  endif
  i = i + 1
enddo

end subroutine command_line_argument_check





!**********************************************************************!
!                                                                      !
!                    process_input_file :: read_input                  !
!                    -----------------------------                     !
!                                                                      !
! SUBROUTINE process_input_file(buff1,stco2,xlatf,xlon0,xlatres,xlonres,  !
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
 site0,sitef,luse(max_years),ilanduse,persum
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
        CALL stripbs(st1,st2)
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

if (stcmp(st1,st2)==1) then
  !1 for daily data 0 for monthly
  day_mnth = 1
  !0 for daily site run 1 for all else
  thty_dys = 1
  !1 for site run 0 for non-site run
  sit_grd = 0
  read(fid,*)
  !Upper left lat and long
  read(fid,*) xlatf,xlon0
  read(fid,*)
  !lat and lon resolution
  read(fid,*) xlatres,xlonres
  read(fid,*)
  !Number of pixels (north to south, east to west)
  read(fid,*) xlatresn,xlonresn
  read(fid,*)
  !Initial and final year of the climate dataset
  read(fid,*) xyear0,xyearf
elseif (stcmp(st1,st3)==1) then
  day_mnth = 0
  thty_dys = 1
  sit_grd = 0
  read(fid,*)
  read(fid,*) xlatf,xlon0
  read(fid,*)
  read(fid,*) xlatres,xlonres
  read(fid,*)
  read(fid,*) xlatresn,xlonresn
  read(fid,*)
  read(fid,*) xyear0,xyearf
elseif (stcmp(st1,st4)==1) then
  day_mnth = 1
  thty_dys = 0
  sit_grd = 1
  read(fid,*)
  read(fid,*) xlatf,xlon0
  read(fid,*)
  read(fid,*) xyear0,xyearf
elseif (stcmp(st1,st5)==1) then
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



if(inp%run%subdaily) then
  write(*,*) ''
  write(*,'('' SDGVM running with '',i3, &
    '' sub-daily photosynthesis time-points'')') par_loops
else
! write(*,*) ''
! write(*,*) 'SDGVM running with 1 sub-daily photosynthesis time-points'
endif


      
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
xspeedc = speedc


!Set the number of years 'nyears' for the run.If starting year and final year
!not set properly then it will set nyears to be spinup
if ((inp%run%yearf<inp%run%year0).or.(inp%run%yearf.eq.0).or. &
 (inp%run%year0.eq.0)) then
  inp%run%year0 = 0
  inp%run%yearf = 0
  nyears = inp%run%spinup_length
else
  nyears = inp%run%spinup_length + inp%run%yearf - inp%run%year0 + 1
endif

! Set nyears if it hasn't been set.
if (inp%output%nyears.eq.0) then
  if (inp%run%year0.gt.0) then
    inp%output%nyears = inp%run%yearf - inp%run%year0 + 1
  else
    inp%output%nyears = inp%run%spinup_cycle_length + 1
  endif
endif

!Assigns the total run years
outyears = min(nyears,inp%output%nyears)

!It stops if the run exceeds 'max_years' which is set to 1000 in dims.f90.
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
if (inp%run%year0>0) then
  yr0 = min(inp%run%spinup_year0,inp%run%year0)
  yrf = max(inp%run%spinup_year0+min(inp%run%spinup_length,inp%run%spinup_cycle_length)-1,inp%run%yearf)
else
  yr0 = inp%run%spinup_year0
  yrf = inp%run%spinup_year0+min(inp%run%spinup_length,inp%run%spinup_cycle_length)-1
endif

yr0m = yr0
yrfm = yrf

if(met_seq) yr0m = min(yr0s,yr0ms)
if(met_seq) yrfm = max(yr0s+min(inp%run%spinup_length,inp%run%spinup_cycle_length)-1,yrfms)

!Checks if the first and last year required by the run are in the climate data
!It won't exit even if they arent in the climate data
if ((yr0m<xyear0).or.(yrfm>xyearf)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,'('' Trying to use '',i4,''-'',i4,'' climate.'')') yr0m,yrfm 
  write(*,'('' Climate database runs from '',i4,''-'',i4,''.'')') xyear0,xyearf
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
      !Assigns the calendar year for each year of spinup.Here is just cycles
      !through the years of the spinup period
      yearv(i) = mod(i-1,inp%run%spinup_cycle_length) + inp%run%spinup_year0
    endif
  else
    !Assigns the calendar year for each year of transient run.
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
call output_options(nomdos,otags,omav,ofmt_daily,ofmt_monthly,ofmt_yearly)

call set_pixel_out(st1,st2,outyears1,nomdos,otagsn,otags,oymd)
outyears1 = min(outyears1,nyears)

st1 = inp%dirs%output

!----------------------------------------------------------------------!
! Determine whether daily or monthly subpixel outputs are required.    !
!----------------------------------------------------------------------!

call set_subpixel_out(st1,st2,outyears2,nomdos,otagsnft,otags,oymdft,&
 out_cov,out_bio,out_bud,out_sen)
outyears2 = min(outyears2,nyears)

!----------------------------------------------------------------------!
! Read in compulsory functional types.                                 !
!----------------------------------------------------------------------!
pft_tab(1)%tag = 'BARE'
pft_tab(2)%tag = 'CITY'
pft_tab(1)%itag = 1
pft_tab(2)%itag = 2

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
!----------------------------------------------------------------------!
grassrc = 0.05
barerc = 0.5
fireres = 1000

!----------------------------------------------------------------------!
! Read in sites.                                                       !
!----------------------------------------------------------------------!
latdel = 0.0
londel = 0.0

if (trim(inp%sites%style)=='list') then
  if (inp%sites%sitef>inp%sites%nlist)  inp%sites%sitef = inp%sites%nlist
  sites = inp%sites%sitef - inp%sites%site0 + 1
  do i=1,sites
    lat_lon(i,1) = inp%sites%list(i+inp%sites%site0-1,1)
    lat_lon(i,2) = inp%sites%list(i+inp%sites%site0-1,2)
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
  inp%sites%site0 = 1
  inp%sites%sitef = sites

elseif (trim(inp%sites%style)=='file') then
  site0 = inp%sites%site0
  sitef = inp%sites%sitef

  recl1 = 17

  open(99,file=trim(inp%sites%filename),STATUS='OLD',FORM='FORMATTED', &
 ACCESS='DIRECT',RECL=recl1,iostat=kode)
  if (kode/=0) then
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) 'File containing list of sites does not exist.'
    write(*,'('' "'',A,''"'')') trim(inp%sites%filename)
    stop
  endIF
  do site=site0,sitef
    read(99,'(F7.3,F9.3)',rec=site,iostat=kode) &
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
call open_diag()

!----------------------------------------------------------------------!
! Check sites against land mask and disregard when no land.            !
!----------------------------------------------------------------------!
call land_site_check(st1,sites,lat_lon,latdel,londel,du)

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

fire_ant(:)    = .false.
harvest_ant(:) = .false.
if (ilanduse==1) then
! Use landuse defined in input file.
  call landuse1(luse,yr0,yrf,fire_ant,harvest_ant)
endif
if ((ilanduse<0).or.(ilanduse>2)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'No landuse defined'
  write(*,*) '0:= Defined from a map.'
  write(*,*) '1:= Defined in the input file.'
  write(*,*) '2:= Natural vegetation.'
  stop
endif
call open_snapshots(snp_no,snpshts)

ssp%nft = nft

end subroutine process_input_file






!**********************************************************************!
!                                                                      !
!                    read_param :: read_input                          !
!                    ------------------------                          !
!                                                                      !
! subroutine read_param(stver)                                         !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "inc/param.dat" file, and io
!! parameters from "inc/screen_output.dat.dat".
!! @details First reads screen_output.dat which holds parameters on how
!! to output on screen during run.
!! Parameters are saved in sop structure defined in screen_output_parameters.f90
!! with more details.
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
! Read 'inc/screen_output.dat'.                                        !
!----------------------------------------------------------------------!
open(newunit=fid,file='inc/screen_output.dat',status='OLD',iostat=kode)
if (kode/=0) then
  write(*,*) ' File does not exist: "screen_output.dat"'
  write(*,*) ' Using screen output options: 1 0. '
  write(*,*) ' Using countries, not regions. '
  sop%site_out = 1
  sop%year_out = 0
  sop%l_regional = .false.
else
  read(fid,*)
  read(fid,*) sop%site_out,sop%year_out
  read(fid,*)
  read(fid,*) sop%l_regional
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
st1 = adjustl(st1)
st2 = stver
call stripbs(st2,st3)
st2 = adjustl(st2)
st1 = trim(st1)
st2 = trim(st2)
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





!**********************************************************************!
!                                                                      !
!                       land_site_check :: read_input                  !
!                       -----------------------------                  !
!                                                                      !
! subroutine land_site_check(st1,st2,sites,lat_lon,latdel,londel,du,   !
! stmask)                                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine land_site_check(st1,sites,lat_lon,latdel,londel,du)
!**********************************************************************!
character(len=str_len) :: st1
integer :: sites,du,i,site,fid
real(dp) :: lat_lon(max_sites,2),latdel,londel,lat,lon
logical :: lc
!----------------------------------------------------------------------!

call fun%open(trim(inp%dirs%output)//'/land_sites.dat',fid)
i = 0
write(*,'(''Number of potential sites = '',i0)') sites
write(*,'(''Land sites'')')
do site=1,sites
  lat = lat_lon(site,1)
  lon = lat_lon(site,2)
  call lorc(du,lat,lon,latdel,londel,lc)
  if (lc) then
  if (mod(site,1)==0)  write(*,'(i6,f8.3,f9.3)') site,lat,lon
    write(fid,'(f7.3,f9.3)') lat,lon
    i = i + 1
    lat_lon(i,1) = lat_lon(site,1)
    lat_lon(i,2) = lat_lon(site,2)
  endif
enddo
call fun%close(fid)

sites = i

write(*,'(''Number of land sites = '',i0)') sites
if (inp%sites%create_land_sites_only) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Created land sites file only.'
  stop
endif

end subroutine land_site_check





!**********************************************************************!
!                                                                      !
!                     lorc :: data                                     !
!                     ------------                                     !
!                                                                      !
! subroutine lorc(du,lat,lon,latdel,londel,xx)                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Determine if the site is a land site
!! @details Using a half an arc-second land sea mask. Determine whether
!! the site is land or sea based on the majority of the mask within
!! the limits of the gridcell.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine lorc(du,lat,lon,latdel,londel,xx)
!**********************************************************************!
logical :: xx
real(dp) :: lat,lon,latdel,londel,del,latf,lon0
character :: outc(6200)
integer :: i,j,k,x(7),ians(43300),sum1,n,col,row,check,nrecl,du,kode
!----------------------------------------------------------------------!

if (du==1) then
!  This works for ftn95
!  nrecl = 6172
  nrecl = 6172 + 2
else
  nrecl = 6172 + 1
endif

open(99,file=trim(inp%dirs%land_mask)//'/land_mask.dat', &
 form='formatted',recl=nrecl,access='direct',status='old',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Land sea mask.'
  write(*,*) 'Either the file doesn''t exist or there is a record length miss match.'
  write(*,*) 'Check that the correct DOS|UNIX switch is being used in the input file.'
  write(*,*) 'Land sea file :',trim(inp%dirs%land_mask)
  stop
endif

del = 1.0d0/60.0d0/2.0d0
latf = 90.0d0 - del/2.0d0
lon0 =-180.0d0 + del/2.0d0

col = int((lon - lon0)/del + 0.5d0)
row = int((latf - lat)/del + 0.5d0)

n = min((latdel/del-1.0d0)/2.0d0,(londel/del-1.0d0)/2.0d0)

!----------------------------------------------------------------------!
! Check nearest pixel for land.                                        !
!----------------------------------------------------------------------!
sum1 = 0
read(99,'(6172a)',rec=row) (outc(j),j=1,6172)

do j=1,6172
  call base72i(outc(j),x)
  do k=1,7
    ians(7*(j-1)+k) = x(k)
  enddo
enddo
if (ians(col)>0) sum1 = sum1 + 1
check = 1

!----------------------------------------------------------------------!
if (n>0) then
!----------------------------------------------------------------------!
! Check outward diagonals for land.                                    !
!----------------------------------------------------------------------!
  do i=1,n
    read(99,'(6172a)',rec=row+i) (outc(j),j=1,6172)
    do j=1,6172
      call base72i(outc(j),x)
      do k=1,7
        ians(7*(j-1)+k) = x(k)
      enddo
    enddo
    if (ians(col+i)>0) sum1 = sum1 + 1
    if (ians(col-i)>0) sum1 = sum1 + 1

    read(99,'(6172a)',rec=row-i) (outc(j),j=1,6172)
    do j=1,6172
      call base72i(outc(j),x)
      do k=1,7
        ians(7*(j-1)+k) = x(k)
      enddo
    enddo
    if (ians(col+i)>0) sum1 = sum1 + 1
    if (ians(col-i)>0) sum1 = sum1 + 1
    check = check + 4
  enddo
!----------------------------------------------------------------------!
endif

close(99)

if (real(sum1)>=(check+1)/2) then
  xx = .true.
else
  xx = .false.
endif

end subroutine lorc





!**********************************************************************!
!                                                                      !
!                          base72I :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine base72i(c,x)                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine base72i(c,x)
!**********************************************************************!
character :: c
integer :: x(7),i
!----------------------------------------------------------------------!

i=ichar(c)-100
x(1)=i/64
i=i-x(1)*64
x(2)=i/32
i=i-x(2)*32
x(3)=i/16
i=i-x(3)*16
x(4)=i/8
i=i-x(4)*8
x(5)=i/4
i=i-x(5)*4
x(6)=i/2
i=i-x(6)*2
x(7)=i

end subroutine base72i





!**********************************************************************!
!                                                                      !
!                          n7 :: data                                  !
!                          ----------                                  !
!                                                                      !
! Returns the value of a seven digit binary number given as a seven    !
! dimensional array of 0's and 1's                                     !
!                                                                      !
!    integer function n7(x)                                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
integer function n7(x)
!----------------------------------------------------------------------!
integer :: x(7)
!----------------------------------------------------------------------!

n7 = 64*x(1)+32*x(2)+16*x(3)+8*x(4)+4*x(5)+2*x(6)+x(7)
 
!**********************************************************************!
end function n7





!**********************************************************************!
!                                                                      !
!                     get_input_filename :: sdgvm1                     !
!                     ----------------------------                     !
!                                                                      !
! subroutine get_input_filename(st1)                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Reads the 2nd argument from the command line which should be
!! the path of the input file
!! @details EPK changed the getarg to get_command_argument
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine get_input_filename(st1)
!**********************************************************************!
character :: st1*80
integer :: iargc
!----------------------------------------------------------------------!

if (IARGC()>0) then
  call get_command_argument(1,st1)
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) ' Input file must be given as an argument.'
  stop
endif

end subroutine get_input_filename







!**********************************************************************!
!                                                                      !
!                       read_climate :: sdgvm1                         !
!                       ----------------------                         !
!                                                                      !
! subroutine read_climate(ststats,lat,lon,xlatf,                       !
! xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,   !
! xcldv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats,siteno,!
! day_mnth,thty_dys,sit_grd,withcloudcover))                           !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_climate(lat,lon,xlatf,xlatres,xlatresn,xlon0,xlonres, &
 xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,xcldv,xswrv,isite,xyear0,xyearf, &
 du,seed1,seed2,seed3,l_clim,l_stats,siteno,day_mnth,thty_dys,sit_grd, &
 withcloudcover)
!**********************************************************************!
real(dp) :: lat,lon,xlatf,xlatres,xlon0,xlonres
integer :: xlatresn,xlonresn,yr0,yrf,isite,xyear0,xyearf,du
integer :: seed1,seed2,seed3,siteno,day_mnth,thty_dys,sit_grd
real(dp), dimension(500,12,31) :: xtmpv,xprcv,xhumv,xswrv
real(dp), dimension(500,12) :: xcldv
integer :: read_par
logical :: l_clim,l_stats,withcloudcover
!----------------------------------------------------------------------!
!Daily climate drivers
if ((day_mnth==1).and.(thty_dys==1)) then
  call EX_CLIM(lat,lon,xlatf,xlatres,xlatresn,xlon0,xlonres, &
 xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,isite,xyear0,xyearf,siteno,du, &
 read_par)
  withcloudcover=.false.
  ssp%latres = xlatres
  ssp%lonres = xlonres
!Monthly climate drivers
elseif ((day_mnth==0).and.(thty_dys==1).and.(sit_grd==0)) then
  call EX_CLIM_WEATHER_GENERATOR(lat,lon,xlatf,xlatres,xlatresn,xlon0, &
 xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,xcldv,xswrv,isite,xyear0, &
 xyearf,du,seed1,seed2,seed3,l_clim,l_stats,inp%run%read_par)
  withcloudcover=.true.
  ssp%latres = xlatres
  ssp%lonres = xlonres
!Daily climate drivers for single site
elseif ((day_mnth==1).and.(thty_dys==0)) then
  call EX_CLIM_SITE(yr0,yrf,xtmpv,xhumv,xprcv,xswrv,xyear0,xyearf)
  withcloudcover=.false.
  siteno = 1
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Error defining climate to read'
  stop
endif

end subroutine read_climate









end module input_methods
