module output_methods

use real_precision
use dims
use data
use weather_generator
use veg_dynamics
use func
use file_class
use file_object
use input_file

implicit none

contains


!**********************************************************************!
!                                                                      !
!                     output_mapping :: output_methods                 !
!                     -------------------------------                  !
!                                                                      !
! subroutine output_mapping(nomdos,otags)                              !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief 
!! @details Checks that the variables requested for outputs exist
!!          Also assigns to the output mapping variable the id of each
!!          requested variable.nomdos is the number of default variables
!!          and otags is a list with the names of the default variables
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine output_mapping(nomdos,otags)
integer :: nomdos,i,j
character(len=str_len), dimension(max_outputs) :: otags

!----------------------------------------------------------------------!
! set up mapping for output files for and place in inp%output%*%map.   !
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! For each variable requested in the output for annual,monthly and     !     
! daily inp%output%tile_vars_*%tag(i), loops through the list of       !
! available variables (otags(j)) and assigns its id in                 !
! inp%output%tile_vars_*%map(i).                                       ! 
!----------------------------------------------------------------------!
do i=1,inp%output%tile_vars_yearly%n
  do j=1,nomdos
    if (trim(inp%output%tile_vars_yearly%tag(i)) == trim(otags(j))) then
      inp%output%tile_vars_yearly%map(i) = j
      exit
    endif
    if (j==nomdos) then
      write(*,*) 'error in pft tags for output of tile_vars_yearly'
      write(*,*) trim(inp%output%tile_vars_yearly%tag(i))
      stop
    endif
  enddo
enddo
do i=1,inp%output%tile_vars_monthly%n
  do j=1,nomdos
    if (trim(inp%output%tile_vars_monthly%tag(i)) == trim(otags(j))) then
      inp%output%tile_vars_monthly%map(i) = j
      exit
    endif
    if (j==nomdos) then
      write(*,*) 'error in pft tags for output of tile_vars_monthly'
      write(*,*) trim(inp%output%tile_vars_monthly%tag(i))
      stop
    endif
  enddo
enddo
do i=1,inp%output%tile_vars_daily%n
  do j=1,nomdos
    if (trim(inp%output%tile_vars_daily%tag(i)) == trim(otags(j))) then
      inp%output%tile_vars_daily%map(i) = j
      exit
    endif
    if (j==nomdos) then
      write(*,*) 'error in pft tags for output of tile_vars_daily'
      write(*,*) trim(inp%output%tile_vars_daily%tag(i))
      stop
    endif
  enddo
enddo

do i=1,inp%output%pft_vars_yearly%n
  do j=1,nomdos
    if (trim(inp%output%pft_vars_yearly%tag(i)) == trim(otags(j))) then
      inp%output%pft_vars_yearly%map(i) = j
      exit
    endif
    if (j==nomdos) then
      write(*,*) 'error in pft tags for output of pft_vars_yearly'
      write(*,*) trim(inp%output%pft_vars_yearly%tag(i))
      stop
    endif
  enddo
enddo
do i=1,inp%output%pft_vars_monthly%n
  do j=1,nomdos
    if (trim(inp%output%pft_vars_monthly%tag(i)) == trim(otags(j))) then
      inp%output%pft_vars_monthly%map(i) = j
      exit
    endif
    if (j==nomdos) then
      write(*,*) 'error in pft tags for output of pft_vars_monthly'
      write(*,*) trim(inp%output%pft_vars_monthly%tag(i))
      stop
    endif
  enddo
enddo
do i=1,inp%output%pft_vars_daily%n
  do j=1,nomdos
    if (trim(inp%output%pft_vars_daily%tag(i)) == trim(otags(j))) then
      inp%output%pft_vars_daily%map(i) = j
      exit
    endif
    if (j==nomdos) then
      write(*,*) 'error in pft tags for output of pft_vars_daily'
      write(*,*) trim(inp%output%pft_vars_daily%tag(i))
      stop
    endif
  enddo
enddo

end subroutine output_mapping



!**********************************************************************!
!                                                                      !
!                     write_lat_lon :: output_methods                  !
!                     -------------------------------                  !
!                                                                      !
! subroutine WRITE_LAT_LON(lat,lon)                                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Write lat & lon for output files.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine write_lat_lon(lat,lon)
!**********************************************************************!
integer :: i
real(dp):: lat,lon
character(len=str_len) :: st1
!----------------------------------------------------------------------!

do i=1,fun%n
  if (fun%opened(i)) then
    call filename_base(fun%names(i),st1)
    if ((st1(1:13)/='site_info.dat').and.(st1(1:13)/='state.dat')) then
      if ((st1(1:8)=='monthly_').or.(st1(1:6)=='daily_')) then
        write(fun%get_id_n(i),'(f7.3,f9.3)') lat,lon
      else
        write(fun%get_id_n(i),'(f7.3,f9.3)',advance='NO') lat,lon
      endif
    endif
  endif
enddo

end subroutine write_lat_lon





!**********************************************************************!
!                                                                      !
!                          outputs :: output_methods                   !
!                          -------------------------                   !
!                                                                      !
! real(dp) function outputs(daily_out,ind,st1)                         !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Compute outputs from daily_out.
!! @details Four forms of outputs are available and are specified by
!!'st1' which can take any one of the values
!! ['Add', 'Average', 'Max', 'Min'].
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
real(dp) function outputs(daily_out,ind,st1)
!**********************************************************************!
real(dp) ::  daily_out(max_outputs,max_cohorts,12,31),omax,ans
integer :: ind,co,mnth,day
character(len=3) :: st1
!----------------------------------------------------------------------!

outputs = 0.0
if ((st1(1:1) == 'A').and.(st1(2:2) == 'd')) then
  do co=1,ssp%cohorts
    ans = 0.0
    do mnth=1,12
      do day=1,30
        ans = ans + daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov
  enddo
elseif ((st1(1:1) == 'A').and.(st1(2:2) == 'v')) then
  do co=1,ssp%cohorts
    ans = 0.0
    do mnth=1,12
      do day=1,30
        ans = ans + daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov/360.0
  enddo
elseif ((st1(1:1) == 'M').or.(st1(2:2) == 'a')) then
  do co=1,ssp%cohorts
    ans = -1000000.0
    do mnth=1,12
      do day=1,30
        if (daily_out(ind,co,mnth,day) > ans) ans = daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov
  enddo
elseif ((st1(1:1) == 'M').or.(st1(2:2) == 'i')) then
  do co=1,ssp%cohorts
    ans = 1000000.0
    do mnth=1,12
      do day=1,30
        if (daily_out(ind,co,mnth,day) < ans) ans = daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov
  enddo
else
  write(*,'(''Error in the third argument to outputs in SDGVM0.'')')
  write(*,'(''Arguments: Add, Average, Max, Min.'')')
  stop
endif

end function outputs





!**********************************************************************!
!                                                                      !
!                       output_options :: output_methods               !
!                       --------------------------------               !
!                                                                      !
! subroutine output_options(nomdos,otags,omav,ofmt)                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine output_options(nomdos,otags,omav,ofmt_daily,ofmt_monthly, &
 ofmt_yearly)
!**********************************************************************!
integer, dimension(max_outputs) :: omav
integer :: nomdos
character(len=str_len), dimension(max_outputs) :: otags,ofmt_daily, &
 ofmt_monthly,ofmt_yearly

!----------------------------------------------------------------------!
! Read in monthly and/or daily output options.                         !
!----------------------------------------------------------------------!
! file tags.                                                           !
!----------------------------------------------------------------------!
otags(1)  = 'lai'
otags(2)  = 'rof'
otags(3)  = 'evt'
otags(4)  = 'trn'
otags(5)  = 'npp'
otags(6)  = 'gpp'
otags(7)  = 'srp'
otags(8)  = 'nep'
otags(9)  = 'tmp'
otags(10) = 'prc'
otags(11) = 'hum'
otags(12) = 'nps'
otags(13) = 'swf'
otags(14) = 'pet'
otags(15) = 'int'
otags(16) = 'bse'
otags(17) = 'ssm'
otags(18) = 'swc'
otags(19) = 'rsp'
otags(20) = 'qdr'
otags(21) = 'qdf'
otags(22) = 'lfn'
otags(23) = 'lfl'
otags(24) = 'cld'
otags(25) = 'fpr'

nomdos = 25

!----------------------------------------------------------------------!
! Averaged output.                                                     !
!----------------------------------------------------------------------!
omav(1) = 1
omav(2) = 0
omav(3) = 0
omav(4) = 0
omav(5) = 0
omav(6) = 0
omav(7) = 0
omav(8) = 0
omav(9) = 1
omav(10) = 0
omav(11) = 1
omav(12) = 1
omav(13) = 1
omav(14) = 0
omav(15) = 0
omav(16) = 0
omav(17) = 1
omav(18) = 1
omav(19) = 0
omav(20) = 1
omav(21) = 1
omav(22) = 1
omav(23) = 1
omav(24) = 1

!----------------------------------------------------------------------!
! Daily output format.                                                 !
!----------------------------------------------------------------------!
ofmt_daily(1)  ='(31f7.2)'
ofmt_daily(2)  ='(31f7.2)'
ofmt_daily(3)  ='(31f7.2)'
ofmt_daily(4)  ='(31f7.2)'
ofmt_daily(5)  ='(31f7.1)'
ofmt_daily(6)  ='(31f7.1)'
ofmt_daily(7)  ='(31f7.1)'
ofmt_daily(8)  ='(31f7.1)'
ofmt_daily(9)  ='(31f7.1)'
ofmt_daily(10)  ='(31f7.1)'
ofmt_daily(11)  ='(31f7.2)'
ofmt_daily(12)  ='(31f7.1)'
ofmt_daily(13)  ='(31f7.3)'
ofmt_daily(14)  ='(31f7.2)'
ofmt_daily(15)  ='(31f7.2)'
ofmt_daily(16)  ='(31f7.2)'
ofmt_daily(17)  ='(31f7.4)'
ofmt_daily(18)  ='(31f7.2)'
ofmt_daily(19)  ='(31f7.1)'
ofmt_daily(20)  ='(31f7.2)'
ofmt_daily(21)  ='(31f7.2)'
ofmt_daily(22)  ='(31f7.3)'
ofmt_daily(23)  ='(31f7.3)'
ofmt_daily(24)  ='(31f7.3)'
ofmt_daily(25)  ='(31f7.3)'

!----------------------------------------------------------------------!
! Monthly output format.                                               !
!----------------------------------------------------------------------!
ofmt_monthly(1)   ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(2)   ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(3)   ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(4)   ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(5)   ='(i4,1x,12f8.1,f10.3)'
ofmt_monthly(6)   ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(7)   ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(8)   ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(9)   ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(10)  ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(11)  ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(12)  ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(13)  ='(i4,1x,12f8.3,f10.3)'
ofmt_monthly(14)  ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(15)  ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(16)  ='(i4,1x,12f8.2,f10.1)'
ofmt_monthly(17)  ='(i4,1x,12f8.2,f10.5)'
ofmt_monthly(18)  ='(i4,1x,12f8.2,f10.3)'
ofmt_monthly(19)  ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(20)  ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(21)  ='(i4,1x,12f8.1,f10.1)'
ofmt_monthly(22)  ='(i4,1x,12f8.3,f10.1)'
ofmt_monthly(23)  ='(i4,1x,12f8.3,f10.1)'
ofmt_monthly(24)  ='(i4,1x,12f8.3,f10.1)'
ofmt_monthly(25)  ='(i4,1x,12f8.3,f10.1)'

!----------------------------------------------------------------------!
! Yearly output format.                                                 !
!----------------------------------------------------------------------!
ofmt_yearly(1)  ='(31f7.2)'
ofmt_yearly(2)  ='(31f7.2)'
ofmt_yearly(3)  ='(31f7.2)'
ofmt_yearly(4)  ='(31f7.2)'
ofmt_yearly(5)  ='(31f7.1)'
ofmt_yearly(6)  ='(31f7.1)'
ofmt_yearly(7)  ='(31f7.1)'
ofmt_yearly(8)  ='(31f7.1)'
ofmt_yearly(9)  ='(31f7.1)'
ofmt_yearly(10)  ='(31f7.1)'
ofmt_yearly(11)  ='(31f7.2)'
ofmt_yearly(12)  ='(31f7.1)'
ofmt_yearly(13)  ='(31f7.3)'
ofmt_yearly(14)  ='(31f7.2)'
ofmt_yearly(15)  ='(31f7.2)'
ofmt_yearly(16)  ='(31f7.2)'
ofmt_yearly(17)  ='(31f7.4)'
ofmt_yearly(18)  ='(31f7.2)'
ofmt_yearly(19)  ='(31f7.1)'
ofmt_yearly(20)  ='(31f7.2)'
ofmt_yearly(21)  ='(31f7.2)'
ofmt_yearly(22)  ='(31f7.3)'
ofmt_yearly(23)  ='(31f7.3)'
ofmt_yearly(24)  ='(31f7.3)'
ofmt_yearly(25)  ='(31f7.3)'

end subroutine output_options





!**********************************************************************!
!                                                                      !
!                       set_pixel_out :: output_methods                !
!                       -------------------------------                !
!                                                                      !
! subroutine set_pixel_out(st1,st2,outyears,nomdos,otagsn,otags,oymd)  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_pixel_out(st1,st2,outyears,nomdos,otagsn,otags,oymd)
!**********************************************************************!
integer, dimension(max_outputs) :: otagsn
integer :: nomdos,i,j,l,ii,oymd,outyears
character(len=str_len), dimension(max_outputs) :: otags
character(len=str_len) :: st1,st2,st3,st4,st5
!----------------------------------------------------------------------!

st1 = adjustl(st1)
ii = n_fields(st1)

if (ii>1) then
  call stripbs(st1,st2)
  st1 = adjustl(st1)
  st3 = 'MONTHLY'
  st4 = 'DAILY'
  st5 = 'ALL'
  if (stcmp(st1,st3)==1) then
    oymd = 1
  elseif (stcmp(st1,st4)==1) then
    oymd = 2
  elseif (stcmp(st1,st5)==1) then
    oymd = 3
  else
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) & 
 'The second field of line 15 of the input file must read MONTHLY DAILY or ALL.'
    write(*,'('' "'',A,''"'')') trim(st1)
    write(*,*) 'Output variable options:'
    write(*,'(1x,20a4)') (otags(j),j=1,15)
    write(*,'(1x,20a4)') (otags(j),j=16,nomdos)
    stop
  endif
  call stripbs(st1,st2)

  ii = n_fields(st1)
  if (ii>=1) then
    call STRIPBN(st1,outyears)
    if (outyears==-9999) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) &
 'The third field of the ''PIXEL'' line must contain an integer, if it exists.'
      stop
    endif
  elseif (ii==0) then
    outyears = 0
  endif
else
  outyears = 0
  oymd = 0
endif

do i=1,50
  otagsn(i) = 0
enddo

if (ii>1) then
  st1 = adjustl(st1)
  do i=1,ii-1
    call stripbs(st1,st2)
    l = ntags(otags,st2)
    if (l/=-9999) then
      otagsn(l) = 1
    else
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) 'Error in tag name on the ''PIXEL'' line.'
      write(*,'('' "'',A,''"'')') trim(st2)
      write(*,*) 'Available tags:'
      write(*,'(1x,20a4)') (otags(j),j=1,15)
      write(*,'(1x,20a4)') (otags(j),j=16,nomdos)
      stop
    endif 
 enddo
else
  do i=1,nomdos
    otagsn(i) = 1
  enddo
endif

end subroutine set_pixel_out





!**********************************************************************!
!                                                                      !
!                      set_subpixel_out :: output_methods              !
!                      ----------------------------------              !
!                                                                      !
! subroutine set_subpixel_out(st1,st2,outyears,nomdos, &               !
! otagsnft,otags,oymdft,out_cov,out_bio,out_bud,out_sen)               !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_subpixel_out(st1,st2,outyears,nomdos, &
 otagsnft,otags,oymdft,out_cov,out_bio,out_bud,out_sen)
!**********************************************************************!
integer,dimension(max_outputs) :: otagsnft
integer :: nomdos,i,j,l,ii,oymdft,outyears
character(len=str_len), dimension(max_outputs) :: otags
character(len=str_len) :: st1,st2,st3,st4,st5
logical :: out_cov,out_bio,out_bud,out_sen
!----------------------------------------------------------------------!

out_cov = .false.
out_bio = .false.
out_bud = .false.
out_sen = .false.

st1 = adjustl(st1)
ii = n_fields(st1)
if (ii>1) then
  call stripbs(st1,st2)
  st1 = adjustl(st1)
  st3 = 'MONTHLY'
  st4 = 'DAILY'
  st5 = 'ALL'
  if (stcmp(st1,st3)==1) then
    oymdft = 1
  elseif (stcmp(st1,st4)==1) then
    oymdft = 2
  elseif (stcmp(st1,st5)==1) then
    oymdft = 3
  else
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) &
 'The second field of line 16 of the input file must read MONTHLY DAILY or ALL.'
    write(*,'('' "'',A,''"'')') trim(st1)
    write(*,*) 'Output variable options:'
    write(*,'(1x,20a4)') (otags(j),j=1,15)
    write(*,'(1x,20a4)') (otags(j),j=16,nomdos),'cov ','bio ','bud ','sen '
    stop
  endif
  call stripbs(st1,st2)

  ii = n_fields(st1)
  if (ii>=1) then
    call STRIPBN(st1,outyears)
    if (outyears==-9999) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) & 
 'The third field of the ''SUB_PIXEL'' line must contain an integer, if it exists.'
      stop
    endif
  elseif (ii==0) then
    outyears = 0
  endif
else
  outyears = 0
  oymdft = 0
endif

do i=1,50
  otagsnft(i) = 0
enddo

if (ii>1) then
  st1 = adjustl(st1)
  do i=1,ii-1
    call stripbs(st1,st2)
    st3 = 'cov'
    if (stcmp(st2,st3)==1) then
      out_cov = .true.
    else
      st3 = 'bio'
      if (stcmp(st2,st3)==1) then
        out_bio = .true.
      else
        st3 = 'bud'
        if (stcmp(st2,st3)==1) then
          out_bud = .true.
        else
          st3 = 'sen'
          if (stcmp(st2,st3)==1) then
            out_sen = .true.
          else
            l = ntags(otags,st2)
            if (l/=-9999) then
              otagsnft(l) = 1
            else
              write(*,'('' PROGRAM TERMINATED'')')
              write(*,*) 'Error in tag name in the ''SUBPIXEL'' line.'
              write(*,'('' "'',A,''"'')') trim(st2)
              write(*,*) 'Available tag names:'
              write(*,'(1x,20a4)') (otags(j),j=1,15)
              write(*,'(1x,20a4)') (otags(j),j=16,nomdos),'cov ','bio ','bud ','sen '
              stop
            endif
          endif
        endif
      endif
    endif
  enddo
else
  do i=1,nomdos
    otagsnft(i) = 1
  enddo
  out_cov = .true.
  out_bio = .true.
  out_bud = .true.
  out_sen = .true.
endif

end subroutine set_subpixel_out



end module output_methods

