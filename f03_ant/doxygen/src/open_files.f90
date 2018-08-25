module open_files

use real_precision
use func
use pft_parameters
use file_class
use file_object
use input_methods

implicit none

contains

!**********************************************************************!
!                                                                      !
!                     open_diag :: open_files                          !
!                     -----------------------                          !
!                                                                      !
! subroutine open_diag()                                               !
!----------------------------------------------------------------------!
!> @brief Open diagnostic file.
!! @details This is used to report any potential errors.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_diag()
!**********************************************************************!
integer :: fid
!----------------------------------------------------------------------!

call fun%open(trim(inp%dirs%output)//'/diag.dat',fid)
write(fid,'(''ERRORS'')')

end subroutine open_diag





!**********************************************************************!
!                                                                      !
!                     open_snapshots :: open_files                     !
!                     ----------------------------                     !
!                                                                      !
! subroutine open_snapshots(snp_no,snpshts)                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open snapshot files file.
!! @details These files contain the entire system state for the final
!! day of the years given in the input files. The files can be used
!! the make site 'cartoons' using the utility sorfware.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_snapshots(snp_no,snpshts)
!----------------------------------------------------------------------!
character(len=str_len) :: st1,st4
integer :: i,fno,snp_no
integer, dimension(max_years) :: snpshts
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Open snapshot files.                                                 !
!----------------------------------------------------------------------!
if (snp_no>0) then
  do i=1,snp_no
    fno = 100 + (i-1)*4
    st1=in2st(snpshts(i))
    call STRIPB(st4)
    open(fno+1,file=trim(inp%dirs%output)//'/initbio_'&
 //trim(st1)//'.dat')
    open(fno+2,file=trim(inp%dirs%output)//'/initcov_'&
 //trim(st1)//'.dat')
    open(fno+3,file=trim(inp%dirs%output)//'/initppm_'&
 //trim(st1)//'.dat')
    open(fno+4,file=trim(inp%dirs%output)//'/inithgt_'&
 //trim(st1)//'.dat')
  enddo
endif

end subroutine open_snapshots






!**********************************************************************!
!                                                                      !
!                     open_default :: open_files                       !
!                     --------------------------                       !
!                                                                      !
! subroutine open_default()                                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open default output files.
!! @details These files  will contain yearly values.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_default()
!**********************************************************************!
integer :: fid
!----------------------------------------------------------------------!

call fun%open(trim(inp%dirs%output)//'/lai.dat',fid)
call fun%open(trim(inp%dirs%output)//'/npp.dat',fid)
call fun%open(trim(inp%dirs%output)//'/scn.dat',fid)
call fun%open(trim(inp%dirs%output)//'/snn.dat',fid)
call fun%open(trim(inp%dirs%output)//'/nep.dat',fid)
call fun%open(trim(inp%dirs%output)//'/swc.dat',fid)
call fun%open(trim(inp%dirs%output)//'/biomass.dat',fid)
call fun%open(trim(inp%dirs%output)//'/bioind.dat',fid)
call fun%open(trim(inp%dirs%output)//'/covind.dat',fid)
call fun%open(trim(inp%dirs%output)//'/rof.dat',fid)
call fun%open(trim(inp%dirs%output)//'/fcn.dat',fid)
call fun%open(trim(inp%dirs%output)//'/nppstore.dat',fid)
call fun%open(trim(inp%dirs%output)//'/stembio.dat',fid)
call fun%open(trim(inp%dirs%output)//'/rootbio.dat',fid)
call fun%open(trim(inp%dirs%output)//'/leafper.dat',fid)
call fun%open(trim(inp%dirs%output)//'/stemper.dat',fid)
call fun%open(trim(inp%dirs%output)//'/rootper.dat',fid)
call fun%open(trim(inp%dirs%output)//'/sresp.dat',fid)
call fun%open(trim(inp%dirs%output)//'/evt.dat',fid)
call fun%open(trim(inp%dirs%output)//'/gpp.dat',fid)
call fun%open(trim(inp%dirs%output)//'/lch.dat',fid)
call fun%open(trim(inp%dirs%output)//'/prc.dat',fid)
call fun%open(trim(inp%dirs%output)//'/nbp.dat',fid)
call fun%open(trim(inp%dirs%output)//'/trn.dat',fid)
call fun%open(trim(inp%dirs%output)//'/fab.dat',fid)
call fun%open(trim(inp%dirs%output)//'/tmp.dat',fid)
call fun%open(trim(inp%dirs%output)//'/hum.dat',fid)

end subroutine open_default





!**********************************************************************!
!                                                                      !
!                     open_optional :: open_files                      !
!                     ----------------------------                     !
!                                                                      !
! subroutine open_optional(nft,out_cov,out_bio,out_bud,out_sen)        !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open optional output files.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_optional(out_cov,out_bio,out_bud,out_sen)
!**********************************************************************!
integer :: ft,nft,iofn,fid
logical :: out_cov,out_bio,out_bud,out_sen

!----------------------------------------------------------------------!
! Open optional yearly cover, biomass, budburst and senescence files.  !
!----------------------------------------------------------------------!
do ft=1,inp%npft
  if (out_cov) then
    call fun%open(trim(inp%dirs%output)//'/cov_'//trim(pft_tab(ft)%tag)//'.dat',fid)
  endif
  if (out_bio) then
    call fun%open(trim(inp%dirs%output)//'/bio_'//trim(pft_tab(ft)%tag)//'.dat',fid)
  endif
  if (out_bud) then
    call fun%open(trim(inp%dirs%output)//'/bud_'//trim(pft_tab(ft)%tag)//'.dat',fid)
  endif
  if (out_sen) then
    call fun%open(trim(inp%dirs%output)//'/sen_'//trim(pft_tab(ft)%tag)//'.dat',fid)
  endif
enddo

end subroutine open_optional



!**********************************************************************!
!                                                                      !
!                     open_tile_pft :: open_files                      !
!                     ---------------------------                      !
!                                                                      !
! subroutine open_tile_pft()                                           !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open input and/or output files for state vector.
!! @details Opens input and/or output files for state vector
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_tile_pft()
!**********************************************************************!
integer :: ft,nft,i,fid

!----------------------------------------------------------------------!
! Open tile and pft output files.                                      !
!----------------------------------------------------------------------!
! Yearly tile files
if (inp%output%tile_yearly>0) then
  do i=1,inp%output%tile_vars_yearly%n
call fun%open(trim(inp%dirs%output)//'/yearly_'//trim(inp%output%tile_vars_yearly%tag(i))//'.dat',fid)
  enddo
endif
! Monthly tile files
if (inp%output%tile_monthly>0) then
  do i=1,inp%output%tile_vars_monthly%n
call fun%open(trim(inp%dirs%output)//'/monthly_'//trim(inp%output%tile_vars_yearly%tag(i))//'.dat',fid)
  enddo
endif
! Daily tile files
if (inp%output%tile_daily>0) then
  do i=1,inp%output%tile_vars_daily%n
call fun%open(trim(inp%dirs%output)//'/daily_'//trim(inp%output%tile_vars_yearly%tag(i))//'.dat',fid)
  enddo
endif

! Yearly pft files
if (inp%output%pft_yearly>0) then
  do i=1,inp%output%pft_vars_yearly%n
    do ft=1,inp%npft
call fun%open(trim(inp%dirs%output)//'/yearly_'//trim(inp%output%pft_vars_yearly%tag(i))//'_'//trim(pft_tab(ft)%tag)//'.dat',fid)
    enddo
  enddo
endif
! Monthly pft files
if (inp%output%pft_monthly>0) then
  do i=1,inp%output%pft_vars_monthly%n
    do ft=1,inp%npft
call fun%open(trim(inp%dirs%output)//'/monthly_'//trim(inp%output%pft_vars_monthly%tag(i))//'_'//trim(pft_tab(ft)%tag)//'.dat',fid)
    enddo
  enddo
endif
! Daily pft files
if (inp%output%pft_daily>0) then
  do i=1,inp%output%pft_vars_daily%n
    do ft=1,inp%npft
call fun%open(trim(inp%dirs%output)//'/daily_'//trim(inp%output%pft_vars_daily%tag(i))//'_'//trim(pft_tab(ft)%tag)//'.dat',fid)
    enddo
  enddo
endif



end subroutine open_tile_pft




!**********************************************************************!
!                                                                      !
!                     open_state :: open_files                         !
!                     ------------------------                         !
!                                                                      !
! subroutine open_state()                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open input and/or output files for state vector.
!! @details Opens input and/or output files for state vector
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_state()
!**********************************************************************!
integer :: fid

!----------------------------------------------------------------------!
! Open input files for state vector.                                   !
!----------------------------------------------------------------------!
if (.not.inp%run%initialise_state) then
  call fun%open_old(trim(inp%dirs%input)//'/init.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/initbio.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/initcov.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/initppm.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/inithgt.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/initleaf.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/initstem.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/initroot.dat',fid)
  call fun%open_old(trim(inp%dirs%input)//'/initmisc.dat',fid)
endif

!----------------------------------------------------------------------!
! Open output files for state vector.                                  !
!----------------------------------------------------------------------!
call fun%open(trim(inp%dirs%output)//'/init.dat',fid)
if (inp%run%output_state) then
  call fun%open(trim(inp%dirs%output)//'/initbio.dat',fid)
  call fun%open(trim(inp%dirs%output)//'/initcov.dat',fid)
  call fun%open(trim(inp%dirs%output)//'/initppm.dat',fid)
  call fun%open(trim(inp%dirs%output)//'/inithgt.dat',fid)
  call fun%open(trim(inp%dirs%output)//'/initleaf.dat',fid)
  call fun%open(trim(inp%dirs%output)//'/initstem.dat',fid)
  call fun%open(trim(inp%dirs%output)//'/initroot.dat',fid)
  call fun%open(trim(inp%dirs%output)//'/initmisc.dat',fid)
endif

end subroutine open_state





!**********************************************************************!
!                                                                      !
!                     open_site_info :: open_files                     !
!                     ----------------------------                     !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Write lat & lon for output files.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_site_info()
!**********************************************************************!
integer :: fid,ft

call fun%open(trim(inp%dirs%output)//'/site_info.dat',fid)
write(fid, &
'(''COUNTRY           ID     lat      lon     co20  co2f    sand   &
 &silt   bulk     wp     fc    swc    dep       '' &
& ,50A13)') (pft_tab(ft)%tag,ft=1,inp%npft)

end subroutine open_site_info





!**********************************************************************!
!                                                                      !
!                     write_lat_lon :: open_files                      !
!                     ---------------------------                      !
!                                                                      !
!                                     !
! subroutine WRITE_LAT_LON(lat,lon,nft,out_cov,out_bio,out_bud,        !
! out_sen,oymdft,otagsn,oymd,otagsnft,outyears1,outyears2)             !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Write lat & lon for output files.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine write_lat_lon(lat,lon,nft,out_cov,out_bio,out_bud,out_sen, &
 oymdft,otagsn,oymd,otagsnft,outyears1,outyears2)
!**********************************************************************!
integer :: i,ft,nft,outyears1,outyears2,iofn,oymdft,oymd
integer, dimension(max_outputs) :: otagsn,otagsnft
real(dp):: lat,lon
logical :: out_cov,out_bio,out_bud,out_sen
character(len=str_len) :: st1

do i=1,fun%n
  if (fun%opened(i)) then
    call filename_base(fun%names(i),st1)
    if (st1(1:13)/='site_info.dat') then
      if ((st1(1:8)=='monthly_').or.(st1(1:6)=='daily_')) then
        write(fun%get_id_n(i),'(f7.3,f9.3)') lat,lon
      else
        write(fun%get_id_n(i),'(f7.3,f9.3)',advance='NO') lat,lon
      endif
    endif
  endif
enddo

end subroutine write_lat_lon



end module open_files



