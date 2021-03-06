module veg_dynamics

use real_precision
use dims
use system_state
use site_parameters
use state_methods
use func
use tuning_parameters
use file_class
use file_object

implicit none

contains

!**********************************************************************!
!                                                                      !
!                     cover :: veg_dynamics                            !
!                      --- ---------------                             !
!                                                                      !
! subroutine COVER(nft,tmp,prc,firec,fireres,fprob,ftprop,             !
! check_closure)                                                       !
!                                                                      !                                                            
!----------------------------------------------------------------------!
!> @brief cover
!! @details ! First it calculates empirical fire probability fprob by  !
!! calling the FIRE subroutine. NEWGROWTH will kill cohorts that       !
!! exceeded the permissible age for the ft and do similar for cohorts  !
!! with low NPP and fire.It will also set ssv%cover to zero for those  !
!! cases which means cover will be available for use.It then           !
!! proportianally assigns new cover to the fts that need to make the   !
!! biggest jump in cover from the previous year to the current.        !
!! ftprop(ft) now describes the new cover to be ADDED to each ft.      !
!!                                                                     !
!! @author Mark Lomas                                                  ! 
!! @date Feb 2006                                                      !
!----------------------------------------------------------------------!
subroutine cover(nft,tmp,prc,firec,fireres,fprob,ftprop,check_closure)
!**********************************************************************!
real(dp) :: npp(max_cohorts),nps(max_cohorts),tmp(12,31),prc(12,31)
real(dp) :: firec,fprob,ftprop(max_cohorts),norm
real(dp) ::ftprop0(max_cohorts),total_carbon,old_total_carbon,mtmp,mprc
integer nft,ft,fireres
logical check_closure

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------!

do ft=1,ssp%cohorts
  npp(ft) = ssv(ft)%npp
  nps(ft) = ssv(ft)%nps
enddo

!----------------------------------------------------------------------!
! Compute the likelihood of fire in the current year 'fprob'.          !
!----------------------------------------------------------------------!
call fire(prc,tmp,fprob)

!----------------------------------------------------------------------!
! Take off cohorts past their age and fractions burnt by fire together !
! with low npp.Save the carbon and nitrogen for redistribution.        !
!----------------------------------------------------------------------!
call newgrowth(fprob,npp,nps,fireres,firec)

! ftprop(ft) here holds the % of cover for each pft I want for this
! year

do ft=1,nft
  ftprop(ft) = ftprop(ft)/100.0
  ftprop0(ft) = 0.0
enddo

! Calculates total current cover for each type by summing cohorts
do ft=1,ssp%cohorts
  ftprop0(pft(ft)%itag)=ftprop0(pft(ft)%itag) + ssv(ft)%cov
enddo

! Calculates ftprop(ft) which here stands for the cover it needs to ADD
! to each ft.
norm = 0.0d0
do ft=1,nft           
  if (ftprop(ft)>0.0d0) then
    ! Subtracts what I want for the new year from what I have.If its
    ! positive it shows how much I want to add.norm adds up all the
    ! additional cover I need for all covers
    ftprop(ft) = ftprop(ft) - ftprop0(ft)
    if (ftprop(ft)<0.0d0) ftprop(ft)=0.0d0 !can't remove cov
    norm = norm + ftprop(ft)
  endif
enddo

! It normalizes each cover I want to add for each ft by the total
! cover I want to add and multiplies by the available new cover
! ssp%new_cov which originates from cohorts that died or burned.
! Each ft that wants to increase will get a fraction of the new
! cover available.ftprop now stands for a fraction of ssp%new_cov
! and should sum up to ssp%new_cov
if (ssp%new_cov > 0.0) then
  do ft=1,nft
    ftprop(ft) = ftprop(ft)/norm*ssp%new_cov
  enddo
endif

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  call sum_carbon(total_carbon,.false.)
  if (abs(total_carbon+firec-old_total_carbon) > 1.0e-6) then
    write(*,*) 'Breach of carbon closure in COVER:',total_carbon+firec-old_total_carbon,' g/m^2.'
  endif
endif

end subroutine COVER





!**********************************************************************!
!                                                                      !
!               INITIALISE_NEW_COHORTS :: veg_dynamics                 !
!                -----------------------------------                   !
!                                                                      !
!     subroutine INITIALISE_NEW_COHORTS(nft,ftprop,check_closure)      !
!                                                                      !                                                            
!----------------------------------------------------------------------!
!> @brief INITIALISE_NEW_COHORTS
!! @details ! It gets ftprop(ft) which holds the cov it needs to add 
!! for each pft for this year from the COVER subroutine.One new cohort 
!! for each ft will carry this cover.It will be initialised relatively 
!! to its cover and depending on what was left in the pools from cohorts
!! that died.
!!
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine initialise_new_cohorts(nft,ftprop,check_closure)
!**********************************************************************!
real(dp) :: ftprop(max_cohorts),sumc,sumftprop,total_carbon,old_total_carbon
integer ft,nft,i,cohort,ierase
logical check_closure
real(dp) :: extraC,amount_taken,prop_taken,prop_needed,sum_prop_taken
real(dp), dimension(8) :: extras_c,extras_n
real(dp), dimension(3) :: extras_minn
real(dp), dimension(4) :: extras_soil_h2o
real(dp) :: extras_snow,extras_l_snow,extras_slc,extras_rlc,extras_sln
real(dp) :: extras_rln,sumslc
integer :: loop_check(max_cohorts),adjust_check

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Check that the soil isn't barren by containing too little soil       !
! carbon. If barren, then set new cover (ftprop) to bare ground for    !
! each pft.                                                            !
!----------------------------------------------------------------------!
do ft=2,nft
  sumc = 0.0
  do i=1,8
    sumc = sumc + ssp%xnew_c(ft,i)
  enddo
  if (ftprop(ft)>0.0) then
    if (sumc/ftprop(ft)<1000.0) then
      ftprop(1) = ftprop(1) + ftprop(ft)
      ftprop(ft) = 0.0
    endif
  endif
enddo

!----------------------------------------------------------------------!
! Set cover arrays for this years ft proportions, take carbon from     !
! litter to provide nppstore and canopy.                               !
!----------------------------------------------------------------------!
sumslc = 0.0
do ft=1,nft
  sumslc = sumslc + ssp%xnew_slc(ft)
enddo

!----------------------------------------------------------------------!
! Make the averaged extraC pool. This is made up from any pfts where   !
! the new proportions (ftprop) are less than the old ssp%new_cov.      !
!----------------------------------------------------------------------!
sumftprop = 0.0
do ft=1,nft
  sumftprop = sumftprop + ftprop(ft)
enddo

extras_c = 0.0
extras_n = 0.0
extras_minn = 0.0
extras_soil_h2o = 0.0
extras_snow = 0.0
extras_l_snow = 0.0
extras_slc = 0.0
extras_rlc = 0.0
extras_sln = 0.0
extras_rln = 0.0

sum_prop_taken = 0.0
prop_taken = 0.0

loop_check = 0
adjust_check = 0

! ssp%xnew_cov(ft) is the cover available for each pft due to age,fire
! or low npp removal and is ftprop(ft) the cover I want to add for this
! year for the pft.sumftprop is the total cover that needs to be added for all pft.

! If ssp%xnew_cov(ft) is greater than the one I want to
! add for the pft, then their will be carbon,nitrogen and water from that
! pfty available.These are added up in the extra_* variables and are 
! subtracted from xnew_* which held how much was available from age,
! fire and low npp removal.

do ft=1,nft
  if (ssp%xnew_cov(ft)-ftprop(ft)>1e-6) then
    loop_check(ft) = 1
    adjust_check = 1
    prop_taken = (ssp%xnew_cov(ft) - ftprop(ft))/sumftprop
    sum_prop_taken = sum_prop_taken + prop_taken

    do i=1,8
      amount_taken = prop_taken*ssp%xnew_c(ft,i)
      ssp%xnew_c(ft,i) = ssp%xnew_c(ft,i) - amount_taken
      extras_c(i) = extras_c(i) + amount_taken
    enddo

    do i=1,8
      amount_taken = prop_taken*ssp%xnew_n(ft,i)
      ssp%xnew_n(ft,i) = ssp%xnew_n(ft,i) - amount_taken
      extras_n(i) = extras_n(i) + amount_taken
    enddo

    do i=1,3
      amount_taken = prop_taken*ssp%xnew_minn(ft,i)
      ssp%xnew_minn(ft,i) = ssp%xnew_minn(ft,i) - amount_taken
      extras_minn(i) = extras_minn(i) + amount_taken
    enddo

    do i=1,4
      amount_taken = prop_taken*ssp%xnew_soil_h2o(ft,i)
      ssp%xnew_soil_h2o(ft,i) = ssp%xnew_soil_h2o(ft,i) - amount_taken
      extras_soil_h2o(i) = extras_soil_h2o(i) + amount_taken
    enddo

    amount_taken = prop_taken*ssp%xnew_snow(ft)
    ssp%xnew_snow(ft) = ssp%xnew_snow(ft) - amount_taken
    extras_snow = extras_snow + amount_taken

    amount_taken = prop_taken*ssp%xnew_l_snow(ft)
    ssp%xnew_l_snow(ft) = ssp%xnew_l_snow(ft) - amount_taken
    extras_l_snow = extras_l_snow + amount_taken

    amount_taken = prop_taken*ssp%xnew_slc(ft)
    ssp%xnew_slc(ft) = ssp%xnew_slc(ft) - amount_taken
    extras_slc = extras_slc + amount_taken

    amount_taken = prop_taken*ssp%xnew_rlc(ft)
    ssp%xnew_rlc(ft) = ssp%xnew_rlc(ft) - amount_taken
    extras_rlc = extras_rlc + amount_taken

    amount_taken = prop_taken*ssp%xnew_sln(ft)
    ssp%xnew_sln(ft) = ssp%xnew_sln(ft) - amount_taken
    extras_sln = extras_sln + amount_taken

    amount_taken = prop_taken*ssp%xnew_rln(ft)
    ssp%xnew_rln(ft) = ssp%xnew_rln(ft) - amount_taken
    extras_rln = extras_rln + amount_taken

  endif
enddo

!----------------------------------------------------------------------!
! Spread the extraC pool amongst the pfts where the new proportions    !
! (ftprop) the averaged pool.                                          !
!----------------------------------------------------------------------!
!if (adjust_check==1) then
prop_needed = 0.0

! If for a pft the new available cover is smaller than the one I want to 
! add then it gets carbon,nitrogen and water pools from the ones in
! excess calculated above as extras_*
do ft=1,nft
  if (ssp%xnew_cov(ft)-ftprop(ft)<-1.e-6) then
!  if (loop_check(ft)==0) then
    prop_needed = (ftprop(ft) - ssp%xnew_cov(ft))/sumftprop/&
      sum_prop_taken

    do i=1,8
      ssp%xnew_c(ft,i) = ssp%xnew_c(ft,i) + prop_needed*extras_c(i)
      ssp%xnew_n(ft,i) = ssp%xnew_n(ft,i) + prop_needed*extras_n(i)
    enddo
    do i=1,3
      ssp%xnew_minn(ft,i) = ssp%xnew_minn(ft,i) + prop_needed*extras_minn(i)
    enddo
    do i=1,4
      ssp%xnew_soil_h2o(ft,i) = ssp%xnew_soil_h2o(ft,i) + prop_needed*extras_soil_h2o(i)
    enddo

    ssp%xnew_snow(ft) = ssp%xnew_snow(ft) + prop_needed*extras_snow
    ssp%xnew_l_snow(ft) = ssp%xnew_l_snow(ft) + prop_needed*extras_l_snow
    ssp%xnew_slc(ft) = ssp%xnew_slc(ft) + prop_needed*extras_slc
    ssp%xnew_rlc(ft) = ssp%xnew_rlc(ft) + prop_needed*extras_rlc
    ssp%xnew_sln(ft) = ssp%xnew_sln(ft) + prop_needed*extras_sln
    ssp%xnew_rln(ft) = ssp%xnew_rln(ft) + prop_needed*extras_rln

  endif
enddo
!endif

sumslc = 0.0
do ft=1,nft
  sumslc = sumslc + ssp%xnew_slc(ft)
enddo

! Add a cohort for all pfts with cover>0
cohort = ssp%cohorts
do ft=1,nft
  if (ftprop(ft)>0.0) then
!    ierase=1
!    IF (pft_tab(ft)%phen==3.and.ssp%co2ftmap(ft,1)==1) THEN
!      IF (ssv(ssp%co2ftmap(ft,2))%sown==1) THEN
!        ierase=0
!        pft(ssp%co2ftmap(ft,2))=pft_tab(ft)
!      ENDIF
!    ENDIF
!    IF (ierase==1) THEN
    cohort = cohort + 1
    ssp%co2ftmap(ft,1) = ssp%co2ftmap(ft,1) + 1
    ssp%co2ftmap(ft,ssp%co2ftmap(ft,1)+1) = cohort

! Set plant functional type parameterisation for the new cohorts.
    pft(cohort) = pft_tab(ft)

! Set initial values of the system state.
    call initialise_state_cohort(cohort)

    ssv(cohort)%stemfr = ssv(1)%stemfr

    ssv(cohort)%nppstore(1) = pft(cohort)%stemx
    ssv(cohort)%nppstore(2) = pft(cohort)%stemx
    ssv(cohort)%nppstore(3) = pft(cohort)%stemx
    if ((pft(cohort)%mort < 5).and.(pft(cohort)%mort > 0)) then
      if (ssp%co2ftmap(ft,1) > 1) then
        ssv(cohort)%nppstore(1) = ssv(ssp%co2ftmap(ft,2))%nppstore(1)
        ssv(cohort)%nppstore(2) = ssv(ssp%co2ftmap(ft,2))%nppstore(2)
        ssv(cohort)%nppstore(3) = ssv(ssp%co2ftmap(ft,2))%nppstore(3)
      endif
    endif

! Set new cohort cover
    ssv(cohort)%cov = ftprop(ft)
    ssv(cohort)%ppm = pft_tab(ft)%ppm0
    ssv(cohort)%hgt = 0.004
    ssv(cohort)%age = 0.5
    ssv(cohort)%age = 1.0
    call set_new_soil_res(cohort,ftprop(ft)/sumftprop)
!    call set_new_soil_res2(cohort,ft,ftprop(ft)/sumftprop)

!----------------------------------------------------------------------!
! Take carbon from litter to balance the storage.
!----------------------------------------------------------------------!
    ssv(cohort)%slc = ssv(cohort)%slc - ssv(cohort)%nppstore(1)*ssv(cohort)%cov

    if (ssv(cohort)%slc<0.0) then
!----------------------------------------------------------------------!
! Make up any shortfall with soil carbon.
!----------------------------------------------------------------------!
      sumc = 0.0
      do i=1,8
         sumc = sumc + ssv(cohort)%c(i)
      enddo
      do i=1,8
        ssv(cohort)%c(i) = ssv(cohort)%c(i)*(1.0+ssv(cohort)%slc/sumc/ssv(cohort)%cov)
      enddo
      ssv(cohort)%slc = 0.0
    endif
!----------------------------------------------------------------------!
!ENDIF
  endif
enddo

ssp%cohorts = cohort

! This resets the ssp%xnew_* since won't be needed any more down the
! annual loop.
call reset_soil_res()

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  call sum_carbon(total_carbon,.false.)
  if (abs(total_carbon-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon closure in INITIALISE_NEW_COHORTS:', &
 total_carbon-old_total_carbon,' g/m^2.'
  endif
endif

end subroutine initialise_new_cohorts




!**********************************************************************!
!                                                                      !
!                          SUBROUTINE GROWTH                           ! 
!                          *****************                           !
!                                                                      !
!----------------------------------------------------------------------!
subroutine GROWTH(nft,lai,stembio,rootbio,check_closure)
!**********************************************************************!
real(dp) :: npp(max_cohorts),lai(max_cohorts),nps(max_cohorts),npr(max_cohorts), &
 evp(max_cohorts),rootbio,slc(max_cohorts),rlc(max_cohorts), &
 sln(max_cohorts),rln(max_cohorts),stembio, &
 total_carbon,old_total_carbon,ans
integer nft,ftmor(max_cohorts),ft,i
logical check_closure

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------!

do ft=1,ssp%cohorts
  npp(ft) = ssv(ft)%npp
  nps(ft) = ssv(ft)%nps
  npr(ft) = ssv(ft)%npr
  evp(ft) = ssv(ft)%evp
enddo

!----------------------------------------------------------------------!
! Initialise litter arrays, and add on leaf litter computed in DOLY.   !
!----------------------------------------------------------------------!
do ft=1,ssp%cohorts
  slc(ft) = 0.0
  rlc(ft) = 0.0
  sln(ft) = 0.0
  rln(ft) = 0.0
enddo

!----------------------------------------------------------------------!
! Thin vegetation where npp is not sufficient to maintain sensible     !
! growth rate.                                                         !
!----------------------------------------------------------------------!
call THIN(nft,npp,lai,nps,evp,slc,check_closure)

!----------------------------------------------------------------------!
!Compute leaf root and stem biomasses.                                 !
!----------------------------------------------------------------------!
stembio = 0.0
rootbio = 0.0
do ft=1,ssp%cohorts
  stembio = stembio + ssv(ft)%bio(1)*ssv(ft)%cov
  rootbio = rootbio + ssv(ft)%bio(2)*ssv(ft)%cov
enddo

do ft=1,ssp%cohorts
  ssv(ft)%npp = npp(ft)
  ssv(ft)%nps = nps(ft)
  ssv(ft)%npr = npr(ft)
  ssv(ft)%evp = evp(ft)
enddo

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  call SUM_CARBON(total_carbon,.false.)
  if (abs(total_carbon-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon closure in GROWTH:',total_carbon-old_total_carbon,' g/m^2.'
  endif
endif

end subroutine GROWTH





!**********************************************************************!
!                                                                      !
!                          SUBROUTINE THIN                             !
!                          ***************                             !
!                                                                      !
!----------------------------------------------------------------------!
subroutine THIN(nft,npp,lai,nps,evp,slc,check_closure)
!**********************************************************************!
real(dp) :: ftmat(max_cohorts),npp(max_cohorts),lai(max_cohorts)
real(dp) :: nps(max_cohorts),evp(max_cohorts),storelit(max_cohorts)
real(dp) :: slc(max_cohorts),pbionew,pbioold,no,pbio,ftcov(max_cohorts)
real(dp) :: covnew(max_age,max_cohorts),ppmnew(max_age,max_cohorts),shv
real(dp) :: lmv,nv,pi,hwv(max_age),emv,totno,totcov,scale(max_age)
real(dp) :: oldbio(max_cohorts,2),dimold,g0,gf,gm,grate(max_age)
real(dp) :: hgtnew(max_age),dimnew,maxhgt,minhgt,hc1,hc2,tcov1,tcov2,sum
real(dp) :: xxx,total_carbon,old_total_carbon
integer nft,ft,i,year,coh,ift
logical check_closure

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------!

pi = 3.1415926

hc1 = 0.05
hc2 = 2.0

do ft=1,ssp%cohorts
  do i=1,2
    oldbio(ft,i) = ssv(ft)%bio(i)
  enddo
  storelit(ft) = 0.0
enddo

do ft=1,ssp%cohorts
  ftmat(ft) = real(pft(ft)%mort)
enddo

!----------------------------------------------------------------------!
! FT loop for thinning and height competition.                         !
!----------------------------------------------------------------------!
do ft=1,nft
  if ((pft_tab(ft)%gr0>0.0).and.(ssp%co2ftmap(ft,1)>0)) then

!----------------------------------------------------------------------!
! Height competition.                                                  !
!----------------------------------------------------------------------!
    minhgt = 10000.0
    maxhgt =-10000.0
    i = 0
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      year = int(ssv(ift)%age+.50001)
      if (ssv(ift)%hgt>0) i = i + 1
      scale(year) = 1.0 - ssv(ift)%cov
      if (scale(year)<minhgt)  minhgt = scale(year)
      if (scale(year)>maxhgt)  maxhgt = scale(year)
    enddo

    if (maxhgt-minhgt<-1000.0) then
      do year=1,pft_tab(ft)%mort
        scale(year) = (scale(year) - minhgt)/(maxhgt - minhgt)
        scale(year) = (scale(year)*hc1 + 1.0 - hc1)**hc2
      enddo

      tcov1 = 0.0
      tcov2 = 0.0
      do coh=2,ssp%co2ftmap(ft,1)+1
        ift = ssp%co2ftmap(ft,coh)
        year = int(ssv(ift)%age+0.5)
        covnew(year,ift) = ssv(ift)%cov*scale(year)
        tcov1 = tcov1 + ssv(ift)%cov
        tcov2 = tcov2 + covnew(year,ift)
      enddo

      do coh=2,ssp%co2ftmap(ft,1)+1
        ift = ssp%co2ftmap(ft,coh)
        year = int(ssv(ift)%age+.5)
        if (covnew(year,ift)>0.0) then
          covnew(year,ift) = covnew(year,ift)*tcov1/tcov2
          do i=1,2
            ssv(ift)%bio(i) = ssv(ift)%bio(i)*ssv(ift)%cov/covnew(year,ift)
            oldbio(ift,i) = oldbio(ift,i)*ssv(ift)%cov/covnew(year,ift)
          enddo
          ssv(ift)%ppm = ssv(ift)%ppm*ssv(ift)%cov/covnew(year,ift)
          ssv(ift)%cov = covnew(year,ift)
        endif
      enddo

    endif
!----------------------------------------------------------------------!

    g0 = pft_tab(ft)%gr0
    gf = pft_tab(ft)%grf
    gm = real(pft_tab(ft)%mort)/10.0
    do year=1,pft_tab(ft)%mort
      grate(year) = (gf - g0)/gm*real(year - 1) + g0
      if (real(year)>=gm)  grate(year) = gf
    enddo

    totno = 0.0
    totcov = 0.0
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      totno = totno + ssv(ift)%cov*ssv(ift)%ppm
      totcov = totcov + ssv(ift)%cov
    enddo
    totno = totno/totcov

!----------------------------------------------------------------------!
! Compute cover and ppm to sustain a minimum growth rate, put these    !
! values in covnew and ppmnew.                                         !
!----------------------------------------------------------------------!
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)

      year = int(ssv(ift)%age+.5)

      shv = ssv(ift)%stem%tot(1)
      emv = evp(ift)/3600.0/1000.0
      lmv = lai(ift)*1.3/(lai(ift) + 3.0)
      lmv = 1.0
      nv = 1.0

      if (ssv(ift)%ppm*ssv(ift)%cov>0.0) then
!----------------------------------------------------------------------!
! Calculate the increase in diameter produced by stem NPP 'nps'.       !
!----------------------------------------------------------------------!

        if (emv>0.0) then
! hwv = theoretical maximum height (hydrolics)
          hwv(year) = pft(ift)%pdif*lmv/nv*(pft(ift)%xyl*shv/emv/pft(ift)%wden/10000.0)**0.5
        else
          hwv(year) = ssv(ift)%hgt*0.9
        endif

!----------------------------------------------------------------------!
! Calculate new height.
!----------------------------------------------------------------------!
        if (hwv(year)>0.0) then
          hgtnew(year) = (hwv(year)-ssv(ift)%hgt)/hwv(year)*0.5
        else
          hgtnew(year) = 0.0
        endif
        if (hgtnew(year)<0.0)  hgtnew(year) = 0.0
        hgtnew(year) = hgtnew(year) + ssv(ift)%hgt

!----------------------------------------------------------------------!
! pbio = g/individual
!----------------------------------------------------------------------!
        pbioold = ssv(ift)%bio(1)/ssv(ift)%ppm
        pbionew = (ssv(ift)%bio(1) + ssv(ift)%stem%tot(1)*(1.0 - &
 stlit(real(year,dp),ftmat(ift))))/ssv(ift)%ppm

!----------------------------------------------------------------------!
! Old diameter
!----------------------------------------------------------------------!
        if (ssv(ift)%hgt>0.0) then
! dimold = 2.0*(pbioold/1000000.0/ssv(ift)%hgt/
          dimold = 2.0*(pbioold/1000000.0/hgtnew(year)/pi/pft(ift)%wden)**0.5
        else
          dimold = 0.0
        endif

!----------------------------------------------------------------------!
! New diameter
!----------------------------------------------------------------------!
        if (hgtnew(year)>0.0) then
          dimnew = 2.0*(pbionew/1000000.0/hgtnew(year)/pi/pft(ift)%wden)**0.5
        else
          dimnew = 0.0
        endif

        if ((dimnew-dimold)/2.0>grate(year)) then
!----------------------------------------------------------------------!
! No thinning required.                                                !
!----------------------------------------------------------------------!
          ppmnew(year,ift) = ssv(ift)%ppm
          covnew(year,ift) = ssv(ift)%cov
        else
!----------------------------------------------------------------------!
! Thinning required.                                                   !
!----------------------------------------------------------------------!
          xxx = (grate(year)+dimold/2.0)**2.0*1000000.0*hgtnew(year)*pi*pft(ift)%wden
          ppmnew(year,ift) = (ssv(ift)%bio(1) + ssv(ift)%stem%tot(1)* &
 (1.0 - stlit(real(year,dp),ftmat(ift)))/100.0)/xxx
          covnew(year,ift) = ssv(ift)%cov

        endif
      else
        pbio = 0.0
        scale(year) = 0.0
        ppmnew(year,ift) = 0.0
        covnew(year,ift) = 0.0
        hgtnew(year) = 0.0
      endif
    enddo

!----------------------------------------------------------------------!
! Correct biomass array, and adjust litter for any thinned trees.      !
!----------------------------------------------------------------------!
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      year = int(ssv(ift)%age+0.5)
      if ((ssv(ift)%ppm>0.0).and.(ssv(ift)%cov>0.0)) then
        no = ssv(ift)%ppm*ssv(ift)%cov - ppmnew(year,ift)*covnew(year,ift)

        do i=1,2
          pbio = ssv(ift)%bio(i)/ssv(ift)%ppm
          slc(ift) = slc(ift) + pbio*no

          if (ssv(ift)%cov>0.0) then
            ssv(ift)%bio(i) = (oldbio(ift,i)*ssv(ift)%cov - pbio*no)/covnew(year,ift)
            ssv(ift)%bio(i) = oldbio(ift,i)*ppmnew(year,ift)/ssv(ift)%ppm
          endif
        enddo
        storelit(ift) = storelit(ift) + ssv(ift)%nppstore(1)/ssv(ift)%ppm*no
        ssv(ift)%cov = covnew(year,ift)
        ssv(ift)%ppm = ppmnew(year,ift)
        ssv(ift)%hgt = hgtnew(year)
      else
        ssv(ift)%cov = 0.0
        ssv(ift)%ppm = 0.0
        ssv(ift)%hgt = 0.0
      endif

    enddo
!----------------------------------------------------------------------!
  else
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      ssv(ift)%ppm = 0.0
      ssv(ift)%hgt = 0.0
    enddo
  endif
!----------------------------------------------------------------------!
! End of ft loop.
!----------------------------------------------------------------------!
enddo

do ft=1,nft
  ftcov(ft) = 0.0
enddo

do ft=1,ssp%cohorts
  ftcov(pft(ft)%itag) = ftcov(pft(ft)%itag) + ssv(ft)%cov
enddo

!----------------------------------------------------------------------!
! Correct nppstore to account for thinning.                            !
!----------------------------------------------------------------------!
do ft=1,ssp%cohorts
  if ((pft(ft)%gr0>0.0).and.(ssv(ft)%nppstore(1)>0.0).and.(ssv(ft)%cov>0.0)) then
    slc(ft) = slc(ft) +  storelit(ft)
    ssv(ft)%nppstore(1) = (ssv(ft)%nppstore(1)*ssv(ft)%cov - storelit(ft))/ssv(ft)%cov
  endif
enddo

do ft=1,ssp%cohorts
  ssv(ft)%slc = slc(ft)
enddo

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  call SUM_CARBON(total_carbon,.false.)
  if (abs(total_carbon-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon closure in THINNING:',total_carbon-old_total_carbon,' g/m^2.'
    stop
  endif
endif

end subroutine THIN





!**********************************************************************!
!                                                                      !
!                             FUNCTION find                            !
!                             *************                            !
!                                                                      !
!----------------------------------------------------------------------!
function find(tmp,prc)
!**********************************************************************!
real(dp) :: find,tmp(12),prc(12)
real(dp) :: tmplim,prclim
integer i
!----------------------------------------------------------------------!

tmplim = -5.0
prclim = 50.0

find = 0.0
do i=1,12
  if (tmp(i)>tmplim) then
    if (prc(i)<prclim) then
      find = find + prc(i)/(12.0*prclim)
    else
      find = find + 1/12.0
    endif
  else
    find = find + 1/12.0
  endif
enddo

end function find





!**********************************************************************!
!                                                                      !
!                             SUBROUTINE c3c4                          !
!                             ***************                          !
!                                                                      !
!----------------------------------------------------------------------!
subroutine c3c42(ftprop,npp,range)
!**********************************************************************!
real(dp) :: ftprop(max_cohorts),npp(max_cohorts),range
real(dp) :: grass,nd

grass = ftprop(2)
nd = npp(2) - npp(3)
ftprop(2) = grass*nd/(2.0*range) + grass/2.0
if (ftprop(2)>grass)  ftprop(2) = grass
if (ftprop(2)<0.0)  ftprop(2) = 0.0
ftprop(3) = grass - ftprop(2)

end subroutine c3c42





!**********************************************************************!
!                                                                      !
!                             SUBROUTINE c3c4                          !
!                             ***************                          !
!                                                                      !
!----------------------------------------------------------------------!
subroutine c3c4(ftprop,c3old,c4old,npp,nps)
!**********************************************************************!
real(dp) :: ftprop(max_cohorts),npp(max_cohorts),nps(max_cohorts),c3old,c4old
real(dp) :: grass,c3p,c4p,adj

if (c3old+c4old>0.0) then
  c3p = c3old/(c3old + c4old)
  c4p = c4old/(c3old + c4old)
else
  c3p = 0.5
  c4p = 0.5
endif

if (npp(2)*nps(2)+npp(3)*nps(3)>0.0) then
  adj = npp(2)*nps(2)/(npp(2)*nps(2) + npp(3)*nps(3)) - 0.5
else
  adj = 0.0
endif

c3p = c3p + adj
c4p = c4p - adj

if (c3p<0.0) then
  c3p = 0.0
  c4p = 1.0
endif

if (c4p<0.0) then
  c4p = 0.0
  c3p = 1.0
endif

grass = ftprop(2) + ftprop(3)

ftprop(2) = grass*c3p
ftprop(3) = grass*c4p

end subroutine c3c4





!**********************************************************************!
!                                                                      !
!                             SUBROUTINE GRASSREC                      !
!                             *******************                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine GRASSREC(nft,ftprop,gold,x,nat_map)
!**********************************************************************!
real(dp) :: ftprop(max_cohorts),gold,x,ntcov,ftt,ftpropo(max_cohorts)
integer nft,ft,nat_map(8)
!----------------------------------------------------------------------!

ntcov = 0.0
ftt = 0.0
do ft=5,8
  ntcov = ntcov + ftprop(nat_map(ft))*ssp%new_cov/100.0
  ftt = ftt + ftprop(nat_map(ft))
enddo

if (ntcov>gold*x) then
  do ft=4,nft
    ftpropo(ft) = ftprop(ft)
    ftprop(ft) = 100.0*gold*x*ftprop(ft)/(ftt*ssp%new_cov)
    ftprop(2) = ftprop(2) + ftpropo(ft) - ftprop(ft)
  enddo

  ntcov = 0.0
  do ft=4,nft
    ntcov = ntcov + ftprop(ft)*ssp%new_cov/100.0
  enddo
  if (abs(ntcov-gold*x)>0.000001) &
 write(fun%get_id('diag.dat'),'(''Treerec subroutine error'')')
endif

end subroutine GRASSREC





!**********************************************************************!
!                                                                      !
!                             SUBROUTINE BAREREC                       !
!                             ******************                       !
!                                                                      !
!----------------------------------------------------------------------!
subroutine BAREREC(ftprop,bpaold,x)
!**********************************************************************!
real(dp) :: ftprop(max_cohorts),bpaold,x
real(dp) :: nbp,cgcov
!----------------------------------------------------------------------!

! Current growth 'cgcov'
cgcov = 1.0 - ssp%new_cov

nbp = ssp%new_cov*ftprop(1)/100.0

if (nbp<bpaold*(1.0-x)) then
  ssp%new_cov = 1.0 - bpaold*(1.0 - x) - cgcov
  ssv(1)%cov = bpaold*(1.0 - x)
else
  ssv(1)%cov = ssp%new_cov*ftprop(1)/100.0
  ssp%new_cov = ssp%new_cov*(1.0 - ftprop(1)/100.0)
endif

end subroutine BAREREC





!**********************************************************************!
!                                                                      !
!                            SUBROUTINE ADDBIO                         !
!                            *****************                         !
! Add on biomass to bio array.                                         !
!----------------------------------------------------------------------!
subroutine ADDBIO(npp,nps,npr)
!**********************************************************************!
real(dp) :: npp(max_cohorts),nps(max_cohorts),npr(max_cohorts),npps,nppr
integer ft
!----------------------------------------------------------------------!

do ft=1,ssp%cohorts
  npps = nps(ft)*npp(ft)/100.0
  nppr = npr(ft)*npp(ft)/100.0
  ssv(ft)%bio(1) = ssv(ft)%bio(1) + npps
  ssv(ft)%bio(2) = ssv(ft)%bio(2) + nppr
enddo

end subroutine ADDBIO





!**********************************************************************!
!                                                                      !
!                            SUBROUTINE MKLIT                          !
!                            ****************                          !
! Make litter.                                                         !
!----------------------------------------------------------------------!
subroutine MKLIT(ftmat,slc,rlc,sln,rln,npp,nps,npr)
!**********************************************************************!
real(dp) :: ftmat(max_cohorts),slc(max_cohorts),rlc(max_cohorts),sln(max_cohorts), &
 rln(max_cohorts),npp(max_cohorts),nps(max_cohorts),npr(max_cohorts),npps,nppr,sl,rl
integer ft
!----------------------------------------------------------------------!

do ft=1,ssp%cohorts
  if ((pft(ft)%gr0>0.0).and.(ssv(ft)%age>1)) then
    npps = nps(ft)*npp(ft)/100.0
    nppr = npr(ft)*npp(ft)/100.0
    sl = stlit(ssv(ft)%age,ftmat(ft))
    rl = 0.8
    slc(ft) = slc(ft) + sl*npps*ssv(ft)%cov
    rlc(ft) = rlc(ft) + rl*nppr*ssv(ft)%cov
    ssv(ft)%bio(1) = ssv(ft)%bio(1) - sl*npps
    ssv(ft)%bio(2) = ssv(ft)%bio(2) - rl*nppr
  endif
enddo

do ft=1,ssp%cohorts
  sln(ft) = 0.0
  rln(ft) = 0.0
enddo

end subroutine MKLIT





!**********************************************************************!
function stlit(age,mat)
!----------------------------------------------------------------------!
real(dp) :: stlit,mat,temp,age
!----------------------------------------------------------------------!

stlit = 0.9*real(age)/mat + 0.1
if (stlit>1.0)  stlit = 1.0

if (age<mat) then
  stlit = 0.0
else
  stlit = 1.0
endif

temp = age/mat
if (temp>1.0) temp = 1.0

stlit = temp**0.5
stlit = stlit*0.6
!stlit = 0.0

end function stlit





!**********************************************************************!
!                                                                      !
!                            SUBROUTINE VEGMAT                         !
!                            *****************                         !
! Compute the years to reach maturity for each of the fts.             !
!----------------------------------------------------------------------!
subroutine VEGMAT(nft,ftmor,ftwd,ftxyl,ftpd,evp,lai,npp,nps,ftmat)
!**********************************************************************!
real(dp) :: ftwd(max_cohorts),ftxyl(max_cohorts),ftpd(max_cohorts),lai(max_cohorts)
real(dp) :: npp(max_cohorts),nps(max_cohorts),evp(max_cohorts),ftmat(max_cohorts),emxv,wd
real(dp) :: pd,lv,fs,shv,pi,lmvt,nvt,hwvt,dvt,massvt,minvt,ppvt
integer ft,nft,ftmor(max_cohorts)
!----------------------------------------------------------------------!

pi = 3.14159
do ft=2,nft
  if (npp(ft)>1.0e-6) then
    emxv = evp(ft)/3600.0/1000.0
    lv = lai(ft)
    fs = nps(ft)/100.0
    shv = fs*npp(ft)/100.0

    lmvt = lv*1.3/(lv + 3.0)
    nvt = 1.0 + 26.0*exp(-0.9*lv)
    nvt = 1.0

    wd = ftwd(ft)
    pd = ftpd(ft)

    hwvt = pd*lmvt*sqrt(ftxyl(ft)*shv/(emxv*wd*10000.0))/nvt
    dvt = 0.0028*hwvt**1.5
    massvt = pi*(dvt/2.0)**2.0*hwvt*wd*nvt
    minvt = pi*((dvt + 0.001/nvt)/2.0)**2.0*hwvt*wd*nvt - massvt
    ppvt = shv/minvt
    ftmat(ft) = massvt*ppvt/shv
    ftmat(ft) = 1.0
  else
    ftmat(ft) = 0.0
  endif
  if (ftmor(ft)<ftmat(ft)) ftmat(ft) = ftmor(ft)
enddo

end subroutine VEGMAT





!**********************************************************************!
!                                                                      !
!                            SUBROUTINE FIRE                           !
!                            ***************                           !
! Compute the fire return interval 'fri' and returns the probability   !
! of a fire in the current year 'fprob'.                               !
!----------------------------------------------------------------------!
subroutine FIRE(dprc,dtmp,fprob)
!**********************************************************************!
real(dp) :: fri,fprob,dtmp(12,31),dprc(12,31),prct(12),prco,lim1,lim2,maxfri
real(dp) :: pow,tlim1,tlim2,totp,indexx,tadj,weight,tmp(12),prc(12)
integer k,no,ind,i,j
!----------------------------------------------------------------------!

do i=1,12
  tmp(i) = 0.0
  prc(i) = 0.0
  do j=1,30
    tmp(i) = tmp(i) + dtmp(i,j)
    prc(i) = prc(i) + dprc(i,j)
  enddo
  tmp(i) = tmp(i)/30
enddo

!----------------------------------------------------------------------!
! Fire model parameters.                                               !
!----------------------------------------------------------------------!
no = 3
lim1 = 150.0
lim2 = 50.0
tlim1 =  0.0
tlim2 = -5.0
weight = 0.5
maxfri = 800.0
pow = 3.0

!----------------------------------------------------------------------!
! Adjust for temperature.                                              !
!----------------------------------------------------------------------!
do i=1,12
  if (tmp(i)>tlim1) then
    tadj = 0.0
  elseif (tmp(i)<tlim2) then
    tadj = lim1
  else
    tadj = (tlim1 - tmp(i))/(tlim1 - tlim2)*lim1
  endif
  prct(i) = min(lim1,prc(i) + tadj)
enddo

!----------------------------------------------------------------------!
! Calculate yearly componet of index.                                  !
!----------------------------------------------------------------------!
totp = 0.0
do k=1,12
  totp = totp + prct(k)/lim1/12.0
enddo

!----------------------------------------------------------------------!
! Calculate monthly component of index.                                !
!----------------------------------------------------------------------!
prco = 0.0
do k=1,no
  ind = minx(prct)
  prco = prco + min(lim2,prct(ind))/real(no)/lim2
  prct(ind) = lim2
enddo

!----------------------------------------------------------------------!
! Compute fire return interval, and convert to probability.            !
!----------------------------------------------------------------------!
indexx = weight*(prco) + (1.0 - weight)*totp

fri = indexx**pow*maxfri

if (fri<2.0)  fri = 2.0
fprob = (1.0 - exp(-1.0/fri))*tgp%p_fprob

end subroutine FIRE





!**********************************************************************!
!                                                                      !
!                     NEWGROWTH :: veg_dynamics                        !
!                      ----------------------                          !
!                                                                      !
!        subroutine NEWGROWTH(fprob,npp,nps,fireres,firec)             !
!                                                                      !                                                            
!----------------------------------------------------------------------!
!> @brief NEWGROWTH
!! @details ! Compute newgrowth and alter cover array accordingly.                 
!! It first checks the age of each cohort.If it has exceeded the limit
!! set by pft%mort,it kills it by calling the ACCUMUATE_DIST_SOIL_RES.
!! It then adds +1 to the age of all other cohorts.Calls COMPRESS_STATE
!! to remove dead/empty cohorts resulting from the previous step.
!! Performs similar process to kill cohorts due to low NPP and fire. 
!!
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine NEWGROWTH(fprob,npp,nps,fireres,firec)
!**********************************************************************!
real(dp) ::  fprob,npp(max_cohorts),nps(max_cohorts),tmor,tmor0,npp0,firec,xfprob
integer ft,fireres
!----------------------------------------------------------------------!
xfprob = fprob
firec = 0.0

do ft=1,max_pfts
  ssp%co2ftmap(ft,1) = 0
enddo

!----------------------------------------------------------------------!
! Take away veg that has died of old age ie > than pft%mort, and       !
! age the veg by one year.                                             !
!----------------------------------------------------------------------!
do ft=1,ssp%cohorts
!  IF (ssv(ft)%sown/=1) THEN 
    if (ssv(ft)%age > real(pft(ft)%mort)-0.1) then
      call ACCUMULATE_DIST_SOIL_RES(ft,1.0_dp)
    else
      ssv(ft)%age = ssv(ft)%age + 1.0
    endif
!  ENDIF
enddo

!----------------------------------------------------------------------!
! Remove dead cohorts and rearranges the remaining ones                !
!----------------------------------------------------------------------!
call COMPRESS_STATE()

!----------------------------------------------------------------------!
! Take away veg that is burnt or has died through a hard year ie small !
! LAI.                                                                 !
!----------------------------------------------------------------------!
do ft=1,ssp%cohorts

!----------------------------------------------------------------------!
! 'tmor' is the mortality rate of the forrest based on 'npp'.          !
!----------------------------------------------------------------------!
  npp0 = 0.2
  tmor0 = 6.0
  if (npp(ft)/100.0<0.1) then
    tmor = 1.0
  elseif (npp(ft)/100.0<3.0) then
    tmor = (0.04 - npp0)/3.0*(npp(ft)/100.0) + npp0
  elseif (npp(ft)/100.0<tmor0) then
    tmor = 0.04/(3.0 - tmor0)*(npp(ft)/100.) + &
 3.0*0.04/(tmor0 - 3.0) + 0.04
  else
    tmor = 0.0
  endif
  if (tmor<0.01)  tmor = 0.01

!IF (npp(ft)*nps(ft)/100.0>10.0) THEN
    tmor = 0.002
!ELSE
!  tmor = 1.0 - npp(ft)*nps(ft)/1000.0
!ENDIF
  if (ssv(ft)%age > real(pft(ft)%mort)/2.0) then
    tmor = (1.0 - 2.0*(pft(ft)%mort - ssv(ft)%age)/pft(ft)%mort)/10.0
  else
    tmor = 0.002
  endif
!  tmor = 0.0
  
!----------------------------------------------------------------------!
! kill off trees with no storage left.                                 !
!----------------------------------------------------------------------!
  if (.not.(ssv(ft)%nppstore(1))>0.1) then
    tmor = 1.0
  endif
  
  if (ssv(ft)%age<real(fireres)) then
    fprob = xfprob
  else
    fprob = 0.0
  endif
  if (fireres<0) fprob = real(-fireres)/1000.0

  call ACCUMULATE_DIST_SOIL_RES(ft,tmor)
  
  ! No fire for crop phenology  
  if (pft(ft)%phen==3) fprob=0.0
  call ACCUMULATE_BURNT_SOIL_RES(ft,fprob,firec)

enddo

call COMPRESS_STATE()

end subroutine NEWGROWTH



end module veg_dynamics

