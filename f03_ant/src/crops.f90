!> @brief Crop subroutines
!! @details
!! @author 
!! @date 

module crops

use real_precision
use pft_parameters
use site_parameters
use func
use light_methods
use system_state
use input_file

implicit none

contains

!**********************************************************************!
!                                                                      !
!                     seasonality :: crops                             !
!                     ---------------------                            !
!                                                                      !
! subroutine seasonality                                               !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Crops Seasonality
!! @details Following Waha et al 2012 (published version of LGJ 
!! van Bussel's thesis (2011)
!! chapter 5), calculate variation coefficients for seasonality, which may 
!! control sowing dates.  Then, get the sowing day based on precipitation or 
!! summer/winter temperature thresholds for each crop.
!
!! Following van Bussel et al 2015, Appendix S1 (van Bussel thesis 2011 fig 6.1a
!! on page 119), we also calculate the optimal photoperiod for each crop based on
!! longest and possibly shortest daylengths of the year.  For example, the number
!! of hours in the summer solstice is the optimum photoperiod for wheat, but for
!! maize it is sensitive to latitude: longest*(1-1/(longest-shortest)). Without
!! hardcoded crop types, we pass a parameter "pscale" (0 or 1) for each crop to 
!! tell us what to do: longest*(1-pscale/(longest-shortest)) 
!! 
!! Also following van Bussel et al 2015, Appendix S1 (van Bussel thesis 2011 
!! chapter 5), calculate vernalisation requirements based on temperature.
!!
!! REFERENCES:
!!
!! Waha et al 2012, Global Ecology and Biogeography 21:247-259
!! van Bussel et al 2015, Global Ecology and Biogeography (earlyview)
!! Wang and Engel 1998 (Agri Sys 58:1-24)
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!!
subroutine seasonality(tmp,prc,cld,thty_dys,nft,year,nn1)
!**********************************************************************!
real(dp) :: tmp(12,31),prc(12,31),cld(12)
integer  :: thty_dys,nft,year

real(dp) :: museas,seastmp,seasprc,hrs(12)
integer  :: iseas,mdoy(12),mmid(12),mnth,day,ft,k,daysiny
integer  :: msow,mcoldest,icoldest,nwarm,mcold,ntot,nvern,i,j,m
integer  :: nydays,slen(12),nn1
real(dp) :: fv,vdays,vegphu(12),repphu(12),totveg,wft,wfp,vrat
real(dp) :: totphu(12),targ,totrep
logical  :: dead
real(dp) :: qdir,qdiff,q(12),pet(12) ! internal, for getwet
!----------------------------------------------------------------------!!

  ! Days in a year
  daysiny=0
  do mnth=1,12
    daysiny=daysiny+int(no_days(year,mnth,thty_dys))
  endDO
  
  ! See Page 94 of van Bussel's thesis, section 2.2.1, eqns 5.1-5.3
  ! or page 249 of Waha et al 2012
  ! Calculates variation coefficient for temperature and precipitation
  museas = SUM(ssp%emnthtmp(:,nn1+1))/12.0d0+273.15
  seastmp = SQRT(1.0d0/11.0d0* &
    SUM((ssp%emnthtmp(:,nn1+1)+273.15-museas)*(ssp%emnthtmp(:,nn1+1)+273.15-museas)))/museas
  museas = SUM(ssp%emnthprc(:,nn1+1))/12.0d0
  seasprc = SQRT(1.0d0/11.0d0* &
  sum((ssp%emnthprc(:,nn1+1)-museas)*(ssp%emnthprc(:,nn1+1)-museas)))/museas
  
  ! Flag the seasonality type to use
  ssp%iseas = 0
  ! Thresholds for seasonality come from van Bussel's thesis page 95
  ! or from Waha et al 2012 page 249
  IF (seastmp>0.01d0.AND.seasprc>0.4d0) THEN
  ! Only use one seasonality type, depending on whether 
  ! minimum mean monthly temperature exceeds 10C
    IF (MINVAL(ssp%emnthtmp(:,nn1+1))>10.0d0) THEN
      ! assume seasons controlled primarily by precipitation
      ssp%iseas = 1
    ELSE
      ! assume seasons controlled primarily by temperature
      ssp%iseas = 2
    ENDIF ! min(ssp%emnthtmp)>10C
  ELSEIF (seastmp>0.01d0) THEN 
    ssp%iseas = 2 ! seasons controlled by temperature
  ELSEIF (seasprc>0.4d0) THEN
    ssp%iseas = 1 ! seasons controlled by precipitation
  ENDIF ! check for combined temp and precip seasonality
  
  
  ! For each month, get day-of-year for the first day and mid-month day length
  mdoy(1) = 1; mmid(1) = int(no_days(year,1,thty_dys)/2)
  hrs(1) = dayl(ssp%lat,15) ! Number of hours of daylight in the specified day
  call pfd(ssp%lat,15,hrs(1),cld(1),qdir,qdiff,q(1)) ! photon flux density
  DO mnth = 2,12
    mdoy(mnth) = mdoy(mnth-1)+no_days(year,mnth-1,thty_dys)
    mmid(mnth) = INT(no_days(year,mnth,thty_dys)/2)+mdoy(mnth)-1
    hrs(mnth) = dayl(ssp%lat,mmid(mnth))
    call pfd(ssp%lat,mmid(mnth),hrs(mnth),cld(mnth),qdir,qdiff,q(mnth))
  ENDDO
  
  ! Recalculate sowday each year
  pft_tab(1:nft)%sowday(nn1+1) = 0
  
  ! Follow van Bussel (section 2.2.2) or Waha et al p 249 to find sowing day
  IF (ssp%iseas==1) THEN ! If sowing depends on precipitation         
      ! Run four-month sums of precip/PET ratios 
      ! mnth is an output variable holding the first month of the wet season
      CALL getwet(hrs,q,mnth,pet,nn1,year,thty_dys)      
      ! Find the first julian day of the first wet-season month which is wet enough
      ! All crops will be sown on the same day in this case
      DO day=1,no_days(year,mnth,thty_dys)
          IF (pft_tab(3)%sowday(nn1+1)==0.and.prc(mnth,day)>0.1d0) &
          pft_tab(3:nft)%sowday(nn1+1)=day-1+mdoy(mnth)
      ENDDO ! day=1,no_days(year,mnth,thty_dys)    
  ELSEIF (ssp%iseas==2) THEN ! sowing depends on temperature
      pet(:)=1.0d0
      DO ft=3,nft ! there will be crop-specific threshold temperatures
          ! Figures whether the crop will have a spring sow date
          ! According to the pft parameters this holds for all summer crops
          ! (see van Bussel Table 5.1 on page 97, or Waha et al Table 1 page 250)
          IF (pft_tab(ft)%sowthresh(2)>pft_tab(ft)%lethal(1)) then
              pft_tab(ft)%sowday(nn1+1)=summerday(pft_tab(ft)%sowthresh(2),mmid,nn1)
          ENDIF ! sowday,avtmp>=sumtsow
          ! This avoids doing summer crops again
          IF (pft_tab(ft)%sowthresh(1)>=pft_tab(ft)%lethal(2)) CYCLE 
          ! Does crops that require vernilization
          ! (see van Bussel Table 5.1 on page 97, or Waha et al Table 1 page 250) 
          IF (pft_tab(ft)%sowthresh(1)>maxval(ssp%emnthtmp(:,nn1+1))) then  
              pft_tab(ft)%sowday(nn1+1)=winterday(maxval(ssp%emnthtmp(:,nn1+1))-0.1,mmid,nn1)
          ELSE
              pft_tab(ft)%sowday(nn1+1)=winterday(pft_tab(ft)%sowthresh(1),mmid,nn1)
          ENDIF
      ENDDO ! ft=1,nft
  ENDIF ! seasonality

  
  ! Estimate heat units,vernalisation days
  DO ft=3,nft
      IF(pft_tab(ft)%phen/=3) CYCLE

      ! Aseasonal climate; just sow on new years day as per van Bussel 2011 p 95 
      IF (pft_tab(ft)%sowday(nn1+1)==0) pft_tab(ft)%sowday(nn1+1)=1

      msow=1; mcoldest=1
      ! Get the sowing month and the coldest month
      DO mnth=1,12
          IF (pft_tab(ft)%sowday(nn1+1)>=mdoy(mnth)) msow=mnth
          IF (ssp%emnthtmp(mnth,nn1+1)<=ssp%emnthtmp(mcoldest,nn1+1)) mcoldest=mnth
      ENDDO ! mnth=1,12
    
    ! Get the midday of the coldest month
    icoldest=mmid(mcoldest);
    IF (ssp%iseas<2) icoldest=0
    
    ! croptype(2)=1 means requires vernalization
    ! croptype(2)=0 does not. So fv will attain a value before the loop
    ! of 0 for vern and 1 for non-vern.The latter means it will be ignored
    ! for non vern plants
    fv=1.0-pft_tab(ft)%croptype(2); vdays=0.0; dead=.FALSE.
    j=mdoy(msow)-1; vegphu(:)=0.0; repphu(:)=0.0
    nwarm=0; mcold=0; ntot=0; totveg=0.0; nvern=0

    ! Cycle through all 12 months starting with the sowing month
    ! k is the number of months added to the sowing month in each iteration
    do k=0,11 
        if (dead) cycle
        !Calendar month
        mnth=msow+k
        ! Starts again from January if it reaches December
        if (mnth>12) mnth=mnth-12
        dead=(ssp%emnthtmp(mnth,nn1+1)<=pft_tab(ft)%lethal(1).or. &
          ssp%emnthtmp(mnth,nn1+1)>=pft_tab(ft)%lethal(2))
        if (dead) CYCLE
        ! After growing season over; no more PHU
        if (mcold>0) CYCLE

        ! Calculate the mean daily heat units for this month (tmp*wft*wfp)
          
        ! wft is the mean temperature response for the month (0-1) and wfp
        ! is the mean photoperiod response (0-1)
        ! This is for the vegetative state for both type of plants
        CALL wangengel(pft_tab(ft)%cardinal(1),pft_tab(ft)%cardinal(2), &
        pft_tab(ft)%cardinal(3),ssp%emnthtmp(mnth,nn1+1),pft_tab(ft)%croptype(1), &
        pft_tab(ft)%photoperiod(1),pft_tab(ft)%photoperiod(2),hrs(mnth),wft,wfp)
        
        ! For crops that require vernalization.If the crop requires vernalization but
        ! it has been met,meaning fv>0.95,skip.
        IF (fv<0.95d0) THEN
            DO i=1,no_days(year,mnth,thty_dys)
                j=j+1
                ! If it has not been sowed yet, cycle so it might get to the day
                ! of sow
                IF (j<pft_tab(ft)%sowday(nn1+1)) CYCLE
                ! If vernilization not complete then calculate vernization response fv
                IF (fv<0.95d0) call streck(pft_tab(ft)%cardinal(4), &
                  pft_tab(ft)%cardinal(5),pft_tab(ft)%cardinal(6), &
                  ssp%emnthtmp(mnth,nn1+1),pft_tab(ft)%croptype(1),pft_tab(ft)%photoperiod(3), &
                  pft_tab(ft)%photoperiod(4),dayl(ssp%lat,j),vdays,fv)
                ! Calculate PHU for vern and vegetative states for the month
                totveg=totveg+ssp%emnthtmp(mnth,nn1+1)*fv*wft*wfp
                CYCLE
            endDO ! i=1,no_days(year,mnth,thty_dys)
            nvern=nvern+1
        ELSE
            ! Vegetative growth PHU for non vern plants or vern plants that have completed 
            ! vern with wft and wfp calculated from the previous call to wangengel
            totveg=totveg+ssp%emnthtmp(mnth,nn1+1)*no_days(year,mnth,thty_dys)*wft*wfp
            ! Reproductive PHU for non for non vern plants or vern plants that have completed 
            ! vern  
            CALL wangengel(pft_tab(ft)%cardinal(7),pft_tab(ft)%cardinal(8), &
              pft_tab(ft)%cardinal(9),ssp%emnthtmp(mnth,nn1+1),pft_tab(ft)%croptype(1), &
              pft_tab(ft)%photoperiod(5),pft_tab(ft)%photoperiod(6),hrs(mnth),wft,wfp)
            repphu(mnth)=repphu(mnth)+ssp%emnthtmp(mnth,nn1+1)*no_days(year,mnth,thty_dys)*wft*wfp
        endIF

        ! The cumulative vegetative heat units achieved this month both for vern and non.
        ! This does not include maturity
        vegphu(mnth)=totveg
        ! If vernalisation complete or is not involved, then we need to count the number of months where we
        ! accumulate reproductive PHU.  These months need to be continuous.
        ! It will keed adding to nwarm for as long as repphu is positive.
        ! When pepphu becomes zero or negative,it will mean that the growing season ended
        ! and mcold will attain a value diff than zero which will stop PHU from accumulating
        ! based on a conditional above.
        IF (fv>0.95d0) THEN
            ! Past the growing season; we've stopped accumulating reproductive PHU
            IF (nwarm>0.AND.repphu(mnth)<=0.0d0) mcold=mnth
            ! Growing season and warm enough for anthesis/fruiting
            IF (mcold==0.AND.repphu(mnth)>0.0d0) nwarm=nwarm+1
        ENDIF ! fv>0.995d0
        ! Total number of growing season months
        IF (mcold==0) ntot=ntot+1
    ENDDO ! k=0,11

    ! Sets the minimum GDD    
    pft_tab(ft)%cropgdd(1,nn1+1)=pft_tab(ft)%croprange(1)
    if (ssp%iseas==2) then ! Temperature controlled
        ! Here, I apply a squared cosine function to get the GDD rather than a
        ! quadratic as Bondeau et al (2007) did.  The cosine provides a smoothly
        ! varying function which can easily be shifted depending on when the coldest
        ! month occurs, avoiding any need for hardcoding time windows.

        ! If the total vegetative units are greater than the minimum
        if (totveg>pft_tab(ft)%cropgdd(1,nn1+1).and. &
          pft_tab(ft)%croprange(1)<pft_tab(ft)%croprange(2)) THEN
            nydays=daysiny
            pft_tab(ft)%cropgdd(1,nn1+1)=pft_tab(ft)%croprange(2)
            pft_tab(ft)%cropgdd(2,nn1+1)=pft_tab(ft)%croprange(4)
            ! The number of days from sow till the midday of the coldest month
            vrat=pft_tab(ft)%sowday(nn1+1)-icoldest
            ! squared cosine function.Involves the min and max gdd
            ! and the days of the growing season compared to the days
            ! of the year.This  calculates cropgdd(1,:) which is a guess
            ! on the crop gdd
            pft_tab(ft)%cropgdd(1,nn1+1)=pft_tab(ft)%croprange(1)+ &
              (pft_tab(ft)%croprange(2)-pft_tab(ft)%croprange(1))* &
              (cos(vrat*3.14159d0/nydays)**2)
        endIF ! (croprange(1,ft)<croprange(2,ft))
    ELSE
        ! Estimate using cumulative sum for vegetative growth (vegphu), max allowed
        ! months for reproductive growth (3), and phenological index where senescence
        ! begins (cropgdd(5,ft))
        ! croptype(2,ft) is zero if no vernalisation required
        ! croprange(1,ft) is min PHU for maturity; ntot is total months accumulating PHU
        ! Restrict ntot to 3 if the crop does not require vernalisation
        if (pft_tab(ft)%croptype(2)==0 &
         .and.maxval(vegphu)>pft_tab(ft)%croprange(1)) ntot=min(3,ntot)
        totrep=0;
        i=0; totphu(:)=0.0d0; slen(:)=0
        ! We are going to try several different combinations of vegetative PHU and
        ! reproductive PHU.  
        ! Start with the month after vernalisation is complete
        ! nvern are the months required for vernalization
        do k=nvern+1,ntot-1
            ! cropgdd(1,ft) is initialised to be croprange(1,ft) but may change
            ! croprange(1,ft) is the min GDD required
            if (pft_tab(ft)%cropgdd(1,nn1+1)>pft_tab(ft)%croprange(1)) CYCLE
            ! Get the current calendar month, after vernalising
            mnth=msow+k
            if (mnth>12) mnth=mnth-12
            ! cropphen(5,ft) is the phenological index when we stop adding LAI and start 
            ! adding carbon to fruiting structures instead: in other words, when we start
            ! accumulating reproductive rather than vegetative PHU
            ! targ is the target PHU for maturity, when we can harvest the plant
            ! targ is always greater than vegphu because we add repphu to it
            targ=vegphu(mnth)/pft_tab(ft)%cropphen(5)
            ! recall vegphu is the PHU accumulated by the end of the current month
            ! totphu is going to be the sum of vegphu and enough months of reproductive PHU
            ! to reach maturity
            totphu(mnth)=vegphu(mnth)
            ! loop over reproductive (fruiting) months
            do j=k+1,ntot
                ! We accumulate reproductive PHU only until we achieve the PHU for harvest 
                if (totphu(mnth)>=targ) CYCLE
                ! ok, I think msow should have been added to m, because we want the calendar
                ! month(s) following "mnth", which was the last month for vegetative PHU.  For
                ! gridcells without temperature-seasonality, however, the error may be small...
                m=j+msow
                if (m>12) m=m-12
                ! Add monthly reproductive PHU to total cumulative PHU
                totphu(mnth)=MAX(targ,totphu(mnth)+repphu(m))
                ! estimate growing season length after vernalisation assuming senescence begins
                ! in month "mnth".
                slen(mnth)=j-nvern
            endDO ! do j=k+1,ntot

            ! If the total PHU is between the min and max allowed, update cropgdd
            ! Also increment totrep, which will allow us to calculate the mean value of
            ! totphu test values that fall within the designated croprange.  The number of 
            ! totphu test values that satisfy this criterion will be held in "i".
            if (totphu(mnth)<=pft_tab(ft)%croprange(2).and. &
              totphu(mnth)>pft_tab(ft)%croprange(1)) THEN
                pft_tab(ft)%cropgdd(1,nn1+1)=totphu(mnth)
                totrep=totrep+totphu(mnth) !*4.0d0/(1+slen(mnth))
            if (slen(mnth)>0) i=i+1 !slen(mnth)+1
            elseif (totphu(mnth)>pft_tab(ft)%croprange(2)) THEN
                pft_tab(ft)%cropgdd(1,nn1+1)=pft_tab(ft)%croprange(2)
            endIF
        endDO ! k=0,11

        ! Get the final mean PHU for maturity, if there are several possibilities
 
        if (i>0) pft_tab(ft)%cropgdd(1,nn1+1)=totrep/i !*i/4.0d0
    endIF ! iseas==2, controlled by temperature or otherwise
endDO ! ft=3,nft

end subroutine seasonality


!**********************************************************************!
!                                                                      !
!                         getwet :: crops                              !
!                     ---------------------                            !
!                                                                      !
! subroutine getwet(hrs,q,wet,pet,nn1,year,thty_dys)                   !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief 
!! @details Calculate ratios of precipitation to evt
!! van Bussel (2011) discusses determination of start of growing period
!! in her chapter 5, section 2.2.2, and the sowing day should be the first
!! day where precip>0.1mm of that wet month.
!! Calculate 4-month sums of prec/evt (petsum(12) and find the highest.
!! Output the starting month of the highest (wet).
!! We depart from van Bussel in that we do not use Priestley-Taylor PET, 
!! but instead lift Penman-Montheith code from of dolyday.
!!
!! REFERENCES:
!!
!! Waha et al 2012, Global Ecology and Biogeography 21:247-259
!! van Bussel et al 2015, Global Ecology and Biogeography (earlyview)
!!
!! hrs(12)      IN    Daylight hours for the midday of each month
!! q(12)        IN    Monthly photon flux density
!! wet          OUT   Month index for start of wet season
!! pet(12)      OUT   Monthly ratio precip/(potential evapotranspiration)
!! 
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!!
subroutine getwet(hrs,q,wet,pet,nn1,year,thty_dys)
!**********************************************************************!    
implicit none
real(dp) :: hrs(12),q(12),pet(12)
integer  :: wet,nn1,year,dday,thty_dys

real(dp) :: ht,windspeed,canga,t,rh,rn
real(dp) :: svp,vpd,lam,rho,s,gam,ee,eemm,petsum(12)
integer :: m    
!*----------------------------------------------------------------------*
  
  ! Always has this value in sdgvm0.f
  ht=25.1338
  ! In m/s
  windspeed= 5.0d0
  canga = 0.168d0*windspeed/log((200.0d0 - &
    0.7d0*ht)/(0.1d0*ht))**2

  pet(:)=0
  do m=1,12
    if (ssp%emnthprc(m,nn1+1)<=0.0d0) CYCLE     
    t=ssp%emnthtmp(m,nn1+1); rh=ssp%mnthhum(m);
    
    rn = 0.96d0*(q(m)*1000000.0d0/4.0d0 + 208.0d0 + 6.0d0*t)
    rn = rn*0.52d0
    ! To be used with the new evt
    !dday = no_day(year,m,15,thty_dys)
    !*----------------------------------------------------------------------*
    !* Penman-Monteith equation for evapotranspiration.                     *
    !* Units of ET = W/m2  (CHECK) N.B. Conductances in moles.              *
    !*----------------------------------------------------------------------*
    svp = 6.108d0*exp((17.269d0*t)/(237.3d0 + t))*100.0d0
    vpd = (1.0d0 - rh/100.0d0)*svp
    lam = 2500.0d0 - 2.367d0*t
    rho = 1288.4d0 - 4.103d0*t
    s = 48.7d0*exp(0.0532d0*t)
    gam = 101325.0d0*1.012d0/(0.622d0*lam)

    ee = (s*rn + rho*1.012d0*canga*vpd)/(s + gam)
    eemm = (ee*3600.0d0*hrs(m))/(lam*1000.0d0)
    ! Divides for each month of the year the exponentialy weighted precip
    ! by the potential evapotranspiration
    pet(m)=ssp%emnthprc(m,nn1+1)/eemm !Ratio

  enddo ! m=1,12

  ! Get the four-month sums  of the above
  petsum(1)=sum(pet(1:4))
  wet=1
  do m=2,9
    petsum(m)=sum(pet(m:m+3))
    if (petsum(m)>petsum(m-1)) wet=m
  endDO ! m=1,12   
  petsum(10)=sum(pet(10:12))+pet(1)
  petsum(11)=sum(pet(11:12))+sum(pet(1:2))
  petsum(12)=pet(12)+sum(pet(1:3))
  do m=10,12
    if (petsum(m)>petsum(m-1)) wet=m
  endDO ! m=1,12   
  return
end subroutine getwet


!**********************************************************************!
!                                                                      !
!                        summerday :: crops                            !
!                     ---------------------                            !
!                                                                      !
! function summerday                                                   !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Find the start day of summer as defined by a threshold temperature
!! @details van Bussel (2011, thesis) discusses determination of sowing dates
!! in her chapter 5, section 2.2.2; here we find the day number when
!! summer begins (temperature>=threshold).  We linearly interpolate
!! between the average monthly temperatures to get the mean day.
!!
!!
!! thresh       IN threshold temperature for summer
!! mmid(12)     IN Midday of month in Julian days
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!!
integer function summerday(thresh,mmid,nn1)
!**********************************************************************!
implicit none
real(dp) :: thresh
integer  :: mmid(12)
integer  :: lasttmp,lastm,lastday,m,thistmp,nn1


  summerday=1
  ! If all monthly temperatures below or above threshold,return
  IF (minval(ssp%emnthtmp(:,nn1+1))>=thresh) RETURN
  IF (maxval(ssp%emnthtmp(:,nn1+1))<thresh) RETURN

  lasttmp = 0; lastm = 12; lastday = -mmid(1)
  IF (ssp%emnthtmp(12,nn1+1)>=thresh) lasttmp = 1

 
  DO m=1,12
    thistmp = 0
    IF (ssp%emnthtmp(m,nn1+1)>=thresh) thistmp = 1
    ! Does the linear interpolation
    IF (thistmp>lasttmp) summerday = &
      !Finds the difference between the threshold temperature and the temperature 
      !of the previous lowest month and multiplies it by
      INT((thresh-ssp%emnthtmp(lastm,nn1+1))* &
        !the rate of day/temperature to find the number of days required to reach
        !threshold plus adding the last day to get julian day
        (mmid(m)-lastday)/(ssp%emnthtmp(m,nn1+1)-ssp%emnthtmp(lastm,nn1+1)))+lastday
    lasttmp = thistmp; lastm = m; lastday = mmid(m)
  ENDDO
  IF (summerday<=0) summerday = mmid(12)-summerday
      
  return 
end function summerday


!**********************************************************************!
!                                                                      !
!                       winterday :: crops                             !
!                     ---------------------                            !
!                                                                      !
! function winterday                                                   !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Find the start day of winter as defined by a threshold temperature
!! @details 
!!
!!
!!
!! thresh       IN threshold temperature for winter
!! mmid(12)     IN Midday of month in Julian days
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!!
integer function winterday(thresh,mmid,nn1)
!**********************************************************************!
implicit none
real(dp) :: thresh
integer  :: mmid(12)
integer  :: lasttmp,lastm,lastday,m,thistmp,nn1

  winterday=1
  if (minval(ssp%emnthtmp(:,nn1+1))>=thresh) RETURN
  if (maxval(ssp%emnthtmp(:,nn1+1))<thresh) RETURN
  lasttmp=0; lastm=12; lastday=-mmid(1)
  if (ssp%emnthtmp(12,nn1+1)<=thresh) lasttmp=1
  do m=1,12
    thistmp=0
    if (ssp%emnthtmp(m,nn1+1)<=thresh) thistmp=1
    if (thistmp>lasttmp) winterday=int((thresh-ssp%emnthtmp(lastm,nn1+1))* &
      (mmid(m)-lastday)/(ssp%emnthtmp(m,nn1+1)-ssp%emnthtmp(lastm,nn1+1)))+lastday
    lasttmp=thistmp; lastm=m; lastday=mmid(m)
  endDO
  if (winterday<=0) winterday=mmid(12)-winterday
      
  RETURN
  end function winterday


!**********************************************************************!
!                                                                      !
!                       wangengel :: crops                             !
!                     ---------------------                            !
!                                                                      !
! subroutine wangengel                                                 !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Get temperature and photoperiod response function values for crops
!! @details These functions are from Wang and Engel 1998 (Agri Sys 58:1-24).
!! The sum of the product of the temperature and photoperiod functions gives 
!! the effective physiological days.
!!
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!!
subroutine wangengel(tmin,topt,tmax,tmp,dtype,popt,pcrit,p,ft,fp)
!**********************************************************************!
implicit none
real(dp) :: tmin,topt,tmax ! Cardinal temperatures
real(dp) :: popt,pcrit,p ! Optimal, critical and current photoperiod
real(dp) :: tmp     ! today's mean temperature, input
real(dp) :: ft,fp   ! temperature and photoperiod functions, output
real(dp) :: alpha  ! exponent in temperature response function, internal
real(dp) :: t1,t2 ! internal
real(dp) :: dtype ! 1 for long-day, -1 for short-day, 0 for day-neutral plants

  
  ! Photoperiod effect.Eq. 11,12,13 from Wang & Engel
  fp=max(0.0,1-exp(-dtype*4.0*(p-pcrit)/abs(popt-pcrit)))
  ! Temperature effect.
  ft=0.0
  if (tmax>tmin) then
    ! NB ft will be zero for tmp outside the range tmin to tmax inclusive
    if (tmp<=tmax.and.tmp>=tmin) then
      ! alpha value from Eq.6 of Wang & Engel
      alpha=log(2.0)/log((tmax-tmin)/(topt-tmin))
      t1=(topt-tmin)**alpha
      t2=(tmp-tmin)**alpha
      ft=(2*t2*t1-(t2**2))/(t1**2)
    endIF ! tmp<=tmax.and.tmp>=tmin
  ELSE
  
  endIF ! vtmax>vtmin

return
  
end subroutine wangengel


!**********************************************************************!
!                                                                      !
!                       streck :: crops                                !
!                     ---------------------                            !
!                                                                      !
! subroutine streck                                                    !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Get vernalisation factors for crops
!! @details  Following Streck et al 2003, Ag For Meteor 115:139-150
!! Equations 3 and 4 (temperature function), 5 (photoperiod function)
!! These functions are from Wang and Engel 1998 (Agri Sys 58:1-24).
!! In some papers including Streck et al's first one (below) the temperature
!! function has a typographical error (a minus sign is missing).
!! The sum of the product of the temperature and photoperiod functions gives 
!! the effective vernalised days VD
!! which one can use to get the vernalization response function (Eq 10)
!! or Eq 2 in Streck et al 2003 Agron J 95:155-159
!!
!! Streck et al's vernalisation function and Wang and Engel's functions
!! are used in a number of models:
!!
!! JULES-SUCROS land surface model 
!!    (van den Hoof et al 2011 Ag For Meteor 151:137-153)
!! DANUBIA crop growth model (Lenz-Wiedemann et al 2010 Ecol Model 221:314-329)
!! FROSTOL wheat model (Bergjord et al 2008 Eur J Agr 28:321-330)
!! SPACSYS C and N cycling model (Bingham and Wu 2011 Eur J Agr 34:181-189) 
!! CANDY-PLUS C, N and biomass model (Kruger et al 2013 Vadose Zone J 12)
!! Liu's pasture legume model (Liu 2007 Field Crops Res 101:331-342)
!! 
!! vtmin,vtopt,vtmax  input Cardinal temperatures for vernalisation
!! tmp                input Current mean daily temperature
!! popt,pcrit         input Optimal and critical photoperiod
!! p                  input Current photoperiod
!! dtype              input -1: short-day, 1: long-day, 0: day-neutral
!! vdays       input/output Accumulated vernalisation days
!! fv                output Current vernalisation day
!!
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!!
subroutine streck(vtmin,vtopt,vtmax,tmp,dtype,popt,pcrit,p,vdays,fv)
!**********************************************************************!
implicit none
real(dp) :: vtmin,vtopt,vtmax
real(dp) :: tmp     ! today's mean temperature, input
real(dp) :: vdays   ! accumulated effective vernalisation days, input/output
real(dp) :: valpha  ! exponent in temperature response function, internal
real(dp) :: dtype ! 1 for long-day, -1 for short-day, 0 for day-neutral plants
real(dp) :: popt,pcrit ! photoperiod parameters, input
real(dp) :: p ! today's photoperiod, input
real(dp) :: vd5,ft,fp   ! internal
real(dp) :: fv   ! vernalisation response function, output

if (vtmax>vtmin) then
! NB increment will be zero for tmp outside the range vtmin to vtmax inclusive
! therefore the heat sum vdays will not be incremented and fv will not change
  if (tmp<=vtmax.and.tmp>=vtmin) then
    call wangengel(vtmin,vtopt,vtmax,tmp,dtype,popt,pcrit,p,ft,fp)
    vdays=vdays+ft*fp
    vd5=vdays**5
    fv=(vd5)/((22.5**5)+vd5)
  endif ! tmp<=vtmax.and.tmp>=vtmin
else
  vdays=0.0d0; fv=1.0d0
endif ! vtmax>vtmin
 
return
end subroutine streck

!**********************************************************************!
!                                                                      !
!                       crop_outputs :: crops                          !
!                     ---------------------                            !
!                                                                      !
! subroutine crop_outputs                                              !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Output files for crops
!! @details  
!!
!! @author EPK 
!! @date Dec 2016
!----------------------------------------------------------------------!!
subroutine crop_outputs(nft,sl)
!**********************************************************************!
integer :: nft,ft,sl,coun,fid_c

coun = 0
fid_c = 5500

do ft=3,nft
  if(pft_tab(ft)%phen/=3) cycle
  select case (sl)
  case (0)
    coun=coun+1
    open(fid_c+coun,file=trim(inp%dirs%output)//'/'// &
      trim(pft_tab(ft)%tag)//'_sow.dat')
    coun=coun+1
    open(fid_c+coun,file=trim(inp%dirs%output)//'/'// &
      trim(pft_tab(ft)%tag)//'_harv.dat')
    coun=coun+1
    open(fid_c+coun,file=trim(inp%dirs%output)//'/'// &
      trim(pft_tab(ft)%tag)//'_gdd.dat')
    coun=coun+1
    open(fid_c+coun,file=trim(inp%dirs%output)//'/'// &
      trim(pft_tab(ft)%tag)//'_seas.dat')
    coun=coun+1
    open(fid_c+coun,file=trim(inp%dirs%output)//'/'// &
      trim(pft_tab(ft)%tag)//'_ryield.dat')
  case (1)
    coun=coun+1
    close(fid_c+coun)
    coun=coun+1
    close(fid_c+coun)
    coun=coun+1
    close(fid_c+coun)
    coun=coun+1
    close(fid_c+coun)
    coun=coun+1
    close(fid_c+coun)
  case (2)
    coun=coun+1
    write(fid_c+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon
    coun=coun+1
    write(fid_c+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon
    coun=coun+1
    write(fid_c+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon
    coun=coun+1
    write(fid_c+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon
    coun=coun+1
    write(fid_c+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon 
  case (3)
    coun=coun+1
    write(fid_c+coun,'(i4)',advance='NO') pft_tab(ft)%sowday(1)
    coun=coun+1
    if(ssp%co2ftmap(ft,1)>0) THEN
      write(fid_c+coun,'('' '',i4)',advance='NO') ssv(ssp%co2ftmap(ft,2))%harvest(2)
    ELSE
      write(fid_c+coun,'('' '',i4)',advance='NO') 0
    endIF
    coun=coun+1
    write(fid_c+coun,'(i5)',advance='NO') pft_tab(ft)%cropgdd(1,1)
    coun=coun+1
    write(fid_c+coun,'(i2)',advance='NO') ssp%iseas
    coun=coun+1
    if(ssp%co2ftmap(ft,1)>0) THEN
      write(fid_c+coun,'('' '',f9.2)',advance='NO') ssv(ssp%co2ftmap(ft,2))%yield
    ELSE
      write(fid_c+coun,'('' '',f9.2)',advance='NO') 0.
    endIF
  case (4)
    coun=coun+1
    write(fid_c+coun,*)
    coun=coun+1
    write(fid_c+coun,*)
    coun=coun+1
    write(fid_c+coun,*)
    coun=coun+1
    write(fid_c+coun,*)
    coun=coun+1
    write(fid_c+coun,*)
  end select

enddo
return
end subroutine crop_outputs

!**********************************************************************!
!                                                                      !
!                      irrigate :: crops                               !
!                     ---------------------                            !
!                                                                      !
! subroutine irrigate()                                                !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Decides irrigation and fills up relative soil layers. 
!! @details sfc is soil field capacity in mm
!! sw  is soil wilting point in mm
!! ssv(ft)%soil_h2o is soil water content in mm
!! We should also consider runoff.Ideally,we would check runoff before
!! irrigating as it would give us the water available for irrigation.
!! We would then subtract the water used for irrigation from runoff.
!! 
!! @author EPK,LLT
!! @date Jan 2017
!----------------------------------------------------------------------!!
SUBROUTINE IRRIGATE(ft,sfc,sw)
!**********************************************************************!
  IMPLICIT NONE

  real(dp) :: sfc(4),sw(4)
  real(dp) :: sumh,sumn,sumirr,fr_irr,sumi
  integer :: ft,i
  !Layers I will check to determine irrigation
  INTEGER,PARAMETER,DIMENSION(1) :: ll1=[2]
  !Layers I am filling up with irrigation
  INTEGER,PARAMETER,DIMENSION(1) :: ll2=[2]
 
  ! Return if its not crop phenology
  if(pft(ft)%phen/=3) RETURN
  ! Return if not sown
  !WRITE(*,*)ssp%day+(ssp%mnth-1)*30,pft(ft)%sowday(3),ssv(ft)%sown,ssv(ft)%harvest&
  !  ,ssv(ft)%phu,ssv(ft)%lai%tot(1)
  if(ssv(ft)%sown==0) RETURN

  ! Sum up the available water in the soil layers defined 
  ! with the parameter ll1.
  sumh=0.
  do i=1,SIZE(ll1)
    sumh=sumh+ssv(ft)%soil_h2o(ll1(i))
  endDO
  
  ! Fraction of the crop that is irrigated as read from file
  fr_irr=pft(ft)%irrig(3)
  
  ! Sum the water in the soil layers defined with the parameter ll1
  ! below which irrigation is triggered.This is controlled by the 
  ! parameter pft(ft)%irrig(1)=-1 and the fraction of crop that is 
  ! irrigated fr_irr.The smaller pft(ft)%irrig(1), the sooner irrigation
  ! kicks in.sw is wilting points for the different soil layers
  ! in mm of water.
  sumn=0.
  do i=1,SIZE(ll1)
    !sumn=sumn+sw(ll1(i))+(1-EXP(pft(ft)%irrig(1)*fr_irr))*(sfc(ll1(i))-sw(ll1(i)))
    sumn=sumn+(1-EXP(pft(ft)%irrig(1)*fr_irr))*sfc(ll1(i))
  endDO
  
  ! Decides irrigation
  if(sumh>sumn) RETURN

  ! If irrigation is happening,irrigate soil layers defined with the
  ! parameter ll2 which can be different that ll1 which is used to
  ! define which soil layers trigg irr.The water added will be again 
  ! controlled by the parameter pft(ft)%irrig(2)=-2 and the fraction
  ! of crop that is irrigated fr_irr.
  sumirr=0.
  do i=1,SIZE(ll2)
    ! This is the water I want in the soil layer after irrigation
    !sumi=sw(ll2(i))+(1-EXP(pft(ft)%irrig(2)*fr_irr))*(sfc(ll2(i))-sw(ll2(i)))
    sumi=(1-EXP(pft(ft)%irrig(2)*fr_irr))*sfc(ll2(i))
   ! If the water I want in greater than the water in the layer
    if(sumi>ssv(ft)%soil_h2o(ll2(i))) THEN
      sumirr=sumi-ssv(ft)%soil_h2o(ll2(i))
      ssv(ft)%soil_h2o(ll2(i))=sumi
    endIF
  endDO 
  
end SUBROUTINE IRRIGATE


!**********************************************************************!
!                                                                      !
!                      fert_crops :: crops                             !
!                     ---------------------                            !
!                                                                      !
! subroutine fert_crops(nft)                                           !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Effect of fertilizer usage on optimal crop LAI
!! @details Finds pft_tab(ft)%optlai
!!
!! @author LLT,EPK
!! @date Jan 2017
!----------------------------------------------------------------------!!
subroutine fert_crops(nft)
!**********************************************************************!
implicit none

integer nft,ft

do ft=3,nft
  if(pft_tab(ft)%phen/=3) cycle
  ! cropphen gives the optimal LAI range: cropphen(1) is without fertiliser
  ! cropphen(2) is with maximum fertiliser
  ! if this crop doesn't normally get fertiliser, then cropphen(2)=cropphen(1)
  pft_tab(ft)%optlai=pft_tab(ft)%cropphen(1)  
  ! If the crop is applied with fertilizer
  if(pft_tab(ft)%cropphen(2)>pft_tab(ft)%cropphen(1)) THEN
    ! get optimal LAI by adding log fert usage (kg/ha) to the nonfertilised LAI
    ! 1g/m2=10kg/hc
    ! kg/ha = 1000g/10000m2 so 10 x g/m2 = kg/ha = g/10m2
    ! Make sure optlai remains in the range specified by cropphen
    pft_tab(ft)%optlai=pft_tab(ft)%cropphen(1)+&
      (pft_tab(ft)%cropphen(2)-pft_tab(ft)%cropphen(1))*(1-EXP(pft_tab(ft)%fert(7)*pft_tab(ft)%fert(1)))
    pft_tab(ft)%optlai=min(pft_tab(ft)%cropphen(2),pft_tab(ft)%optlai)
  endIF
endDO


end subroutine fert_crops

!**********************************************************************!
!                                                                      !
!                      READ_FERTILIZERS :: crops                       !
!                     --------------------------                       !
!                                                                      !
! subroutine READ_FERTILIZERS(du,yr0,yrf,lat,lon,nft,cfert)            !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Read fertilizer data following land cover read data approach
!! @details Output is cfert(nft,years,NPK).Units are kg/ha,scale is 0.1
!!
!! @author EPK
!! @date Feb 2017
!----------------------------------------------------------------------!!
SUBROUTINE READ_FERTILIZERS (du,yr0,yrf,lat,lon,nft,cfert)
                             

REAL(dp) :: latf,lon0,latr,lonr,rrow,rcol,ynorm,xnorm,lat,lon,ans,xx(4,4)
REAL(dp) :: ftprop(max_cohorts,3),cfert(max_cohorts,max_years,3)
INTEGER :: latn,lonn,kode,n,years(1000),nclasses,classes(1000),ift
INTEGER :: du,nrecl,yr0,yrf,j,i,k,l,col,row,recn,nft,x,ii,jj,indx(4,4)
CHARACTER(len=str_len) :: st1,st2,st3
CHARACTER(len=1),dimension(3) :: f_typ=(/'N','P','K'/)

!----------------------------------------------------------------------!!
! read in the readme file 'readme.dat'.                                !
!----------------------------------------------------------------------!!
open(99,file=trim(inp%dirs%fert)//'/readme.dat',&
 status='old',iostat=kode)
if(kode/=0) THEN
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Fertilizer file does not exist.'
  write(*,'('' "'',A,''/readme.dat"'')') trim(inp%dirs%fert)
  stop
endIF

READ(99,*) st1
st2='CONTINUOUS'
if(stcmp(st1,st2)==0) THEN
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Fertilizers is not a continuous field ?'
  write(*,*) 'readme.dat should begin with CONTINUOUS'
  stop
endIF

READ(99,*)
READ(99,*) latf,lon0
READ(99,*)
READ(99,*) latr,lonr
READ(99,*)
READ(99,*) latn,lonn
READ(99,*)
READ(99,'(A)') st1
n = n_fields(st1)
CALL ST2ARR(st1,years,1000,n)
READ(99,*)
READ(99,'(A)') st1
close(99)
nclasses = n_fields(st1)
! nclasses is the number of fts provided and classes the ft id
CALL ST2ARR(st1,classes,1000,nclasses)

!----------------------------------------------------------------------!!
if(du==1) THEN
!  This works for ftn95
!  nrecl = 3
  nrecl = 5
ELSE
  nrecl = 4
endIF

! We need to have a map for the first year of the run
if((n>1).and.(yr0<years(1))) THEN
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Can''t start running in ',yr0,&
 ' since fertilizer map begin in ',years(1)
  stop
endIF


!----------------------------------------------------------------------!!
! Find the real(dp) :: row col corresponding to lat and lon.           !
!----------------------------------------------------------------------!!
! rrow and rcol is the gridcell number for these lats and lons,in decimal
rrow = (latf - lat)/latr
rcol = (lon - lon0)/lonr

! ynorm and xnorm is the remainder of the lat and lon
ynorm = rrow - real(int(rrow))
xnorm = rcol - real(int(rcol))
!----------------------------------------------------------------------!!

ftprop=0.

j=1
! For each counting year of the run
do i=1,yrf-yr0+1
  !years holds the fertilizer years available 
  !If you are looking at the first year or any year when we have data,read file
  if((i==1).or.((i+yr0-1)==years(j))) then
    st2=in2st(years(j))
    st2 = adjustl(st2)
    j=j+1
    ! For each of the classes I have fertilizer data,look for file
    do k=1,nclasses

      st3=in2st(classes(k))
      st3 = adjustl(st3)
      ! For each of the 3 fertilizers,look for file 
      do l=1,size(f_typ,1)
        open(99,file= &
   trim(inp%dirs%fert)//'/cont_fert_'//f_typ(l)//'_'//trim(st3)//'_'//st2(1:4)//'.dat', & 
   status='old',form='formatted',access='direct',recl=nrecl,iostat=kode)
        if(kode/=0) THEN
          write(*,'('' PROGRAM TERMINATED'')')
          write(*,*) 'Fertilizer data-base.'
          write(*,*) 'File does not exist:'
          write(*,*) trim(inp%dirs%fert)//'/cont_fert_'//f_typ(l)//'_'//trim(st3),'_',st2(1:4),'.dat'
          stop
        endIF
        
              
        do ii=1,4
          do jj=1,4
            row = INT(rrow)+jj-1
            col = INT(rcol)+ii-1
            if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) then
              recn = (row-1)*lonn + col
              READ(99,'(i3)',rec=recn) x 
              xx(ii,jj) = REAL(x)
              if (x<200) then
                indx(ii,jj) = 1
              ELSE
                indx(ii,jj) = 0
              endIF
            ELSE
              indx(ii,jj) = -1
            endIF
          endDO
        endDO

        CALL bi_lin(xx,indx,xnorm,ynorm,ans)

        x = INT(ans+0.5)
        
        ftprop(classes(k),l) = ans
        close(99)


      endDO !End of loop over the fertilizer types
    endDO ! End of loop over the classes
  endIF ! Finished reading files
  
  ! Assign fertilizer to cfert(ft,years,fert)
  ! If a year doesn't exist in the map files,it will use the value of the previous year
  do ift=1,nft
    cfert(ift,i,:) = ftprop(ift,:)
  endDO
  
endDO ! End of year loop





end SUBROUTINE READ_FERTILIZERS



!**********************************************************************!
!                                                                      !
!                      READ_IRRIGATION :: crops                        !
!                     --------------------------                       !
!                                                                      !
! subroutine READ_IRRIGATION(du,yr0,yrf,lat,lon,nft,cirr)              !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Read irrigation data following land cover read data approach 
!! @details Output is cirr(nft,years).Units are % of land irrigated
!!
!! @author EPK
!! @date Feb 2017
!----------------------------------------------------------------------!!
SUBROUTINE READ_IRRIGATION (du,yr0,yrf,lat,lon,nft,cirr)

REAL(dp) :: latf,lon0,latr,lonr,rrow,rcol,ynorm,xnorm,lat,lon,xx(4,4),ans
REAL(dp) :: ftprop(max_cohorts),cirr(max_cohorts,max_years)
INTEGER :: latn,lonn,kode,n,years(1000),nclasses,classes(1000),ift
INTEGER :: du,nrecl,yr0,yrf,j,i,k,col,row,x,recn,nft,ii,jj,indx(4,4)
CHARACTER(len=str_len) :: st1,st2,st3

!----------------------------------------------------------------------!!
! read in the readme file 'readme.dat'.                                !
!----------------------------------------------------------------------!!
open(99,file=trim(inp%dirs%irri)//'/readme.dat',&
 status='old',iostat=kode)
if(kode/=0) THEN
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Irrigation file does not exist.'
  write(*,'('' "'',A,''/readme.dat"'')') trim(inp%dirs%irri)
  stop
endIF

READ(99,*) st1
st2='CONTINUOUS'
if(stcmp(st1,st2)==0) THEN
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Irrigation is not a continuous field ?'
  write(*,*) 'readme.dat should begin with CONTINUOUS'
  stop
endIF

READ(99,*)
READ(99,*) latf,lon0
READ(99,*)
READ(99,*) latr,lonr
READ(99,*)
READ(99,*) latn,lonn
READ(99,*)
READ(99,'(A)') st1
n = n_fields(st1)
CALL ST2ARR(st1,years,1000,n)
READ(99,*)
READ(99,'(A)') st1
close(99)
nclasses = n_fields(st1)
! nclasses is the number of fts provided and classes the ft id
CALL ST2ARR(st1,classes,1000,nclasses)

!----------------------------------------------------------------------!!
if(du==1) THEN
!  This works for ftn95
!  nrecl = 3
  nrecl = 5
ELSE
  nrecl = 4
endIF

! We need to have a map for the first year of the run
if((n>1).and.(yr0<years(1))) THEN
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Can''t start running in ',yr0,&
 ' since irrigation map begins in ',years(1)
  stop
endIF


!----------------------------------------------------------------------!!
! Find the real(dp) :: row col corresponding to lat and lon.           !
!----------------------------------------------------------------------!!
! rrow and rcol is the gridcell number for these lats and lons,in decimal
rrow = (latf - lat)/latr
rcol = (lon - lon0)/lonr

! ynorm and xnorm is the remainder of the lat and lon
ynorm = rrow - real(int(rrow))
xnorm = rcol - real(int(rcol))
!----------------------------------------------------------------------!!

ftprop=0.

j=1
! For each counting year of the run
do i=1,yrf-yr0+1
  ! years holds the irrigation years available 
  ! If you are looking at the first year or any year when we have data,read file
  if((i==1).or.((i+yr0-1)==years(j))) then
    st2=in2st(years(j))
    st2 = adjustl(st2)
    j=j+1
    ! For each of the classes I have irrigation data,look for file
    do k=1,nclasses

      st3=in2st(classes(k))
      st3 = adjustl(st3)
      open(99,file= &
 trim(inp%dirs%irri)//'/cont_irr_'//trim(st3)//'_'//st2(1:4)//'.dat', & 
 status='old',form='formatted',access='direct',recl=nrecl,iostat=kode)
      if(kode/=0) THEN
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'Irrigation data-base.'
        write(*,*) 'File does not exist:'
        write(*,*) trim(inp%dirs%irri),'/cont_irr_',trim(st3),'_',st2(1:4),'.dat'
        stop
      endIF
      
      do ii=1,4
        do jj=1,4
          row = INT(rrow)+jj-1
          col = INT(rcol)+ii-1
          if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) THEN
            recn = (row-1)*lonn + col
            READ(99,'(i3)',rec=recn) x 
            xx(ii,jj) = REAL(x)
            if (x<200) then
              indx(ii,jj) = 1
            ELSE
              indx(ii,jj) = 0
            endIF
          ELSE
            indx(ii,jj) = -1
          endIF
        endDO
      endDO

      CALL bi_lin(xx,indx,xnorm,ynorm,ans)
      
      x = INT(ans+0.5)
          
      ftprop(classes(k)) = ans
      close(99)


    endDO ! End of loop over the classes
  endIF ! Finished reading files
  
  ! Assign irrigation to cirr(ft,years)
  ! If a year doesn't exist in the map files,it will use the value of the previous year
  do ift=1,nft
    cirr(ift,i) = ftprop(ift)
  endDO
  
endDO ! End of year loop

end SUBROUTINE READ_IRRIGATION

!**********************************************************************!
!                                                                      !
!                          harv :: crops                               !
!                     ---------------------                            !
!                                                                      !
! subroutine harv(lit)                                                 !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief Harvest of crops
!! @details 
!!
!! @author EPK
!! @date Mar 2017
!----------------------------------------------------------------------!!
SUBROUTINE HARV(lit,yielit)
!**********************************************************************!
IMPLICIT NONE

REAL(dp) :: lit,hi,yielit,frc
INTEGER :: co,i

  co = ssp%cohort

  ! Calculates harvest index
  !hi=0.9*pft(co)%harvindx+0.1*pft(co)%harvindx*pft(co)%fert(1)/pft(co)%fert(5)
  !hi=MIN(hi,pft(co)%harvindx)
  ! Since you can have more than 2 yields in a calendar year,it is
  ! calculated here cumulative and is being reset to 0 at the start of each
  ! calendar year in phenology3 sub
  ! Roots which are bio(2) do not account in crop yield

  frc = pft(co)%harvip(1)*pft(co)%harvip(3)

  ssv(co)%yield=ssv(co)%yield+frc* &
    (ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
    ssv(co)%root%tot(1) + ssv(co)%stem%tot(1) + & 
    ssv(co)%bio(1) + ssv(co)%bio(2))/(0.45*pft(co)%harvip(2))

!  ! Finds yield by adding canopy,nppstore,stem and dead steam multiplied by hi
!  ssv(co)%yield=ssv(co)%yield+hi* &
!   (ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
!    ssv(co)%stem%tot(1) + ssv(co)%bio(1))
  
  yielit = frc* &
    (ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
    ssv(co)%root%tot(1) + ssv(co)%stem%tot(1) + & 
    ssv(co)%bio(1) + ssv(co)%bio(2))

!  !Same as above only gets its own variable
!  yielit=hi* &
!   (ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
!    ssv(co)%stem%tot(1) + ssv(co)%bio(1))

  ! Adds up the carbon pools that will send to litter
  lit = 0.0

  ! Sends some root carbon to litter
  if (ssv(co)%root%no > 0) then
    do i=1,ssv(co)%root%no
      lit=lit+(1-frc)*ssv(co)%root%c(i,1)%val
      ssv(co)%root%c(i,1)%val = 0.0
      ssv(co)%root%c(i,1)%age = 0.0
    endDO
    ssv(co)%root%no=0
    ssv(co)%root%tot(1)=0.0
  endIF


!  ! Sends all root carbon to litter
!  if (ssv(co)%root%no > 0) then
!    do i=1,ssv(co)%root%no
!      lit=lit+ssv(co)%root%c(i,1)%val
!      ssv(co)%root%c(i,1)%val = 0.0
!      ssv(co)%root%c(i,1)%age = 0.0
!    endDO
!    ssv(co)%root%no=0
!    ssv(co)%root%tot(1)=0.0
!  endIF

  
  ! Sends (1-hi) stem carbon to litter
  if (ssv(co)%stem%no > 0) then
    do i=1,ssv(co)%stem%no
      lit=lit+(1-frc)*ssv(co)%stem%c(i,1)%val
      ssv(co)%stem%c(i,1)%val = 0.0
      ssv(co)%stem%c(i,1)%age = 0.0
    endDO
    ssv(co)%stem%no=0
    ssv(co)%stem%tot(1)=0.0
  endIF

  
!  ! Sends (1-hi) stem carbon to litter
!  if (ssv(co)%stem%no > 0) then
!    do i=1,ssv(co)%stem%no
!      lit=lit+(1-hi)*ssv(co)%stem%c(i,1)%val
!      ssv(co)%stem%c(i,1)%val = 0.0
!      ssv(co)%stem%c(i,1)%age = 0.0
!    endDO
!    ssv(co)%stem%no=0
!    ssv(co)%stem%tot(1)=0.0
!  endIF

 
  ! Sends (1-hi) leaf carbon to litter
  if (ssv(co)%lai%no > 0) then
    do i=1,ssv(co)%lai%no
      lit=lit+(1-frc)*ssv(co)%lai%c(i,1)%val*12.0/pft(co)%sla/25.0
      ssv(co)%lai%c(i,1)%val = 0.0
      ssv(co)%lai%c(i,1)%age = 0.0
    endDO
    ssv(co)%lai%no=0
    ssv(co)%lai%tot(1)=0.0
  endIF

 
!  ! Sends (1-hi) leaf carbon to litter
!  if (ssv(co)%lai%no > 0) then
!    do i=1,ssv(co)%lai%no
!      lit=lit+(1-hi)*ssv(co)%lai%c(i,1)%val*12.0/pft(co)%sla/25.0
!      ssv(co)%lai%c(i,1)%val = 0.0
!      ssv(co)%lai%c(i,1)%age = 0.0
!    endDO
!    ssv(co)%lai%no=0
!    ssv(co)%lai%tot(1)=0.0
!  endIF

  
  ! Sends (1-hi) nppstore carbon to litter
  lit=lit+(1-frc)*ssv(co)%nppstore(1)
  ssv(co)%nppstore(1)=0.
  ssv(co)%nppstore(2)=0.


  
!  ! Sends (1-hi) nppstore carbon to litter
!  lit=lit+(1-hi)*ssv(co)%nppstore(1)
!  ssv(co)%nppstore(1)=0.
!  ssv(co)%nppstore(2)=0.
  
  
  ! Sends (1-hi) dead stem to litter and all dead roots
  lit=lit+(1-frc)*(ssv(co)%bio(1)+ssv(co)%bio(2))  
  ssv(co)%bio(1)=0.
  ssv(co)%bio(2)=0.


  
!  ! Sends (1-hi) dead stem to litter and all dead roots
!  lit=lit+(1-hi)*ssv(co)%bio(1)+ssv(co)%bio(2)  
!  ssv(co)%bio(1)=0.
!  ssv(co)%bio(2)=0.

!  IF(lit>50) THEN
!    ssv(co)%nppstore(1)=50.
!    lit=lit-50
!  ELSE
!    ssv(co)%nppstore(1)=lit
!    lit=0.
!  ENDIF

end SUBROUTINE harv

!**********************************************************************!
!                                                                      !
!                      READ_OPT_PAR :: crops                           !
!                     ---------------------                            !
!                                                                      !
! subroutine READ_OPT_PAR()                                            !
!                                                                      !
!----------------------------------------------------------------------!!
!> @brief
!! @details
!!
!! @author EPK
!! @date MAR 2017
!----------------------------------------------------------------------!!
SUBROUTINE READ_OPT_PAR()
!**********************************************************************!
IMPLICIT NONE

open(UNIT=669,FILE='opt_par.dat',STATUS='OLD')
READ(669,*)pft_tab(10)%fert(4),pft_tab(10)%fert(6)
close(669)


end SUBROUTINE READ_OPT_PAR

end module crops

