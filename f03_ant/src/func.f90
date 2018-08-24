module func

use real_precision
use dims
use pft_parameters

implicit none

contains

!**********************************************************************!
!                                                                      !
!                     check_ft_grow :: func                            !
!                     ---------------------                            !
!                                                                      !
! integer function check_ft_grow(tmp,xchill,xdschill,ft)               !
!                                                                      !
!----------------------------------------------------------------------!
integer function check_ft_grow(tmp,xchill,xdschill,ft)
!**********************************************************************!
real(dp) :: tmp(12,31),bbsum
integer :: jday,day,mnth,xchill,xdschill,chill,dschill,ft,i
!----------------------------------------------------------------------!

chill = xchill
dschill = xdschill

!----------------------------------------------------------------------!
! Check for budburst using degree days.                                !
!----------------------------------------------------------------------!
if (pft_tab(ft)%bbmem>0) then
jday = 0
10 jday = jday + 1

!----------------------------------------------------------------------!
! Check for chilling.                                                  !
!----------------------------------------------------------------------!
  if (chill==0) then
    bbsum = 0.0
    do i=1,20
      mnth = (jday+i-1)/30 + 1
      day = jday + i - (mnth - 1)*30
      mnth = mod(mnth-1,12) + 1
      if (tmp(mnth,day)<-5.0)  bbsum = bbsum + &
 max(-10.0,tmp(mnth,day)+5.0)
    enddo
    if (bbsum<-100) then
      chill = 1
      dschill = 1
    endif
  endif
  if (chill==1) then
    dschill = dschill + 1
  endif
  if (dschill>260) then
    chill = 0
    dschill = 0
  endif

  bbsum = 0.0
  do i=1,pft_tab(ft)%bbmem
    mnth = (jday+i-1)/30 + 1
    day = jday + i - (mnth - 1)*30
    mnth = mod(mnth-1,12) + 1
    if (tmp(mnth,day)>pft(ft)%bb0)  bbsum = bbsum + min(pft_tab(ft)%bblim,tmp(mnth,day)-pft_tab(ft)%bb0)
  enddo

  if (real(bbsum)>=pft_tab(ft)%bblim*exp(-0.01*real(dschill))) then
    check_ft_grow = 1
    goto 20
  endif
  if ((jday>=360+pft_tab(ft)%bbmem)) then
    check_ft_grow = 0
    goto 20
  endif
goto 10
20 continue
else
  check_ft_grow = 1
endif

end function check_ft_grow





!**********************************************************************!
!                                                                      !
!                     minx :: func                                     !
!                     ------------                                     !
!                                                                      !
!  integer function minx(x)                                            !
!                                                                      !
!----------------------------------------------------------------------!
integer function minx(x)
!**********************************************************************!
real(dp) :: x(12),minxx
integer :: i
!----------------------------------------------------------------------!

minx = -10
minxx = 100000000.0
do i=1,12
  if (x(i)<minxx) then
    minxx = x(i)
    minx = i
  endif
enddo


end function minx





!**********************************************************************!
!                                                                      !
!                     no_days :: func                                  !
!                     ---------------                                  !
!                                                                      !
! integer function no_days(year,mnth,thty_dys)                         !
!                                                                      !
!----------------------------------------------------------------------!
integer function no_days(year,mnth,thty_dys)
!**********************************************************************!
integer :: year,mnth,thty_dys
!----------------------------------------------------------------------!

if (thty_dys==0) then
  if (mnth==2) then
    if (mod(year,4)==0) then
      no_days = 29
    else
      no_days = 28
    endif
    if (mod(year,100)==0) no_days = 28
    if (mod(year,400)==0) no_days = 29
  else
    if ((mnth==4).or.(mnth==6).or.(mnth==9).or.(mnth==11)) then
      no_days = 30
    else
      no_days = 31
    endif
  endif
else
  no_days = 30
endif

end function no_days





!**********************************************************************!
!                                                                      !
!                          no_day:: func                               !
!                          -------------                               !
!                                                                      !
! integer function no_day(year,mnth,day,thty_dys)                      !
!                                                                      !
!----------------------------------------------------------------------!
integer function no_day(year,mnth,day,thty_dys)
!**********************************************************************!
integer :: year,mnth,day,i,thty_dys
!----------------------------------------------------------------------!

no_day = day
if (mnth>1) then
  do i=1,mnth-1
    no_day = no_day + no_days(year,i,thty_dys)
  enddo
endif

end function no_day





!**********************************************************************!
!                                                                      !
!                     st2arr :: func                                   !
!                     --------------                                   !
!                                                                      !
! Strips leading blanks and first integer from the string which is its !
! argument.                                                            !
!                                                                      !
! subroutine st2arr(st1,x,ix,no)                                       !
!                                                                      !
!----------------------------------------------------------------------!
subroutine st2arr(st1,x,ix,no)
!**********************************************************************!
integer :: x(max_years),ix,i,no
character(len=str_len) :: st1
!----------------------------------------------------------------------!

do i=1,max_years
  call STRIPBN(st1,x(i))
enddo

no = 0
do i=1,max_years
  if (x(i)>0)  no = no + 1
enddo

end subroutine st2arr





!**********************************************************************!
!                                                                      !
!                     stripbn :: func                                  !
!                     ---------------                                  !
!                                                                      !
! Strips leading blanks and first integer from the string which is its !
! argument.                                                            !
!                                                                      !
! subroutine stripbn(st1,ans)                                          !
!                                                                      !
!**********************************************************************!
subroutine stripbn(st1,ans)
!**********************************************************************!
character(len=str_len) :: st1
integer :: nb,i,ans,check
!----------------------------------------------------------------------!

nb = blanks(st1)
if (nb>0) then
  do i=1,str_len-nb
    st1(i:i) = st1(i+nb:i+nb)
  enddo
  do i=str_len-nb+1,str_len
    st1(i:i) = ' '
  enddo
endif

check = 0
nb = blank(st1)
ans = 0
do i=1,nb
  if ((ichar(st1(i:i))<48).or.(ichar(st1(i:i))>57)) check = 1
  ans = ans + 10**(nb-i)*(ichar(st1(i:i))-48)
enddo

if (nb>0) then
  do i=1,str_len-nb
    st1(i:i) = st1(i+nb:i+nb)
  enddo
  do i=str_len-nb+1,str_len
    st1(i:i) = ' '
  enddo
endif

if (check==1) ans = -9999

end subroutine stripbn





!**********************************************************************!
!                                                                      !
!                     stripbs :: func                                  !
!                     ---------------                                  !
!                                                                      !
! Strips leading blanks and first string from the string which is its  !
! argument.                                                            !
!                                                                      !
! subroutine stripbs(st1,st2)                                          !
!                                                                      !
!----------------------------------------------------------------------!
subroutine stripbs(st1,st2)
!**********************************************************************!
character(len=str_len) :: st1,st2
integer :: nb,i
!----------------------------------------------------------------------!

st2 = ''
nb = blanks(st1)
do i=1,str_len-nb
  st1(i:i) = st1(i+nb:i+nb)
enddo
do i=str_len-nb+1,str_len
  st1(i:i) = ' '
enddo

nb = blank(st1)
do i=1,nb
  st2(i:i) = st1(i:i)
enddo
st2(nb+1:nb+1) = ' '

if (nb>0) then
  do i=1,str_len-nb
    st1(i:i) = st1(i+nb:i+nb)
  enddo
  do i=str_len-nb+1,str_len
    st1(i:i) = ' '
  enddo
endif

end subroutine stripbs





!**********************************************************************!
!                                                                      !
!                     blank :: func                                    !
!                     -------------                                    !
!                                                                      !
! 'blank' returns the number of characters before the space in the     !
! string which is its argument.                                        !
!                                                                      !
! integer function blank(st1)                                          !
!                                                                      !
!----------------------------------------------------------------------!
integer function blank(st1)
!**********************************************************************!
character(len=str_len) :: st1
!----------------------------------------------------------------------!

blank = 0
10    continue
blank = blank + 1
if ((blank<=str_len).and.(ichar(st1(blank:blank))/=32).and.(ichar(st1(blank:blank))/=9))  goto 10
!Space and tabulation
blank = blank - 1

end function blank





!**********************************************************************!
!                                                                      !
!                      blanks :: func                                  !
!                     ---------------                                  !
!                                                                      !
! 'blanks' returns the number of leading blanks.                       !
!                                                                      !
! integer function blanks(st1)                                         !
!                                                                      !
!----------------------------------------------------------------------!
integer function blanks(st1)
!**********************************************************************!
character(len=str_len) :: st1
!----------------------------------------------------------------------!

blanks = 0
10    continue
blanks = blanks + 1
if (blanks<=str_len) then
  if ((ichar(st1(blanks:blanks))==32).or. &
 (ichar(st1(blanks:blanks))==9)) goto 10
endif
blanks = blanks - 1

end function blanks





!**********************************************************************!
!                                                                      !
!                     n_fields :: func                                 !
!                     ----------------                                 !
!                                                                      !
! Reurns the number of space delimeted fields.                         !
!                                                                      !
! integer function n_fields(st1)                                       !
!                                                                      !
!----------------------------------------------------------------------!
integer function n_fields(st1)
!**********************************************************************!
character(len=str_len) :: st1,st2,st3
!----------------------------------------------------------------------!

n_fields = 0
st2 = st1

10    continue
if (blanks(st2)<str_len) then
  call STRIPBS(st2,st3)
  n_fields = n_fields + 1
  goto 10
endif

end function n_fields





!**********************************************************************!
!                                                                      !
!                     in2st :: func                                    !
!                     -------------                                    !
!                                                                      !
! 'in2st' returns the string value of its integer argument.            !
!                                                                      !
! character(str_len) function in2st(int1)                              !
!                                                                      !
!----------------------------------------------------------------------!
character(str_len) function in2st(int1)
!**********************************************************************!
integer :: loc,int,int1
!----------------------------------------------------------------------!

int = int1

do loc=1,str_len
 in2st(loc:loc) = ' '
enddo

loc = str_len
30    in2st(loc:loc) = char(ichar('0') + mod(int,10))
int = int/10
loc = loc - 1
if (int>0)  goto 30

end function in2st





!**********************************************************************!
!                                                                      !
!                     ntags :: func                                    !
!                     -------------                                    !
!                                                                      !
! integer function ntags(st1,st2)                                      !
!                                                                      !
!----------------------------------------------------------------------!
integer function ntags(st1,st2)
!**********************************************************************!
character(len=str_len), dimension(max_outputs) :: st1
character(len=str_len) :: st2,st3
integer :: i,j
!----------------------------------------------------------------------!

i = 0
10    continue
  i = i + 1
  st3 = st1(i)
  j = stcmp(st2,st3)
if ((j==0).and.(i<max_pfts)) goto 10

if (j==1) then
  ntags = i
else
  ntags = -9999
endif

end function ntags





!**********************************************************************!
!                                                                      !
!                     stcmp :: func                                    !
!                     -------------                                    !
!                                                                      !
!  integer function stcmp(st1,st2)                                     !
!                                                                      !
!----------------------------------------------------------------------!
integer function stcmp(st1,st2)
!**********************************************************************!
character(len=str_len) :: st1,st2
integer :: k
!----------------------------------------------------------------------!

stcmp = 1
k = 0
10    continue
  k = k + 1
  if (((ichar(st1(k:k))==0).or.(ichar(st1(k:k))==32).or.(ichar(st1(k:k))==9)) &
 .and.((ichar(st2(k:k))==0).or.(ichar(st2(k:k))==32).or.(ichar(st2(k:k))==9))) then
    stcmp = 1
    goto 20
  endif
  if (st1(k:k)/=st2(k:k)) then
    stcmp = 0
    goto 20
  endif
goto 10 
20 continue

end function stcmp





!**********************************************************************!
!                                                                      !
!                     randomv :: func                                  !
!                     ---------------                                  !
!                                                                      !
! subroutine randomv(years,yr0,yrf,idum)                               !
!                                                                      !
!----------------------------------------------------------------------!
subroutine randomv(years,yr0,yrf,idum)
!**********************************************************************!
integer :: yr0,yrf,years(yrf-yr0+1),idum,ind(max_years),i,i1,j
!----------------------------------------------------------------------!

do i=1,yrf-yr0+1
  ind(i) = i
enddo

do i=1,yrf-yr0+1
  i1 = ranv(1,yrf-yr0+2-i,idum)
  years(i) = ind(i1)+yr0-1
  do j=i1,yrf-yr0+1-i
    ind(j) = ind(j+1)
  enddo
enddo

end subroutine randomv





!**********************************************************************!
!                                                                      !
!                     ranv :: func                                     !
!                     ------------                                     !
!                                                                      !
! integer function ranv(i,j,idum)                                      !
!                                                                      !
!----------------------------------------------------------------------!
integer function ranv(i,j,idum)
!**********************************************************************!
integer :: i,j,idum

ranv = int(ran1(idum)*(j-i+1))+i

end function ranv







!**********************************************************************!
!                                                                      !
!                          bi_lin :: data                              !
!                          --------------                              !
!                                                                      !
!                                                                      !
!  subroutine bi_lin(xx,indx,xnorm,ynorm,ans)                          !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Performs bilinear interpolation between four points.
!! @details Performs bilinear interpolation between four points, the
!! normalised
!! distances from the point (1,1) are given by 'xnorm' and 'ynorm'.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine bi_lin(xx,indx,xnorm,ynorm,ans)
!----------------------------------------------------------------------!
real(dp) :: xx(4,4),xnorm,ynorm,ans,av
integer :: indx(4,4),iav,ii,jj
!**********************************************************************!

!----------------------------------------------------------------------!
! Fill in averages if necessary.                                       !
!----------------------------------------------------------------------!
do ii=2,3
  do jj=2,3
    if (indx(ii,jj)/=1) then
      av = 0.0
      iav = 0
      if (indx(ii+1,jj)==1) then
        av = av + xx(ii+1,jj)
        iav = iav + 1
      endif
      if (indx(ii-1,jj)==1) then
        av = av + xx(ii-1,jj)
        iav = iav + 1
      endif
      if (indx(ii,jj+1)==1) then
        av = av + xx(ii,jj+1)
        iav = iav + 1
      endif
      if (indx(ii,jj-1)==1) then
        av = av + xx(ii,jj-1)
        iav = iav + 1
      endif
      if (indx(ii+1,jj+1)==1) then
        av = av + xx(ii+1,jj+1)
        iav = iav + 1
      endif
      if (indx(ii-1,jj-1)==1) then
        av = av + xx(ii-1,jj-1)
        iav = iav + 1
      endif
      if (indx(ii+1,jj-1)==1) then
        av = av + xx(ii+1,jj-1)
        iav = iav + 1
      endif
      if (indx(ii-1,jj+1)==1) then
        av = av + xx(ii-1,jj+1)
        iav = iav + 1
      endif
      if (iav>0) then
        xx(ii,jj) = av/real(iav)
      endif
    endif
  enddo
enddo

!----------------------------------------------------------------------!
! Bilinear interpolation.                                              !
!----------------------------------------------------------------------!
ans = xx(2,2)*(1.0-xnorm)*(1.0-ynorm) + xx(3,2)*xnorm*(1.0-ynorm) + &
 xx(2,3)*(1.0-xnorm)*ynorm + xx(3,3)*xnorm*ynorm

!----------------------------------------------------------------------!
! nearest pixel.                                                       !
!----------------------------------------------------------------------!
! ans = xx(int(xnorm+2.5),int(ynorm+2.5))
!----------------------------------------------------------------------!

end subroutine bi_lin





!**********************************************************************!
!                                                                      !
!                     ran1 :: func                                     !
!                     ------------                                     !
!                                                                      !
! real(dp) function ran1(idum)                                         !
!                                                                      !
!----------------------------------------------------------------------!
real(dp) function ran1(idum)
!**********************************************************************!
integer :: idum,IA,IM,IQ,IR,NTAB,NDIV
real(dp) :: AM,EPS,RNMX
parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
integer :: j,k,iv(NTAB),iy
save iv,iy
data iv /NTAB*0/, iy /0/
!----------------------------------------------------------------------!
 
if (idum<=0.or.iy==0) then
  idum = max(-idum,1)
  do j=NTAB+8,1,-1
    k = idum/IQ
    idum = IA*(idum-k*IQ)-IR*k
    if (idum<0)  idum = idum + IM
    if (j<=NTAB)  iv(j) = idum
  enddo
  iy = iv(1)
endif
k = idum/IQ
idum = IA*(idum - k*IQ) - IR*k
if (idum<0)  idum = idum + IM
j = 1 + iy/NDIV
iy = iv(j)
iv(j) = idum
ran1 = min(AM*iy,RNMX)

end function ran1



end module func
