!> @brief Collection of subroutines to read in data.
!! @details
!! @author Mark Lomas
!! @date July 2016

module data

use real_precision
use dims
use func
use input_file
use site_parameters
use system_state

implicit none

contains

!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!


!**********************************************************************!
!                                                                      !
!                          ex_clim :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine ex_clim(lat,lon,xlatf,xlatres,xlatresn,xlon0,     !
! xlonres,xlonresn,yr0,yrf,tmpv,humv,prcv,isite,year0,yearf,siteno,du) !
!                                                                      !                                                            !
!----------------------------------------------------------------------!
!> @brief Extract climate data from the climate database for the
!! nearest site.
!! @details ! Extract climate data from the climate database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database.
!!
!!            UNIX                DOS
!!
!!            ii(4)               ii(5)   for beginning of binary file
!!                                        records
!!            recl = 728          recl = 730     for binary climate
!!            recl = 577          recl = 578     for text map
!!
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_clim(lat,lon,xlatf,xlatres,xlatresn,xlon0, &
 xlonres,xlonresn,yr0,yrf,tmpv,humv,prcv,isite,year0,yearf,siteno,du, &
 read_par)
!**********************************************************************!
real(dp) :: lat,lon,xlon0,xlatf,xlatres,xlonres,ans(12)
real(dp), dimension(500,12,31) :: swrv
real(dp) :: tmpv(300,12,31),humv(300,12,31),prcv(300,12,31)
integer :: year,year0,yearf,nrec,ncol,ans2(1000),siteno,i,du
integer :: nyears,yr0,mnth,day,isite,yrf,recl1,recl2
integer :: xlatresn,xlonresn,fno,read_par
character :: ii(4),jj(5),num*3
character(len=str_len) :: fname1,fname2,fname3,fname4
!----------------------------------------------------------------------!

if (du==1) then
  recl1 = 730
  recl2 = 6*xlonresn
else
  recl1 = 728
  recl2 = 6*xlonresn + 1
endif

nyears = yearf - year0 + 1

nrec = int((xlatf - lat)/xlatres + 1.0)
ncol = int((lon - xlon0)/xlonres + 1.0)

if ((nrec<=xlatresn).and.(ncol<=xlonresn)) then

lat = xlatf - xlatres/2.0 - real(nrec - 1)*xlatres
lon = xlon0 + xlonres/2.0 + real(ncol - 1)*xlonres

fno = 90

open(fno+1,file=trim(inp%dirs%climate)//'/maskmap.dat', &
 access='DIRECT',recl=recl2,form='formatted',status='OLD')

read(fno+1,'(96i6)',rec=nrec) (ans2(i),i=1,xlonresn)
close(fno+1)
siteno = ans2(ncol)

isite = 0
if (siteno>0) then

  isite = 1

  write(num,'(i3.3)') (siteno-1)/100
  write(fname1,'(100a)') trim(inp%dirs%climate)//'/tmp_',num
  write(fname2,'(100a)') trim(inp%dirs%climate)//'/hum_',num
  write(fname3,'(100a)') trim(inp%dirs%climate)//'/prc_',num
  write(fname4,'(100a)') trim(inp%dirs%climate)//'/swr_',num

  siteno = mod(siteno-1,100) + 1

  open(fno+1,file=fname1,access='direct',recl=recl1,form='unformatted',status='old')
  open(fno+2,file=fname2,access='direct',recl=recl1,form='unformatted',status='old')
  open(fno+3,file=fname3,access='direct',recl=recl1,form='unformatted',status='old')
  open(fno+4,file=fname4,access='direct',recl=recl1,form='unformatted',status='old')

  if (du==1) then
    do year=yr0,yrf
      read(fno+1,rec=(siteno-1)*nyears+year-year0+1) jj, &
 ((tmpv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+2,rec=(siteno-1)*nyears+year-year0+1) jj, &
 ((humv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+3,rec=(siteno-1)*nyears+year-year0+1) jj, &
 ((prcv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      do mnth=1,12
        do day=1,30
          prcv(year-yr0+1,mnth,day) = int(real(prcv(year-yr0+1,mnth,day))/10.0 + 0.5)
        enddo
      enddo
    enddo
  else
    do year=yr0,yrf
      read(fno+1,rec=(siteno-1)*nyears+year-year0+1) ii, &
 ((tmpv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+2,rec=(siteno-1)*nyears+year-year0+1) ii, &
 ((humv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+3,rec=(siteno-1)*nyears+year-year0+1) ii, &
 ((prcv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      do mnth=1,12
        do day=1,30
          prcv(year-yr0+1,mnth,day) = int(real(prcv(year-yr0+1,mnth,day))/10.0 + 0.5)
        enddo
      enddo
    enddo
  endif

  close(fno+1)
  close(fno+2)
  close(fno+3)
  close(fno+4)

endif

else
  siteno = 0
endif

do mnth=1,12
  ans(mnth) = 0.0
  do day=1,30
    ans(mnth) = ans(mnth) + real(tmpv(1,mnth,day))/100.0
  enddo
enddo

end subroutine ex_clim





!**********************************************************************!
!                                                                      !
!                          ex_clim_site :: data                        !
!                          --------------------                        !
!                                                                      !
!                                                                      !
! subroutine ex_clim_site(yr0,yrf,tmpv,humv,prcv,year0,swrv,yearf)     !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Extract climate data from the climate database for the
!! nearest site.
!! @details Extract climate data from the climate database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database.
!!
!!            UNIX                DOS
!!
!!            ii(4)               ii(5)   for beginning of binary file
!!                                        records
!!            recl = 728          recl = 730     for binary climate
!!            recl = 577          recl = 578     for text map
!! @author Mark Lomas,EPK added swrv,l_clim,l_stats flags,300->500
!! @date Oct 2018
!----------------------------------------------------------------------!
subroutine ex_clim_site(yr0,yrf,tmpv,humv,prcv,swrv,year0,yearf, &
 l_clim,l_stats)
!----------------------------------------------------------------------!
real(dp) :: tmp,prc,hum,swr
real(dp) :: tmpv(500,12,31),humv(500,12,31),prcv(500,12,31),swrv(500,12,31)
integer :: year,year0,yearf,yr0,mnth,day,yrf,iyear,imnth,iday
logical :: l_clim,l_stats
!**********************************************************************!

open(91,file=trim(inp%dirs%climate)//'/site.dat')

do year=year0,yearf
  do mnth=1,12
    do day=1,no_days(year,mnth,0)
      read(91,*) iyear,imnth,iday,tmp,prc,hum,swr
      if ((iday/=day).or.(imnth/=mnth).or.(iyear/=year)) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'Error in climate data file',year,mnth,day
        stop
      endif
       
      if ((year<=yrf).and.(year>=yr0)) then
        tmpv(year-yr0+1,mnth,day) = tmp
        prcv(year-yr0+1,mnth,day) = prc
        humv(year-yr0+1,mnth,day) = hum
        swrv(year-yr0+1,mnth,day) = swr
      endif
    enddo
  enddo
enddo

!Flags that the climate has been read.Checked in sdgvm.f90
l_clim=.true.
l_stats=.true.
close(91)

end subroutine ex_clim_site





!**********************************************************************!
!                                                                      !
!                          ex_clu :: data                              !
!                          --------------                              !
!                                                                      !
!  subroutine EX_CLU(lat,lon,nft,lutab,cluse,yr0,yrf,du,l_lu)   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details  Extract land use from land use database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database. The data contains the percent (in byte format) for
!! each pft from 2 to nft.
!! First reads the readme.dat from the land use dataset with info such 
!! as resolution,years available and number of classes.Description of
!! classes and proportion assigned to model ft are ignored.
!! Land use must be written per class per year in vector format (i3)
!! 0-100 with 255 for water.Direction is West to East,North to South.
!! lutab(landcover classes,nft) is read from setup file.
!! classprop(landcover classes) The % of each landcover class in the grid
!! ftprop(nft) the fraction of each ft in the gridcell
!! cluse(year,nft) the fraction of each ft in the gridcell per year
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_clu(lat,lon,nft,lutab,cluse,yr0,yrf,du,l_lu)
!**********************************************************************!
real(dp) :: lat,lon,lon0,latf,latr,lonr,classprop(255)
real(dp) :: cluse(max_cohorts,max_years),lutab(255,100),ans
real(dp) :: ftprop(max_cohorts),rrow,rcol,xx(4,4),xnorm,ynorm
integer :: i,n,j,du,latn,lonn,row,col,recn,k,x,nft,ift
integer :: ii,jj,indx(4,4),yr0,yrf,years(1000),nrecl
integer :: classes(1000),nclasses,kode
character(len=str_len) :: st1,st2,st3
logical :: l_lu
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! read the 'readme.dat' file in the directory where land use files are !
!----------------------------------------------------------------------!
open(99,file=trim(inp%dirs%land_use)//'/readme.dat',&
 status='old',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Land use file does not exist.'
  write(*,'('' "'',A,''/readme.dat"'')') trim(inp%dirs%land_use)
  stop
endif

read(99,*) st1
st2='CONTINUOUS'
if (stcmp(st1,st2)==0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'landuse is not a continuous field ?'
  write(*,*) 'readme.dat should begin with CONTINUOUS'
  stop
endif
read(99,*)
read(99,*) latf,lon0
read(99,*)
read(99,*) latr,lonr
read(99,*)
read(99,*) latn,lonn
read(99,*)
read(99,'(A)') st1
n = n_fields(st1)
call ST2ARR(st1,years,1000,n)
read(99,*)
read(99,'(A)') st1
close(99)
nclasses = n_fields(st1)
call ST2ARR(st1,classes,1000,nclasses)

!----------------------------------------------------------------------!
if (du==1) then
!  This works for ftn95
!  nrecl = 3
  nrecl = 5
else
  nrecl = 4
endif

if ((n>1).and.(yr0<years(1))) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Can''t start running in ',yr0,&
 ' since landuse map begin in ',years(1)
  stop
endif

!     look for the first year
j=1

 10   continue
if ((j<n).and.(years(j)<yr0)) then
   j = j + 1
   goto 10
endif

!----------------------------------------------------------------------!
! Find the real(dp) :: row col corresponding to lat and lon.                  !
!----------------------------------------------------------------------!
rrow = (latf - lat)/latr
rcol = (lon - lon0)/lonr

ynorm = rrow - real(int(rrow))
xnorm = rcol - real(int(rcol))
!----------------------------------------------------------------------!

j=1
do i=1,yrf-yr0+1
  if ((i==1).or.((i+yr0-1)==years(j))) then
    st2 = in2st(years(j))
    st2 = adjustl(st2)
    j = j + 1

    do k=1,nclasses
      classprop(classes(k)) = 0

      st3 = in2st(classes(k))
      st3 = adjustl(st3)
      open(99,file= &
 trim(inp%dirs%land_use)//'/cont_lu-'//trim(st3)//'-'//st2(1:4)//'.dat', & 
 status='old',form='formatted',access='direct',recl=nrecl,iostat=kode)
      if (kode/=0) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'Land Use data-base.'
        write(*,*) 'File does not exist:'
        write(*,*) trim(inp%dirs%land_use)//'/cont_lu-',trim(st3),'-',st2(1:4),'.dat'
        stop
      endif

      do ii=1,4
        do jj=1,4
          row = int(rrow)+jj-1
          col = int(rcol)+ii-1
          if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) then
            recn = (row-1)*lonn + col
            read(99,'(i3)',rec=recn) x 
            xx(ii,jj) = real(x)
            if (x<200) then
              indx(ii,jj) = 1
            else
              indx(ii,jj) = 0
            endif
          else
            indx(ii,jj) = -1
          endif
        enddo
      enddo

      call bi_lin(xx,indx,xnorm,ynorm,ans)
      
      x = int(ans+0.5)
      
      !This allows only specific crops
      !IF (k/=29) ans=0.
      
      classprop(classes(k)) = ans
      close(99)

    enddo ! end of loop over the classes
    
!----------------------------------------------------------------------!
! Now calculate the ftprop.
!----------------------------------------------------------------------!
    do ift=2,nft
      ftprop(ift)=0.0
      do k=1,nclasses
        ftprop(ift)=ftprop(ift)+lutab(classes(k),ift)*classprop(classes(k))/100.0
      enddo
    enddo

!----------------------------------------------------------------------!
! Calculate the bare soil.
!----------------------------------------------------------------------!
    if ((ftprop(2)<=100).and.(ftprop(2)>=0)) then
      ftprop(1)=100
      do ift=2,nft
        ftprop(1)=ftprop(1)-ftprop(ift)
      enddo
    endif

  endif ! finished reading

  do ift=1,nft
    cluse(ift,i) = ftprop(ift)
  enddo

enddo

!! TRENDY SDGVM method - assumes smoothed change in land-use between years specified in input dataset
!DO i=1,years(n)-yr0a
!  IF ( (i==1).or.((i+yr0a-1)==years(j1)) ) THEN
!    ij  = i
!    ij1 = ij + years(j1+1) - years(j1) 
!    j1  = j1 + 1
!  ELSE
!    DO ift=1,nft
!      cluse(ift,i) = cluse(ift,ij) + ( (real(i)-real(ij))/(real(ij1)-real(ij))* & 
! (cluse(ift,ij1) - cluse(ift,ij)) )
!    ENDDO
!  ENDIF
!ENDDO

if ((indx(2,2)==1).or.(indx(2,3)==1).or.(indx(3,2)==1).or.(indx(3,3)==1)) then
  l_lu = .true.
else
  l_lu = .false.
endif

end subroutine ex_clu



!**********************************************************************!
!                                                                      !
!                          readCO2 :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine readCO2(yr0,yrf,co2)                                      !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Reads co2 values for each year
!! @details It will read from the co2 file the concentration vector 
!! with starting (yr0) and end year (yrf) of the run as produced by the 
!! read_input_file sub.
!! If the co2 file has only one value then it will be used for all the
!! years meaning constant co2 concentrations.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine readco2(yr0,yrf,co2,daily_co2,nyears)
!----------------------------------------------------------------------!
real(dp) :: co2(max_years),ca
integer :: yr0,yrf,norecs,year,const,kode
integer :: yra,yrfa,prev_year,nyears
integer :: mnth,yr0a,daily_co2,day
logical :: co2spin
!**********************************************************************!

open(98,file=trim(inp%dirs%co2),status='OLD',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Co2 file does not exist.'
  write(*,'('' "'',A,''"'')') trim(inp%dirs%co2)
  stop
endif

norecs = 0
const = 0
10    continue
  read(98,*,end=99) year,ca
    if ((year>=yr0).and.(year<=yrf)) then
    co2(norecs+1) = ca
    norecs = norecs + 1
  endif
  const = const + 1
goto 10
99 continue
close(98)

!If it has found a single value in the co2 file then use it for all years
if (const==1) then
  do year=yr0,yrf
    co2(year-yr0+1) = ca
  enddo
endif

end subroutine readco2




!**********************************************************************!
!                                                                      !
!                         landuse1 :: data                             !
!                         ----------------                             !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine landuse1(luse,yr0,yrf,fire_ant,harvest_ant)
!**********************************************************************!
integer :: yr0,yrf,year,early,rep,luse(max_years),i,use
logical :: fire_ant(max_years),harvest_ant(max_years)
!----------------------------------------------------------------------!

early = yrf+max_years

do i=1,max_years
  luse(i) = 0
enddo

do i=1,inp%land_use%n
  year = inp%land_use%year(i)
  if (year-yr0+1>0)  luse(year-yr0+1) = int(inp%land_use%map(i))
  if (year<early) rep = int(inp%land_use%map(i))
enddo

do i=1,max_years
  if (luse(i)>0) rep = luse(i)
  luse(i) = rep
enddo

return



!do i=1,max_years
!  luse(i) = 0
!enddo

!i = 1
!10    continue
!  read(fid_input,*,end=20) year,use
!  print*,use
!  if (year-yr0+1>0)  luse(year-yr0+1) = use
!  if (year<early) rep = use
!  i = i+1
!goto 10
!20    continue

!do i=1,max_years
!  if (luse(i)>0) rep = luse(i)
!  luse(i) = rep
!enddo


!return

! Anthony's updated version this seems to cause a problem.
!      early = yrf+max_years
!      luse(:) = 1000
      
!      i = 1
!11    CONTINUE
!        READ(98,*,end=21) year,use
!        IF (year-yr0+1>0)  luse(year-yr0+1) = use
!        IF ((year<early).and.(use>0)) rep = use
!        !early = year
!        i = i+1
!      GOTO 11
!21    CONTINUE

!      if(rep<=0) then
!        print*, 'ERROR:: Land use mapping equal to or below 0'
!        print*, 'you must specifiy at least one year with a cover type'
!      endif

!      DO i=1,max_years
!        If ((luse(i)>0).and.(luse(i)<1000)) rep = luse(i)
!        If (luse(i)<0) fire_ant(i)    = .TRUE.
!        If (luse(i)==0) harvest_ant(i) = .TRUE.
!        luse(i) = rep
!      ENDDO

end subroutine landuse1




!**********************************************************************!
!                                                                      !
!                          ex_lu :: data                               !
!                          -------------                               !
!                                                                      !
! subroutine ex_lu(lat,lon,luse,yr0,yrf,du)                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Extract land use.
!! @details Extract land use from land use soils database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_lu(lat,lon,luse,yr0,yrf,du)
!**********************************************************************!
real(dp) :: lat,lon,lon0,latf,latr,lonr,xlat,xlon
integer :: i,n,j,du,latn,lonn,row,col,recn,kode
integer :: luse(max_years),yr0,yrf,rep,years(1000),lu(1000),nrecl
character(len=str_len) :: st1
!**********************************************************************!

open(99,file=trim(inp%dirs%land_use)//'/readme.dat')
read(99,*)
read(99,*) latf,lon0
read(99,*)
read(99,*) latr,lonr
read(99,*)
read(99,*) latn,lonn
read(99,*)
read(99,'(A)') st1
close(99)

n = n_fields(st1)
call ST2ARR(st1,years,1000,n)

if (du==1) then
  nrecl = 16+n*3
else
  nrecl = 16+n*3+1
endif

lu(1) = 0
open(99,file=trim(inp%dirs%land_use)//'/landuse.dat',status='old', &
 form='formatted',access='direct',recl=nrecl,iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'File does not exist:',trim(inp%dirs%land_use)//'/landuse.dat'
  stop
endif

row = int((latf - lat)/latr + 1.0)
col = int((lon - lon0)/lonr + 1.0)

recn = (row-1)*lonn + col
if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=recn) xlat,xlon,(lu(i),i=1,n)

if (lu(1)==0) then
  row = int((latf - (lat + latr))/latr + 1.0)
  col = int(((lon) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat - latr))/latr + 1.0)
  col = int(((lon) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat))/latr + 1.0)
  col = int(((lon + lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat))/latr + 1.0)
  col = int(((lon - lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat + latr))/latr + 1.0)
  col = int(((lon + lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat - latr))/latr + 1.0)
  col = int(((lon + lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat + latr))/latr + 1.0)
  col = int(((lon - lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat - latr))/latr + 1.0)
  col = int(((lon - lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

close(99)

rep = lu(1)
j = 1
do i=1,n
  if (yr0>=years(j+1)) then
    rep = lu(i)
    j = i
  endif
enddo

do i=1,yrf-yr0+1
  if ((i+yr0-1>=years(j+1)).and.(j<n)) then
    j = j+1
    rep = lu(j)
  endif
  luse(i) = rep
enddo

end subroutine ex_lu



!**********************************************************************!
!                                                                      !
!                          ex_soil :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine ex_soil(lat,lon,sol_chr2,du,l_soil)                !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Extract % sand % silt bulk density and depth from soils
!! database.
!! @details Extract % sand % silt bulk density and depth from soils
!! database.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_soil(lat,lon,sol_chr2,du,l_soil)
!----------------------------------------------------------------------!
real(dp) :: lat,lon,lon0,latf,latr,lonr,xlat,xlon,sol_chr2(10)
integer :: row,col,recn,recl1,du,latn,lonn,i,ii,jj,kode
integer, dimension(4,4) :: indx1,indx2,indx3,indx4,indx5
integer, dimension(4,4) :: indx6,indx7,indx8,indx9,indx10
real(dp),dimension(4,4) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8,xx9,xx10
real(dp) :: ynorm,xnorm,rrow,rcol
real(dp) :: ans
logical ::  l_soil(20)

if (du==1) then
! This works for ftn95
!  recl1 = 16+8*10
  recl1 = 16+10*10
else
  recl1 = 16+10*10+1
endif

open(99,file=trim(inp%dirs%soil)//'/readme.dat',status='old',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Soil data file does not exist.'
  write(*,'('' "'',A,''/readme.dat"'')') trim(inp%dirs%soil)
  stop
endif

read(99,*)
read(99,*) latf,lon0
read(99,*)
read(99,*) latr,lonr
read(99,*)
read(99,*) latn,lonn
close(99)

!----------------------------------------------------------------------!
! Find the real(dp) :: row col corresponding to lat and lon.                  !
!----------------------------------------------------------------------!
rrow = (latf - lat)/latr
rcol = (lon - lon0)/lonr

ynorm = rrow - real(int(rrow))
xnorm = rcol - real(int(rcol))
!----------------------------------------------------------------------!

open(99,file=trim(inp%dirs%soil)//'/data.dat',status='old', &
 form='formatted',access='direct',recl=recl1,iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Soil data-base.'
  write(*,*) 'File does not exist:',trim(inp%dirs%soil),'/data.dat'
  write(*,*) 'Or record length missmatch.'
  stop
endif

do ii=1,4
  do jj=1,4
    row = int(rrow)+jj-1
    col = int(rcol)+ii-1
    if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) then
      recn = (row-1)*lonn + col
      read(99,'(f7.3,f9.3,10f10.4)',rec=recn) xlat,xlon,(sol_chr2(i),i=1,10)
      xx1(ii,jj) = sol_chr2(1)
      xx2(ii,jj) = sol_chr2(2)
      xx3(ii,jj) = sol_chr2(3)
      xx4(ii,jj) = sol_chr2(4)
      xx5(ii,jj) = sol_chr2(5)
      xx6(ii,jj) = sol_chr2(6)
      xx7(ii,jj) = sol_chr2(7)
      xx8(ii,jj) = sol_chr2(8)
      xx9(ii,jj) = sol_chr2(9)
      xx10(ii,jj) = sol_chr2(10)
      if (sol_chr2(1)<0.0) then
        indx1(ii,jj) = 0
      else
        indx1(ii,jj) = 1
      endif
      if (sol_chr2(2)<0.0) then
        indx2(ii,jj) = 0
      else
        indx2(ii,jj) = 1
      endif
      if (sol_chr2(3)<0.0) then
        indx3(ii,jj) = 0
      else
        indx3(ii,jj) = 1
      endif
      if (sol_chr2(4)<0.0) then
        indx4(ii,jj) = 0
      else
        indx4(ii,jj) = 1
      endif
      if (sol_chr2(5)<0.0) then
        indx5(ii,jj) = 0
      else
        indx5(ii,jj) = 1
      endif
      if (sol_chr2(6)<0.0) then
        indx6(ii,jj) = 0
      else
        indx6(ii,jj) = 1
      endif
      if (sol_chr2(7)<0.0) then
        indx7(ii,jj) = 0
      else
        indx7(ii,jj) = 1
      endif
      if (sol_chr2(8)<0.0) then
        indx8(ii,jj) = 0
      else
        indx8(ii,jj) = 1
      endif
      if (sol_chr2(9)<0.0d0) then
        indx9(ii,jj) = 0
      ELSE
        indx9(ii,jj) = 1
      endIF
      if (sol_chr2(10)<0.0d0) then
        indx10(ii,jj) = 0
      ELSE
        indx10(ii,jj) = 1
      endIF
    else
      indx1(ii,jj) = -1
      indx2(ii,jj) = -1
      indx3(ii,jj) = -1
      indx4(ii,jj) = -1
      indx5(ii,jj) = -1
      indx6(ii,jj) = -1
      indx7(ii,jj) = -1
      indx8(ii,jj) = -1
      indx9(ii,jj) = -1
      indx10(ii,jj) = -1
    endif
  enddo
enddo

close(99)

if ((indx1(2,2)==1).or.(indx1(2,3)==1).or.(indx1(3,2)==1).or. &
 (indx1(3,3)==1)) then
  call bi_lin(xx1,indx1,xnorm,ynorm,ans)
  sol_chr2(1) = ans
  l_soil(1) = .true.
else
  l_soil(1) = .false.
endif
if ((indx2(2,2)==1).or.(indx2(2,3)==1).or.(indx2(3,2)==1).or. &
 (indx2(3,3)==1)) then
  call bi_lin(xx2,indx2,xnorm,ynorm,ans)
  sol_chr2(2) = ans
  l_soil(2) = .true.
else
  l_soil(2) = .false.
endif
if ((indx3(2,2)==1).or.(indx3(2,3)==1).or.(indx3(3,2)==1).or. &
 (indx3(3,3)==1)) then
  call bi_lin(xx3,indx3,xnorm,ynorm,ans)
  sol_chr2(3) = ans
  l_soil(3) = .true.
else
  l_soil(3) = .false.
endif
if ((indx4(2,2)==1).or.(indx4(2,3)==1).or.(indx4(3,2)==1).or. &
 (indx4(3,3)==1)) then
  call bi_lin(xx4,indx4,xnorm,ynorm,ans)
  sol_chr2(4) = ans
  l_soil(4) = .true.
else
  l_soil(4) = .false.
endif
if ((indx5(2,2)==1).or.(indx5(2,3)==1).or.(indx5(3,2)==1).or. &
 (indx5(3,3)==1)) then
  call bi_lin(xx5,indx5,xnorm,ynorm,ans)
  sol_chr2(5) = ans
  l_soil(5) = .true.
else
  l_soil(5) = .false.
endif
if ((indx6(2,2)==1).or.(indx6(2,3)==1).or.(indx6(3,2)==1).or. &
 (indx6(3,3)==1)) then
  call bi_lin(xx6,indx6,xnorm,ynorm,ans)
  sol_chr2(6) = ans
  l_soil(6) = .true.
else
  l_soil(6) = .false.
endif
if ((indx7(2,2)==1).or.(indx7(2,3)==1).or.(indx7(3,2)==1).or. &
 (indx7(3,3)==1)) then
  call bi_lin(xx7,indx7,xnorm,ynorm,ans)
  sol_chr2(7) = ans
  l_soil(7) = .true.
else
  l_soil(7) = .false.
endif
if ((indx8(2,2)==1).or.(indx8(2,3)==1).or.(indx8(3,2)==1).or. &
 (indx8(3,3)==1)) then
  call bi_lin(xx8,indx8,xnorm,ynorm,ans)
  sol_chr2(8) = ans
  l_soil(8) = .true.
else
  l_soil(8) = .false.
endif
if ((indx9(2,2)==1).or.(indx9(2,3)==1).or.(indx9(3,2)==1).or. &
 (indx9(3,3)==1)) THEN
  CALL BI_LIN(xx9,indx9,xnorm,ynorm,ans)
  sol_chr2(9) = ans
  l_soil(9) = .true.
ELSE
  l_soil(9) = .false.
endIF
if ((indx10(2,2)==1).or.(indx10(2,3)==1).or.(indx10(3,2)==1).or. &
 (indx10(3,3)==1)) THEN
  CALL BI_LIN(xx10,indx10,xnorm,ynorm,ans)
  sol_chr2(10) = ans
  l_soil(10) = .true.
ELSE
  l_soil(10) = .false.
endIF

end subroutine ex_soil




!**********************************************************************!
!                                                                      !
!                          read_soil :: sdgvm1                         !
!                          -------------------                         !
!                                                                      !
! subroutine read_soil(lat,lon,soil_chr,soil_chr2,du,l_soil)           !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_soil(lat,lon,soil_chr,soil_chr2,du,l_soil)
!**********************************************************************!
real(dp) :: lat,lon,soil_chr(10),soil_chr2(10)
integer :: du
logical :: l_soil(20)
!----------------------------------------------------------------------!

call EX_SOIL(lat,lon,soil_chr2,du,l_soil)
if (soil_chr(1)>0.01) then
  ssp%sand = soil_chr(1)
  ssp%silt = soil_chr(2)
  l_soil(1) = .true.
  l_soil(2) = .true.
else
  ssp%sand = soil_chr2(1)
  ssp%silt = soil_chr2(2)
endif
ssp%clay = 100.0 - ssp%sand - ssp%silt

if (soil_chr(3)>0.01) then
  ssp%bulk  = soil_chr(3)
  l_soil(3) = .true.
else
  ssp%bulk = soil_chr2(3)
endif

if (soil_chr(4)>0.01) then
  ssp%orgc  = soil_chr(4)
  l_soil(4) = .true.
else
  ssp%orgc = soil_chr2(4)
endif

if (soil_chr(5)>0.01) then
  ssp%wilt  = soil_chr(5)
  ssp%field = soil_chr(6)
  ssp%sat   = soil_chr(7)
  l_soil(5) = .true.
  l_soil(6) = .true.
  l_soil(7) = .true.
else
  ssp%wilt = soil_chr2(5)
  ssp%field = soil_chr2(6)
  ssp%sat = soil_chr2(7)
endif

if (soil_chr(8)>0.01) then
  ssp%soil_depth = soil_chr(8)
  l_soil(8) = .true.
else
  ssp%soil_depth = soil_chr2(8)
endif

end subroutine read_soil


!**********************************************************************!
!                                                                      !
!                          co2_0_f :: sdgvm1                           !
!                          -----------------                           !
!                                                                      !
! subroutine co2_0_f(co20,co2f,yearv,yr0,co2,nyears)                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Extract initial and final C02 values
!! 
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine co2_0_f(co20,co2f,yearv,yr0,co2,nyears)
!**********************************************************************!
integer :: iyear,year,nyears,yr0,yearv(max_years)
real(dp) :: co20,co2f,co2(max_years)
!----------------------------------------------------------------------!

iyear = 1
year = yearv(iyear)
if ((inp%run%spinup_length>0).and.(iyear>inp%run%spinup_length)) then
  co20 = co2(year-yr0+1)
else
  if (inp%run%co2_constant>0.0) then
    co20 = inp%run%co2_constant
  else
    co20 = co2(year-yr0+1)
  endif
endif

iyear = nyears
year = yearv(iyear)
if ((inp%run%spinup_length>0).and.(iyear>inp%run%spinup_length)) then
  co2f = co2(year-yr0+1)
else
  if (inp%run%co2_constant>0.0) then
    co2f = inp%run%co2_constant
  else
    co2f = co2(year-yr0+1)
  endif
endif

end subroutine co2_0_f


!**********************************************************************!
!                                                                      !
!                          set_landuse :: sdgvm1                       !
!                          ---------------------                       !
!                                                                      !
! subroutine set_landuse(ftprop,ilanduse,tmp,prc,nat_map,nft,cluse,    !
! year,yr0)                                                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Set landuse by allocating the cover of new ground.
!! @details Assign gound available for new growth, through disturbance
!! and mortality, given by ftprop. 
!! For the specific year,check whether each ft can grow.
!! If it can then ftprop(ft) acquires the value of cluse array which holds
!! the desired cover for each ft and year as read from cover file.
!! If ftprop(1)<0 then it sets it to 0 and proportionally reduced the 
!! cover of the ofther fts to add up to 100.
!!
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_landuse(ftprop,tmp,prc,nat_map,nft,cluse,year,yr0)
!**********************************************************************!
integer :: ft,nft,nat_map(8)
integer :: year,yr0
real(dp) :: tmp(12,31),prc(12,31),cluse(max_cohorts,max_years)
real(dp) :: ftprop(max_cohorts)
!----------------------------------------------------------------------!

ftprop(1) = 100.0
do ft=2,nft
  if (check_ft_grow(tmp,ssv(1)%chill,ssv(1)%dschill,ft)==1) then
    ftprop(ft) = cluse(ft,year-yr0+1)
    ftprop(1) = ftprop(1) - ftprop(ft)
  else
    ftprop(ft) = 0.0
  endif
enddo
if (ftprop(1)<0.0) then
  do ft=2,nft
    ftprop(ft) = ftprop(ft)*100.0/(100.0 - ftprop(1))
  enddo
  ftprop(1) = 0.0
endif

end subroutine set_landuse



!**********************************************************************!
!                                                                      !
!                          set_climate :: sdgvm1                       !
!                          ---------------------                       !
!                                                                      !
! subroutine set_climate(xtmpv,xprcv,xhumv,xcldv,withcloudcover,yearv, !
! iyear,tmp,prc,hum,cld,thty_dys,yr0,year)                             !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Set 'tmp' 'hum' 'prc' 'swr' and 'cld'.
!! @details Assign climate from the raw data which is held in 'x???v'
!! arrays.
!! It also calculates the average temperature,precipitation for a month
!! and humidity.It uses that info to calculate the exponentially-weighted
!! 20-year monthly means which are assigned to site structure
!! variable ssp.Those values are required for the crop processes.
!! Note that here I want to calculate the 20-year monthly mean twice,one
!! for this year and one for the next which is the reason why this sub
!! is in a loop.Careful,the index nn1 that goes from the loop in the sub
!! goes from 1 to 0. 
!! @author Mark Lomas,EPK
!! @date Oct 2016
!----------------------------------------------------------------------!
subroutine set_climate(xtmpv,xprcv,xhumv,xcldv,xswrv,withcloudcover,yearv,&
 iyear,tmp,prc,hum,cld,swr,thty_dys,yr0,year,nyears,nn1)
!**********************************************************************!
real(dp), dimension(500,12,31) :: xtmpv,xprcv,xhumv,xswrv
real(dp), dimension(500,12) :: xcldv
integer :: yearv(max_years),iyear,mnth,day,thty_dys,yr0,year,nn1,nyears,nn2
logical :: withcloudcover
real(dp) :: tmp(12,31),prc(12,31),hum(12,31),cld(12),swr(12,31)
!----------------------------------------------------------------------!

!nn2 plays no role unless we are in the last year of the run where it
!ensures that we wont be reading outside the array.
nn2=0
if(iyear==nyears.and.nn1==1) nn2=-1

ssp%mnthtmp(:)=0.
ssp%mnthprc(:)=0.
ssp%mnthhum(:)=0.

do mnth=1,12
  do day=1,no_days(year,mnth,thty_dys)
    tmp(mnth,day) = real(xtmpv(yearv(iyear+nn1+nn2)-yr0+1,mnth,day))/100.0
    prc(mnth,day) = real(xprcv(yearv(iyear+nn1+nn2)-yr0+1,mnth,day))/10.0
    hum(mnth,day) = real(xhumv(yearv(iyear+nn1+nn2)-yr0+1,mnth,day))/100.0
    swr(mnth,day) = real(xswrv(yearv(iyear+nn1+nn2)-yr0+1,mnth,day))
    !Mean month temperature
    ssp%mnthtmp(mnth) = ssp%mnthtmp(mnth)+tmp(mnth,day)/no_days(year+nn1+nn2,mnth,thty_dys)
    !Mean month precipitation
    ssp%mnthprc(mnth) = ssp%mnthprc(mnth)+prc(mnth,day)
    if (withcloudcover) then
      cld(mnth) = real(xcldv(yearv(iyear+nn1+nn2)-yr0+1,mnth))/1000.0
    else
      cld(mnth) = 0.5
    endif
  enddo
enddo

do mnth=1,12
  do day=1,no_days(year+nn1+nn2,mnth,thty_dys)
    if (hum(mnth,day)<30.0)  hum(mnth,day) = 30.0
    if (hum(mnth,day)>95.0)  hum(mnth,day) = 95.0
    !Mean month humidity
    ssp%mnthhum(mnth) = ssp%mnthhum(mnth)+hum(mnth,day)/no_days(year+nn1+nn2,mnth,thty_dys)
  enddo
enddo


IF(iyear==1) THEN
  ssp%emnthtmp(:,nn1+1) = ssp%mnthtmp(:)
  ssp%emnthprc(:,nn1+1) = ssp%mnthprc(:)
  ssp%emnthhum(:,nn1+1) = ssp%mnthhum(:)
ELSE
  !Running averages for this year ssp%*(1:12,1) and the next ssp%*(1:12,2)
  ssp%emnthtmp(:,nn1+1) = 0.95*ssp%emnthtmp(:,nn1+1)+0.05*ssp%mnthtmp(:)
  ssp%emnthprc(:,nn1+1) = 0.95*ssp%emnthprc(:,nn1+1)+0.05*ssp%mnthprc(:)
  ssp%emnthhum(:,nn1+1) = 0.95*ssp%emnthhum(:,nn1+1)+0.05*ssp%mnthhum(:)
ENDIF

end subroutine set_climate





!**********************************************************************!
!                                                                      !
!                          set_co2 :: data                             ! 
!                          -----------------                           !
!                                                                      !
! subroutine set_co2(ca,iyear,speedc,co2,year,yr0)                     !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Set CO2 value 'ca'.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_co2(ca,iyear,speedc,co2,year,yr0)
!**********************************************************************!
real(dp) :: ca,co2(max_years)
integer :: iyear,year,yr0
logical :: speedc
!----------------------------------------------------------------------!

if ((inp%run%spinup_length>0).and.(iyear>inp%run%spinup_length)) then
  speedc = .false.
  ca = co2(year-yr0+1)
else
  if (inp%run%co2_constant>0.0) then
    ca = inp%run%co2_constant
  else
    ca = co2(year-yr0+1)
  endif
endif

end subroutine set_co2





!**********************************************************************!
!                                                                      !
!                        read_landuse :: sdgvm1                        !
!                        ----------------------                        !
!                                                                      !
! subroutine read_landuse(ilanduse,yr0,yrf,du,nft,lat,lon,lutab,       !
! luse,cluse,l_lu)                                                     !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details ilanduse:0 reads from map file,1 from the input file
!! 2 natural vegetation ;)
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_landuse(ilanduse,yr0,yrf,du,nft,lat,lon,lutab,&
 luse,cluse,l_lu)
!**********************************************************************!
integer :: ilanduse,icontinuouslanduse,yr0,yrf,du,year,ft,nft,&
 luse(max_years),nat_map(8)
real(dp) :: lat,lon,cluse(max_cohorts,max_years),lutab(255,100),sum
character(len=str_len) :: st2
character(len=str_len), dimension(max_outputs) :: fttags

logical :: l_lu
!----------------------------------------------------------------------!

icontinuouslanduse = 1

if (ilanduse==0) then
  if (icontinuouslanduse==0) then
    call EX_LU(lat,lon,luse,yr0,yrf,du)
    ! Create the continuous land use (cluse)
    do year=yr0,yrf
      do ft=1,nft
        cluse(ft,year-yr0+1) = lutab(luse(year-yr0+1),ft)
      enddo
    enddo
  else
    ! in data.f90
    ! Reads landuse from data map
    call EX_CLU(lat,lon,nft,lutab,cluse,yr0,yrf,du,l_lu)
  endif
elseif (ilanduse==1) then
  ! Reads landuse from input file
  l_lu = .true.
  do year=yr0,yrf
    do ft=1,nft
      cluse(ft,year-yr0+1) = lutab(luse(year-yr0+1),ft)
    enddo
  enddo
elseif (ilanduse==2) then
  write(*,*) 'Checking natural vegetation types exist:'
  write(*,*) 'BARE CITY C3 C4 Ev_Bl Ev_Nl Dc_Bl Dc_Nl.'
  l_lu = .true.
  st2 = 'BARE'
  nat_map(1) = ntags(fttags,st2)
  st2 = 'CITY'
  nat_map(2) = ntags(fttags,st2)
  st2 = 'C3'
  nat_map(3) = ntags(fttags,st2)
  st2 = 'C4'
  nat_map(4) = ntags(fttags,st2)
  st2 = 'Ev_Bl'
  nat_map(5) = ntags(fttags,st2)
  st2 = 'Ev_Nl'
  nat_map(6) = ntags(fttags,st2)
  st2 = 'Dc_Bl'
  nat_map(7) = ntags(fttags,st2)
  st2 = 'Dc_Nl'
  nat_map(8) = ntags(fttags,st2)
  do year=yr0,yrf
    do ft=1,nft
      cluse(ft,year-yr0+1) = 0.0
    enddo
  enddo
  cluse(1,1) = 100.0
endif

end subroutine read_landuse





!**********************************************************************!
!                                                                      !
!                          country :: sdgvm1                           !
!                          -----------------                           !
!                                                                      !
! subroutine country(lat,lon,country_name,country_id,l_regional)       !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Gets country id for the gridcell
!! 
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine country(lat,lon,country_name,country_id,l_regional)
!**********************************************************************!
real(dp) :: lat,lon,adj_lat,adj_lon,xlat,xlon
integer :: country_id,ilat,ilon,ans(360),x,i
character(len=str_len) :: country_name
logical :: l_regional
!----------------------------------------------------------------------!

ilat = int(91.0 - lat)
ilon = int(lon + 181.0)

open(99,file=trim(inp%dirs%land_mask)//'/country.dat',status='old')
  if (ilat>1) then
    do i=1,ilat-1
      read(99,*)
    enddo
  endif
  read(99,*) ans
  country_id = ans(ilon)
close(99)

!----------------------------------------------------------------------!
! If OCEAN was found try the nearest adjacent squares.                 !
!----------------------------------------------------------------------!
if (country_id==0) then

! Nearest lateral.

  xlat = lat+500-real(int(lat+500.0))
  if (xlat>0.5) then
    adj_lat = 1.0
  else
    adj_lat = -1.0
  endif

  xlon = lon+500-real(int(lon+500.0))
  if (xlon>0.5) then
    adj_lon = 1.0
  else
    adj_lon = -1.0
  endif

  if (abs(xlat-0.5)>abs(xlon-0.5)) then
    adj_lon = 0.0
  else
    adj_lat = 0.0
  endif

  ilat = int(91.0 - (lat + adj_lat))
  ilon = int((lon + adj_lat) + 181.0)

  if (ilon==0) ilon = 360
  if (ilon==361) ilon = 1
  if (ilat==0) ilat = 1
  if (ilat==91) ilat = 90

  open(99,file=trim(inp%dirs%land_mask)//'/country.dat',status='old')
    if (ilat>1) then
      do i=1,ilat-1
        read(99,*)
      enddo
    endif
    read(99,*) ans
    country_id = ans(ilon)
  close(99)

endif

if (country_id==0) then

! Next nearest lateral.

  if (abs(adj_lat)<0.5) then
    xlat = lat+500-real(int(lat+500.0))
    if (xlat>0.5) then
      adj_lat = 1.0
    else
      adj_lat = -1.0
    endif
  endif

  if (abs(adj_lon)<0.5) then
    xlon = lon+500-real(int(lon+500.0))
    if (xlon>0.5) then
      adj_lon = 1.0
    else
      adj_lon = -1.0
    endif
    adj_lat = 0.0
  endif

  if (abs(adj_lat)>0.5)  adj_lon = 0.0

  ilat = int(91.0 - (lat + adj_lat))
  ilon = int((lon + adj_lat) + 181.0)

  if (ilon==0) ilon = 360
  if (ilon==361) ilon = 1
  if (ilat==0) ilat = 1
  if (ilat==91) ilat = 90

  open(99,file=trim(inp%dirs%land_mask)//'/country.dat',status='old')
  if (ilat>1) then
    do i=1,ilat-1
      read(99,*)
    enddo
  endif
  read(99,*) ans
  country_id = ans(ilon)
  close(99)

endif

if (country_id==0) then

! Nearest diagonal

  xlat = lat+500-real(int(lat+500.0))
  if (xlat>0.5) then
    adj_lat = 1.0
  else
    adj_lat = -1.0
  endif

  xlon = lon+500-real(int(lon+500.0))
  if (xlon>0.5) then
    adj_lon = 1.0
  else
    adj_lon = -1.0
  endif

  ilat = int(91.0 - (lat + adj_lat))
  ilon = int((lon + adj_lat) + 181.0)

  if (ilon==0) ilon = 360
  if (ilon==361) ilon = 1
  if (ilat==0) ilat = 1
  if (ilat==91) ilat = 90

  open(99,file=trim(inp%dirs%land_mask)//'/country.dat',status='old')
    if (ilat>1) then
      do i=1,ilat-1
        read(99,*)
      enddo
    endif
    read(99,*) ans
    country_id = ans(ilon)
  close(99)

endif

!----------------------------------------------------------------------!
! Regional or country switch.                                          !
!----------------------------------------------------------------------!
if (.not.(l_regional)) country_id = country_id-mod(country_id,100)
!----------------------------------------------------------------------!

open(99,file=trim(inp%dirs%land_mask)//'/country_id.dat',status='old')
10    continue
  read(99,'(i6,5x,a15)') x,country_name
  if (x==country_id) goto 20
goto 10
20    continue
close(99)

end subroutine country







end module data
