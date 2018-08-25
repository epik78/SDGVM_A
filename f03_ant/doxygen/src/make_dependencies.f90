program make_dependencies

implicit none

integer, parameter :: nfiles = 3
character(len=100), dimension(nfiles) :: filenames
character(len=100), dimension(nfiles,100) :: uses
character(len=100) :: line
integer :: i,kode

filenames(1) = 'dims.f90'
filenames(2) = 'real_precision.f90'
filenames(3) = 'input_methods.f90'

i = 3
  open(10,file=trim(filenames(i)),status='old',iostat=kode)
if (kode==0) then
do
  read(10,'(a)',iostat=kode) line
  if (kode/=0) exit
  line = trim(line)
  if ((line(1:4)=='use ').or.(line(1:4)=='USE ')) then
    print*,trim(line)
  endif
enddo
else
  print*,'file does not exist'
  print*,trim(filenames(i))
endif




end program make_dependencies
