module file_class

  use dims

  implicit none
!  integer :: str_len
!  parameter(str_len = 1000)

  integer :: max_files
  parameter(max_files = 10000)

  type NoIds
    integer :: n
    integer, dimension(max_files) :: ids
    character(len=str_len), dimension(max_files) :: names
    logical,dimension(max_files) :: opened
	
    contains
      procedure :: open_old
      procedure :: openu_old
      procedure :: open
      procedure :: openu
      procedure :: close      ! by id
      procedure :: close_name
      procedure :: close_n
      procedure :: print
      procedure :: get_id     ! get id corresponing to a filename
      procedure :: get_id_n   ! get the n^th id
  end type NoIds

  contains

!**********************************************************************!
    subroutine filename_base(st1,st2)
!**********************************************************************!
    character(len=*) :: st1,st2
    integer :: ii

    ii = index(st1,'/',.true.)
    if (ii>0) then
      st2 = st1(ii+1:len(st1))
    else
      st2 = st1
    endif
    end subroutine filename_base





!**********************************************************************!
    function get_id(this,filename)
!**********************************************************************!
    class(NoIds) :: this
    integer :: get_id,kode,i
    character(len=*) :: filename
    character(len=200) :: st1,st2

    call filename_base(filename,st2)

    i = 0
    do
      i = i + 1
      if (i>this%n) then
        write(*,*) 'Error: Can''t find file with base:'
        write(*,*) trim(filename)
        stop
      endif
      if (this%opened(i)) then
        call filename_base(trim(this%names(i)),st1)
        if (trim(st1)==trim(st2)) then
          get_id = this%ids(i)
          exit
        endif
      endif
    enddo
    end function get_id





!**********************************************************************!
    function get_id_n(this,n)
!**********************************************************************!
    class(NoIds) :: this
    integer :: get_id_n,n

    get_id_n = this%ids(n)

    end function get_id_n





!**********************************************************************!
    subroutine open_old(this,filename,fid)
!**********************************************************************!
    class(NoIds) :: this
    integer :: fid,kode
    character(len=*) :: filename

    this%n = this%n + 1
    if (this%n>max_files) then
      write(*,*) 'Maximum number of files has been exceeded:'
      write(*,*) max_files
      stop
    endif

    open(newunit=fid, file=trim(filename), status='old', iostat=kode)
    if (kode/=0) then
      write(*,*) 'Error: file does not exist (open_old:file_io_methods):'
      write(*,*) trim(filename)
      stop
    endif
    this%ids(this%n) = fid
    this%names(this%n) = trim(filename)
    this%opened(this%n) = .true.
    end subroutine open_old





!**********************************************************************!
    subroutine openu_old(this,filename,fid)
!**********************************************************************!
    class(NoIds) :: this
    integer :: fid,kode
    character(len=*) :: filename

    this%n = this%n + 1
    if (this%n>max_files) then
      write(*,*) 'Maximum number of files has been exceeded:'
      write(*,*) max_files
      stop
    endif

    open(newunit=fid, file=trim(filename), status='old', &
 form='unformatted', iostat=kode)
    if (kode/=0) then
      write(*,*) 'Error: file does not exist (open_old:file_io_methods):'
      write(*,*) trim(filename)
      stop
    endif
    this%ids(this%n) = fid
    this%names(this%n) = trim(filename)
    this%opened(this%n) = .true.
    end subroutine openu_old





!**********************************************************************!
    subroutine open(this,filename,fid)
!**********************************************************************!
    class(NoIds) :: this
    integer :: kode,fid
    character(len=*) :: filename

    this%n = this%n + 1
    if (this%n>max_files) then
      write(*,*) 'Maximum number of files has been exceeded:'
      write(*,*) max_files
      stop
    endif

    open(newunit=fid, file=trim(filename), status='unknown', iostat=kode)
    if (kode/=0) then
      write(*,*) 'Error: file does not exist (open:file_io_methods):'
      write(*,*) trim(filename)
      stop
    endif
    this%ids(this%n) = fid
    this%names(this%n) = trim(filename)
    this%opened(this%n) = .true.
    end subroutine open




	
!**********************************************************************!
    subroutine openu(this,filename,fid)
!**********************************************************************!
    class(NoIds) :: this
    integer :: kode,fid
    character(len=*) :: filename

    this%n = this%n + 1
    if (this%n>max_files) then
      write(*,*) 'Maximum number of files has been exceeded:'
      write(*,*) max_files
      stop
    endif

    open(newunit=fid, file=trim(filename), status='unknown', &
 form='unformatted', iostat=kode)
    if (kode/=0) then
      write(*,*) 'Error: file does not exist (open:file_io_methods):'
      write(*,*) trim(filename)
      stop
    endif
    this%ids(this%n) = fid
    this%names(this%n) = trim(filename)
    this%opened(this%n) = .true.
    end subroutine openu




	
!**********************************************************************!
    subroutine close_name(this,filename)
!**********************************************************************!
    class(NoIds) :: this
    integer :: i,kode
    character(len=*) :: filename

    i = 0
    do
      i = i + 1
      if (i>this%n) then
        write(*,*) 'Error: Can''t find file to close:'
        write(*,*) trim(filename)
        stop
      endif
      if (this%opened(i)) then
        if (trim(this%names(i))==trim(filename)) then
          close(this%ids(i), iostat=kode)
          if (kode/=0) then
            write(*,*) 'Error: Can''t close file:'
            write(*,*) trim(filename)
            stop
          endif
          this%opened(i) = .false.
          exit
        endif
      endif
    enddo
    end subroutine close_name





!**********************************************************************!
    subroutine close_n(this,n_id)
!**********************************************************************!
    class(NoIds) :: this
    integer :: n_id,kode
    character(len=str_len) :: filename

    close(this%ids(n_id),iostat=kode)
    if (kode/=0) then
      write(*,*) 'Error: Can''t close file:'
      write(*,*) trim(filename)
      stop
    endif
    this%opened(n_id) = .false.
    end subroutine close_n





!**********************************************************************!
    subroutine close(this,id)
!**********************************************************************!
    class(NoIds) :: this
    integer :: id,kode,i

    i = 0
    do
      i = i + 1
      if (i>this%n) then
        write(*,*) 'file_io_methods : close : line 180'
        write(*,*) 'Error: Can''t find file to close:'
        stop
      endif
      
      if (this%ids(i)==id) then
        close(id,iostat=kode)
        if (kode/=0) then
          write(*,*) 'file_io_methods : close : line 188'
          write(*,*) 'Error: Can''t close file:'
          write(*,*) this%ids(i),trim(this%names(i))
          stop
        endif
        this%opened(i) = .false.
        exit
      endif
    enddo
    end subroutine close





!**********************************************************************!
    subroutine print(this)
!**********************************************************************!
    class(NoIds) :: this
    integer :: i
    do i=1,this%n
      if (this%opened(i)) then
        write(*,'(i5,''  open   '',i6,''  '',a)') i,this%ids(i),trim(this%names(i))
      else
        write(*,'(i5,''  closed '',i6,''  '',a)') i,this%ids(i),trim(this%names(i))
      endif
    enddo
    end subroutine print

end module file_class



module file_object
use file_class
!
!type Units
!
!  type(NoIds) :: default_outputs
!  type(NoIds) :: daily_outputs
!  type(NoIds) :: monthly_outputs
!  type(NoIds) :: yearly_outputs
!  type(NoIds) :: snapshots
!  type(NoIds) :: misc
!  integer :: next
!
!end type Units
!
type(Noids) :: fun, funMisc
end module file_object
