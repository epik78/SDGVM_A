module fnames_class

  use dims
  implicit none

  type fmt

    character(len=str_len) :: name = 'None'
    logical                :: averaged = .false.
    logical                :: output_tile = .false.
    logical                :: output_subtile = .false.
    character(len=str_len) :: daily   = '(31f7.2)'
    character(len=str_len) :: monthly = '(i4,1x,12f8.2,f10.1)'
    character(len=str_len) :: yearly  = '(31f7.2)'
  end type fmt


  type Fnames
    integer :: n
    type(fmt), dimension(max_outputs) :: fmt

    contains 
      procedure :: set_names
      procedure :: set_outputs
  end type Fnames

  contains


    subroutine set_outputs(this,st1)
    class(Fnames) :: this
    character(len=str_len) :: st1
    integer :: i

    do i=1,max_outputs
!      print*,trim(otags(i))
    enddo
    print*,'in set outputs'


    end subroutine set_outputs


    subroutine set_names(this)
    class(Fnames) :: this

    this%n = 25

    this%fmt(1)%name = 'lai'
    this%fmt(1)%averaged = .true.
    this%fmt(2)%name = 'rof'
    this%fmt(3)%name = 'evt'
    this%fmt(4)%name = 'trn'
    this%fmt(5)%name = 'npp'
    this%fmt(6)%name = 'gpp'
    this%fmt(7)%name = 'srp'
    this%fmt(8)%name = 'nep'
    this%fmt(9)%name = 'tmp'
    this%fmt(9)%averaged = .true.
    this%fmt(10)%name = 'prc'
    this%fmt(11)%name = 'hum'
    this%fmt(11)%averaged = .true.
    this%fmt(12)%name = 'nps'
    this%fmt(12)%averaged = .true.
    this%fmt(13)%name = 'swf'
    this%fmt(13)%averaged = .true.
    this%fmt(14)%name = 'pet'
    this%fmt(15)%name = 'int'
    this%fmt(16)%name = 'bse'
    this%fmt(17)%name = 'ssm'
    this%fmt(17)%averaged = .true.
    this%fmt(18)%name = 'swc'
    this%fmt(18)%averaged = .true.
    this%fmt(19)%name = 'rsp'
    this%fmt(20)%name = 'qdr'
    this%fmt(20)%averaged = .true.
    this%fmt(21)%name = 'qdf'
    this%fmt(21)%averaged = .true.
    this%fmt(22)%name = 'lfn'
    this%fmt(22)%averaged = .true.
    this%fmt(23)%name = 'lfl'
    this%fmt(23)%averaged = .true.
    this%fmt(24)%name = 'cld'
    this%fmt(24)%averaged = .true.
    this%fmt(25)%name = 'fpr'
    this%fmt(25)%averaged = .true.


    end subroutine set_names

end module fnames_class

module fnames_object
  use fnames_class

  type(Fnames) :: fnms

end module fnames_object
