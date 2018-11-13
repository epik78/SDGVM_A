module input_file
! To add a new section/tag add a new tag_list() string containing the name of 
! the new section in the read_input_file subroutine, and increase nTags by 1.
! Add a new type definition here using the section name as the type name
! and populate it with values that you require.
! Finally add the read method to the subroutine set_input_parameters 
! using the tags that you will use within the section to select the read method.

use dims
use real_precision
implicit none

private :: read_input_file

type, private :: DefIntYears
  integer :: n = 0
  integer, dimension(max_years) :: val = 0
end type DefIntYears

type, private :: DefChrPfts
  integer :: n = 0
  character(len=str_len), dimension(max_outputs) :: tag = ''
  integer, dimension(max_outputs) :: map = 0
end type DefChrPfts

type, private :: Land_use
  logical :: read_from_landuse_dir = .true.
  integer :: n = 0
  integer, dimension(max_years) :: year
  real(dp), dimension(max_years) :: map
end type Land_use

type, private :: Dirs
  logical :: input_argument  = .false.
  logical :: output_argument = .false.
  character(len=str_len) :: input
  character(len=str_len) :: output
  character(len=str_len) :: climate
  character(len=str_len) :: stats
  character(len=str_len) :: soil
  character(len=str_len) :: land_use
  character(len=str_len) :: co2
  character(len=str_len) :: land_mask
  character(len=str_len) :: fert
  character(len=str_len) :: irri
  character(len=str_len) :: params 
end type Dirs

type, private :: Run
  real(dp) :: co2_constant = 29.62
  logical :: s070607 = .true.
  logical :: s140129 = .false.
  logical :: subdaily = .true.
  logical :: read_daily_co2 = .false. 
  logical :: read_par = .true.
  logical :: par_scaling = .true.
  logical :: read_clump = .false.
  logical :: calc_zen = .false.
  logical :: no_soil_water_limitation = .false.
  integer :: cstype = 0
  integer :: ncalc_type = 0
  integer :: ttype = 0
  integer :: vcmax_type = 0
  integer :: soilp_map = 0
  integer :: gs_func = 0
  
  logical :: read_in_state = .false.
  logical :: output_state = .false.
  integer :: random_seed = 1
    
  integer :: spinup_length = 0
  integer :: spinup_year0 = 0
  integer :: spinup_cycle_length = 0
  integer :: year0 = 1901
  integer :: yearf = 2000
end type Run

type, private :: Pft_parameters
  character(len=str_len) :: tag = ''
  logical    :: c3c4    = .true.
  integer    :: phen    = 0
  real(dp)   :: mix     = 0.0
  real(dp)   :: crop    = 0.0
  integer    :: d2h     = 0   ! days to harvest
  integer    :: mort    = 1
  real(dp)   :: wden    = 0.0
  real(dp)   :: xyl     = 0.0
  real(dp)   :: pdif    = 0.0
  real(dp)   :: sla     = 0.0 ! specific leaf area ()
  integer    :: lls     = 0   ! Leaf life span (days)
  integer    :: sls     = 0   ! Stem life span (days)
  integer    :: rls     = 0   ! Root life span (days)
  real(dp)   :: lmor    = 0.0 ! Leaf mortality
  real(dp)   :: lrat    = 0.0 ! Leaf rate
  integer    :: bbmem   = 0   ! Budburst memory
  real(dp)   :: bb0     = 0.0 ! Budburst lower limit
  real(dp)   :: bbmax   = 0.0 ! Budburst upper limit
  real(dp)   :: bblim   = 0.0 ! Budburst limit
  integer    :: senm    = 0   ! Senescence memory
  integer    :: sens    = 0   ! Senescence sensitivity
  real(dp)   :: senlim  = 0.0 ! Senescence limit
  integer    :: chill   = 0   ! Chilling
  integer    :: dschill = 0
  real(dp)   :: stemx   = 0.0
  real(dp)   :: gr0     = 0.0
  real(dp)   :: grf     = 0.0
  real(dp)   :: ppm0    = 0   ! Plant density

  real(dp)   :: can_clump = 0.0
  real(dp)   :: vna       = 0.0
  real(dp)   :: vnb       = 0.0
  real(dp)   :: jva       = 0.0
  real(dp)   :: jvb       = 0.0
  real(dp)   :: g0        = 0.0
  real(dp)   :: g1        = 0.0

  real(dp),dimension(2)   :: sowthresh = 0.0
  real(dp),dimension(2)   :: lethal = 0.0
  real(dp),dimension(9)   :: cardinal = 0.0
  real(dp),dimension(2)   :: croptype = 0.0
  real(dp),dimension(6)   :: photoperiod = 0.0
  real(dp),dimension(4)   :: croprange = 0.0
  real(dp),dimension(6)   :: cropphen = 0.0
  real(dp),dimension(3)   :: irrig = 0.0
  real(dp),dimension(3)   :: sowday = 0.0
  real(dp),dimension(2,3) :: cropgdd = 0.0
  real(dp),dimension(6)   :: fert = 0.0
  real(dp)                :: optlai = 0.0
  real(dp)                :: harvindx = 0.0
  real(dp)                :: limdharv = 0.0
end type Pft_parameters

type, private :: Soil
  real(dp) :: sand  = -1.0
  real(dp) :: silt  = -1.0
  real(dp) :: bulk  = -1.0
  real(dp) :: orgc  = -1.0
  real(dp) :: wilt  = -1.0
  real(dp) :: field = -1.0
  real(dp) :: sat   = -1.0
  real(dp) :: topsl = -1.0
  real(dp) :: depth = -1.0
end type Soil

type, private :: Pft_map
  integer :: i = 0
  integer :: n = 0
  real(dp), dimension(max_pftps) :: percent = 0.0
  character(len=str_len), dimension(max_pftps) :: pft = ''
end type Pft_map

type, private :: Output
  integer :: nyears = 10

  integer :: tile_yearly  = 0
  integer :: tile_monthly = 0
  integer :: tile_daily   = 0

  integer :: pft_yearly  = 0
  integer :: pft_monthly = 0
  integer :: pft_daily   = 0


  type(DefChrPfts) :: tile_vars_yearly
  type(DefChrPfts) :: tile_vars_monthly
  type(DefChrPfts) :: tile_vars_daily

  type(DefChrPfts) :: pft_vars_yearly
  type(DefChrPfts) :: pft_vars_monthly
  type(DefChrPfts) :: pft_vars_daily

  type(DefIntYears) :: state_years
end type Output


type, private :: Sites
  logical :: create_land_sites_only = .false.
  logical :: n_range_argument = .false.
  logical :: filename_argument = .false.
  character(len=str_len) :: style
  integer :: site0 = 1
  integer :: sitef = 100000
  integer :: nlist = 0 ! not set in input file
  real(dp), dimension(max_sites,2) :: lat_lon
  real(dp), dimension(max_sites,2) :: list
  real(dp), dimension(4) :: box
  real(dp), dimension(2) :: res
  character(len=str_len) :: filename
end type Sites


type, private :: Input
  type(Dirs) :: dirs
  type(Pft_map), dimension(max_pftps) :: pft_mapping
  type(Run) :: run
  type(Pft_parameters), dimension(max_pftps) :: pft
  type(Output) :: output
  type(Soil) :: soil
  type(Sites) :: sites
  type(Land_use) :: land_use
  integer :: i
  integer :: npft = 0
  integer :: npft_mapping = 0
contains 
  procedure :: read_input_file
end type Input

type(Input) :: inp

contains
  
!**********************************************************************!
  subroutine read_input_file(this,input_filename)
!**********************************************************************!
    class(input) :: this
    character(len=*) :: input_filename
    character(len=str_len) :: open_tag,close_tag,param_name
    character(len=str_len) :: param_values
    integer :: out1 = 1, i1,kode = 99
    character(len=str_len),dimension(max_pftps) :: tag_list
    integer :: param_or_closing_tag,iostatus
    integer, parameter :: nTags = 8
!----------------------------------------------------------------------!
    tag_list(1) = 'dirs'
    tag_list(2) = 'pft'
    tag_list(3) = 'pft_mapping'
    tag_list(4) = 'run'
    tag_list(5) = 'output'
    tag_list(6) = 'soil'
    tag_list(7) = 'sites'
    tag_list(8) = 'land_use'

    open(11,file=input_filename,status='old',iostat=kode)
    if (kode/=0) then
      write(*,'(''Input file does not exist'')')
      write(*,*) trim(input_filename)
      stop
    endif

! Loop through the tags
    do    
      call find_next_open_tag(open_tag,iostatus)
      if (iostatus /= 0) exit
!      print*,'---------- opening tag ',trim(open_tag)
      if (.not.check_tag_exists(open_tag,nTags,tag_list)) then
        write(*,*) 'Error opening tag doesn''t exist!!!'
        write(*,*) '<',trim(open_tag),'>'
        stop
      endif
! Loop through parameters within the tag
      do
        call read_tag_parameter(param_name,param_values,open_tag,param_or_closing_tag)
        if (param_or_closing_tag == 0) exit
        call set_input_parameters(this,open_tag,param_name,param_values)
      enddo
    enddo
    
    close(11)
    
  end subroutine read_input_file
  



!**********************************************************************!
  subroutine set_input_parameters(this,tag,param_name,param_values)
!**********************************************************************!
  class(input) :: this
  character(len=str_len) :: tag,param_name
  character(len=str_len) :: param_values
  integer :: i,j,ii,jj,n
  real(dp), dimension(max_pftps) :: percent
  character(len=str_len), dimension(max_pftps) :: pft
  integer, dimension(max_pftps), save :: pftind
  character(len=str_len), dimension(max_pftps), save :: pfttag
  integer, save :: pftn = 0, nlist = 0
  logical :: match
!----------------------------------------------------------------------!
  
  select case (trim(tag))
    case ('land_use')
      select case (trim(param_name))
        case('read_from_landuse_dir')
          read(param_values,*) this%land_use%read_from_landuse_dir
        case('list')
          if (nfields(param_values)==2) then
            this%land_use%n = this%land_use%n + 1
            read(param_values,*) this%land_use%year(this%land_use%n), &
 this%land_use%map(this%land_use%n)
          else
            write(*,*) 'Error in case sites!!!'
            write(*,*) 'Need 2 inputs for list.'
            write(*,*) trim(param_values)
          endif

        case default    
          write(*,*) 'Error in case land_use!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
          stop
      end select
    case ('sites')
      select case (trim(param_name))
        case('create_land_sites_only')
          read(param_values,*) this%sites%create_land_sites_only
        case('style')
          read(param_values,*) this%sites%style
        case('n_range')
          if (nfields(param_values)==2) then
            read(param_values,*) this%sites%site0,this%sites%sitef
            if (this%sites%site0<1) then
              write(*,*) 'Error in case sites!!!'
              write(*,*) 'First site number of range must be >0.'
              write(*,*) this%sites%site0
              stop
            endif
            if (this%sites%sitef<this%sites%site0) then
              write(*,*) 'Error in case sites!!!'
              write(*,*) 'Second site number of range must be greater or equal to the first.'
              write(*,*) this%sites%site0,this%sites%sitef
              stop
            endif
          elseif (trim(param_values)=='ARGUMENT') then
            this%sites%n_range_argument = .true.
          elseif (nfields(param_values)/=0) then
            write(*,*) 'Error in case sites!!!'
            write(*,*) 'Need either 2 integers or ''ARGUMENT'' for n_range.'
            write(*,*) trim(param_values)
            stop
          endif
        case('list')
          nlist = nlist + 1
          this%sites%nlist = nlist
          read(param_values,*) this%sites%list(nlist,1),this%sites%list(nlist,2)
        case('box')
          if (nfields(param_values)==4) then
            read(param_values,*) this%sites%box
          else
            write(*,*) 'Error in case sites!!!'
            write(*,*) 'Need 4 inputs for box.'
            write(*,*) trim(param_values)
            stop
          endif
        case('res')
          read(param_values,*) this%sites%res
        case('filename')
          if (trim(param_values)=='ARGUMENT') then
            this%sites%filename_argument = .true.
          else
            read(param_values,'(a)') this%sites%filename
          endif
        case default    
          write(*,*) 'Error in case sites!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
          stop
      end select
    case ('soil')
      select case (trim(param_name))
        case('sand')
          read(param_values,*) this%soil%sand
        case('silt')
          read(param_values,*) this%soil%silt
        case('bulk')
          read(param_values,*) this%soil%bulk
        case('orgc')
          read(param_values,*) this%soil%orgc
        case('wilt')
          read(param_values,*) this%soil%wilt
        case('field')
          read(param_values,*) this%soil%field
        case('sat')
          read(param_values,*) this%soil%sat
        case('topsl')
          read(param_values,*) this%soil%topsl
        case('depth')
          read(param_values,*) this%soil%depth

        case default    
          write(*,*) 'Error in case soil!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
          stop
      end select

    case ('pft')
      select case (trim(param_name))
        case('NAMES')
          pftn = nfields(param_values)
          if (pftn>0) then
            read(param_values,*) (pfttag(i),i=1,pftn)
            do i=1,pftn
              match = .false.
              do j=1,this%npft
                if (trim(this%pft(j)%tag)==trim(pfttag(i))) then
                  match = .true.
                  jj = j
                  exit
                endif
              enddo
              if (match) then
                pftind(i) = j
                pfttag(i) = trim(pfttag(jj))
              else
                this%npft = this%npft + 1
                this%pft(this%npft)%tag = trim(pfttag(i))
                pftind(i) = this%npft
                pfttag(i) = trim(pfttag(i))
              endif
            enddo
          else
            write(*,*) 'No pft tags given.'
            stop
          endif
        case('mix')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%mix,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''mix'''
            stop
          endif
        case('c3c4')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%c3c4,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''c3c4'''
            stop
          endif
        case('phen')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%phen,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''phen'''
            stop
          endif
        case('crop')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%crop,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''crop'''
            stop
          endif
        case('d2h') ! days to harvest
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%d2h,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''d2h'''
            stop
          endif
        case('mort')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%mort,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''mort'''
            stop
          endif
        case('wden')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%wden,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''wden'''
            stop
          endif
        case('xyl')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%xyl,j=1,pftn)
          else
            write(*,*) 'Not enough fields in ''xyl'''
            stop
          endif
        case('pdif')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%pdif,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''pdif'''
            stop
          endif
        case('sla') ! specific leaf area ()
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sla,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sla'''
            stop
          endif
        case('lls') ! Leaf life span (days)
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%lls,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''lls'''
            stop
          endif
        case('sls') ! Stem life span (days)
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sls,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sls'''
            stop
          endif
        case('rls') ! Root life span (days)
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%rls,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''rls'''
            stop
          endif
        case('lmor') ! Leaf mortality
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%lmor,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''lmor'''
            stop
          endif
        case('lrat') ! Leaf rate
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%lrat,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''lrat'''
            stop
          endif
        case('bbmem') ! Budburst memory
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%bbmem,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''bbmem'''
            stop
          endif
        case('bb0')  ! Budburst lower limit
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%bb0,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''bb0'''
            stop
          endif
        case('bbmax') ! Budburst upper limit
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%bbmax,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''bbmax'''
            stop
          endif
        case('bblim') ! Budburst limit
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%bblim,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''bblim'''
            stop
          endif
        case('senm') ! Senescence memory
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%senm,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''senm'''
            stop
          endif
        case('sens')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sens,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sens'''
            stop
          endif
        case('senlim')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%senlim,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''senlim'''
            stop
          endif
        case('chill')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%chill,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''chill'''
            stop
          endif
        case('dschill')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%dschill,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''dschill'''
            stop
          endif
        case('stemx')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%stemx,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''stemx'''
            stop
          endif
        case('gr0')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%gr0,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''gr0'''
            stop
          endif
        case('grf')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%grf,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''grf'''
            stop
          endif
        case('ppm0')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%ppm0,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''ppm0'''
            stop
          endif


        case('can_clump')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%can_clump,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''can_clump'''
            stop
          endif
        case('vna')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%vna,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''vna'''
            stop
          endif
        case('vnb')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%vnb,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''vnb'''
            stop
          endif
        case('jva')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%jva,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''jva'''
            stop
          endif
        case('jvb')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%jvb,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''jvb'''
            stop
          endif
        case('g0')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%g0,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''g0'''
            stop
          endif
        case('g1')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%g1,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''g1'''
            stop
          endif

        case('sowthresh(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sowthresh(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sowthresh(1)'''
            stop
          endif
        case('sowthresh(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sowthresh(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sowthresh(2)'''
            stop
          endif
        case('lethal(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%lethal(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''lethal(1)'''
            stop
          endif
        case('lethal(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%lethal(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''lethal(2)'''
            stop
          endif
        case('cardinal(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(1)'''
            stop
          endif
        case('cardinal(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(2)'''
            stop
          endif
        case('cardinal(3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(3)'''
            stop
          endif
        case('cardinal(4)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(4),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(4)'''
            stop
          endif
        case('cardinal(5)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(5),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(5)'''
            stop
          endif
        case('cardinal(6)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(6),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(6)'''
            stop
          endif
        case('cardinal(7)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(7),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(7)'''
            stop
          endif
        case('cardinal(8)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(8),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(8)'''
            stop
          endif
        case('cardinal(9)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cardinal(9),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cardinal(9)'''
            stop
          endif
        case('croptype(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%croptype(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''croptype(1)'''
            stop
          endif
        case('croptype(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%croptype(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''croptype(2)'''
            stop
          endif
        case('photoperiod(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%photoperiod(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''photoperiod(1)'''
            stop
          endif
        case('photoperiod(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%photoperiod(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''photoperiod(2)'''
            stop
          endif
        case('photoperiod(3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%photoperiod(3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''photoperiod(3)'''
            stop
          endif
        case('photoperiod(4)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%photoperiod(4),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''photoperiod(4)'''
            stop
          endif
        case('photoperiod(5)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%photoperiod(5),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''photoperiod(5)'''
            stop
          endif
        case('photoperiod(6)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%photoperiod(6),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''photoperiod(6)'''
            stop
          endif
        case('croprange(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%croprange(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''croprange(1)'''
            stop
          endif
        case('croprange(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%croprange(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''croprange(2)'''
            stop
          endif
        case('croprange(3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%croprange(3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''croprange(3)'''
            stop
          endif
        case('croprange(4)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%croprange(4),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''croprange(4)'''
            stop
          endif
        case('cropphen(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropphen(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropphen(1)'''
            stop
          endif
        case('cropphen(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropphen(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropphen(2)'''
            stop
          endif
        case('cropphen(3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropphen(3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropphen(3)'''
            stop
          endif
        case('cropphen(4)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropphen(4),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropphen(4)'''
            stop
          endif
        case('cropphen(5)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropphen(5),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropphen(5)'''
            stop
          endif
        case('cropphen(6)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropphen(6),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropphen(6)'''
            stop
          endif
        case('irrig(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%irrig(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''irrig(1)'''
            stop
          endif
        case('irrig(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%irrig(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''irrig(2)'''
            stop
          endif
        case('irrig(3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%irrig(3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''irrig(3)'''
            stop
          endif
        case('sowday(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sowday(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sowday(1)'''
            stop
          endif
        case('sowday(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sowday(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sowday(2)'''
            stop
          endif
        case('sowday(3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%sowday(3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''sowday(3)'''
            stop
          endif
        case('cropgdd(1,1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropgdd(1,1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropgdd(1,1)'''
            stop
          endif
        case('cropgdd(1,2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropgdd(1,2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropgdd(1,2)'''
            stop
          endif
        case('cropgdd(1,3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropgdd(1,3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropgdd(1,3)'''
            stop
          endif
        case('cropgdd(2,1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropgdd(2,1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropgdd(2,1)'''
            stop
          endif
        case('cropgdd(2,2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropgdd(2,2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropgdd(2,2)'''
            stop
          endif
        case('cropgdd(2,3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%cropgdd(2,3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''cropgdd(2,3)'''
            stop
          endif
        case('fert(1)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%fert(1),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''fert(1)'''
            stop
          endif
        case('fert(2)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%fert(2),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''fert(2)'''
            stop
          endif
        case('fert(3)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%fert(3),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''fert(3)'''
            stop
          endif
        case('fert(4)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%fert(4),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''fert(4)'''
            stop
          endif
        case('fert(5)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%fert(5),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''fert(5)'''
            stop
          endif
        case('fert(6)')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%fert(6),j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''fert(6)'''
            stop
          endif
        case('optlai')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%optlai,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''optlai'''
            stop
          endif
        case('harvindx')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%harvindx,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''harvindx'''
            stop
          endif
        case('limdharv')
          i = nfields(param_values)
          if (i==pftn) then
            read(param_values,*) (this%pft(pftind(j))%limdharv,j=1,pftn)
          else
            write(*,*) 'Wrong number of fields in ''limdharv'''
            stop
          endif

        case default    
          write(*,*) 'Error in case pft!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
          stop
      end select
    
    case ('output')
      select case (trim(param_name))
        case('nyears')
          read(param_values,*) this%output%nyears

        case('tile_yearly')
          read(param_values,*) this%output%tile_yearly
        case('tile_monthly')
          read(param_values,*) this%output%tile_monthly
        case('tile_daily')
          read(param_values,*) this%output%tile_daily

        case('tile_vars_yearly')
          i = nfields(param_values)
          if (i>0) then
            this%output%tile_vars_yearly%n = i
            read(param_values,*) (this%output%tile_vars_yearly%tag(j),j=1,i)
          endif
        case('tile_vars_monthly')
          i = nfields(param_values)
          if (i>0) then
            this%output%tile_vars_monthly%n = i
            read(param_values,*) (this%output%tile_vars_monthly%tag(j),j=1,i)
          endif
        case('tile_vars_daily')
          i = nfields(param_values)
          if (i>0) then
            this%output%tile_vars_daily%n = i
            read(param_values,*) (this%output%tile_vars_daily%tag(j),j=1,i)
          endif

        case('pft_yearly')
          read(param_values,*) this%output%pft_yearly
        case('pft_monthly')
          read(param_values,*) this%output%pft_monthly
        case('pft_daily')
          read(param_values,*) this%output%pft_daily

        case('pft_vars_yearly')
          i = nfields(param_values)
          if (i>0) then
            this%output%pft_vars_yearly%n = i
            read(param_values,*) (this%output%pft_vars_yearly%tag(j),j=1,i)
          endif
        case('pft_vars_monthly')
          i = nfields(param_values)
          if (i>0) then
            this%output%pft_vars_monthly%n = i
            read(param_values,*) (this%output%pft_vars_monthly%tag(j),j=1,i)
          endif
        case('pft_vars_daily')
          i = nfields(param_values)
          if (i>0) then
            this%output%pft_vars_daily%n = i
            read(param_values,*) (this%output%pft_vars_daily%tag(j),j=1,i)
          endif

        case('state_years')
          i = nfields(param_values)
          if (i>0) then
            this%output%state_years%n = i
            read(param_values,*) (this%output%state_years%val(j),j=1,i)
          endif
        case default    
          write(*,*) 'Error in case output!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
          stop
      end select

    case ('run')
      select case (trim(param_name))
        case ('co2_constant')
          read(param_values,*) this%run%co2_constant
        case('s070607')
          read(param_values,*) this%run%s070607
        case('s140129')
          read(param_values,*) this%run%s140129
        case('subdaily')
          read(param_values,*) this%run%subdaily
        case('read_daily_co2')
          read(param_values,*) this%run%read_daily_co2
        case('read_par')
          read(param_values,*) this%run%read_par
        case('par_scaling')
          read(param_values,*) this%run%par_scaling
        case('read_clump')
          read(param_values,*) this%run%read_clump
        case('calc_zen')
          read(param_values,*) this%run%calc_zen
          
        case('no_soil_water_limitation')
          read(param_values,*) this%run%no_soil_water_limitation
        case('cstype')
          read(param_values,*) i
          if (i>=10) then
            this%run%cstype = int(real(i)/10.0)
            i = i - 10*this%run%cstype
          endif
          if (i>-1)  this%run%ncalc_type = i 
        case('ttype')
          read(param_values,*) i
          if (i>=10) then
            this%run%ttype = int(real(i)/10.0)
            i = i - 10*this%run%ttype
          endif
          if (i>-1)  this%run%vcmax_type = i
        case('soilp_map')
          read(param_values,*) this%run%soilp_map
        case('gs_func')
          read(param_values,*) this%run%gs_func
          
        case('read_in_state')
          read(param_values,*) this%run%read_in_state
        case('output_state')
          read(param_values,*) this%run%output_state
        case('random_seed')
          read(param_values,*) this%run%random_seed

        case('spinup_length')
          read(param_values,*) this%run%spinup_length
        case('spinup_year0')
          read(param_values,*) this%run%spinup_year0
        case('spinup_cycle_length')
          read(param_values,*) this%run%spinup_cycle_length
        case('year0')
          read(param_values,*) this%run%year0
        case('yearf')
          read(param_values,*) this%run%yearf

        case default    
          write(*,*) 'Error in case run!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
          stop
      end select
      
    case ('dirs')
      select case (trim(param_name))
        case ('input')
          if (trim(param_values)=='ARGUMENT') then
            this%dirs%input_argument = .true.
          else
            read(param_values,'(a)') this%dirs%input
          endif
        case ('output')
          if (trim(param_values)=='ARGUMENT') then
            this%dirs%output_argument = .true.
          else
            read(param_values,'(a)') this%dirs%output
          endif
        case ('climate')
          read(param_values,'(a)') this%dirs%climate
        case ('stats')
          read(param_values,'(a)') this%dirs%stats
        case ('soil')
          read(param_values,'(a)') this%dirs%soil
        case ('land_use')
          read(param_values,'(a)') this%dirs%land_use
        case ('co2')
          read(param_values,'(a)') this%dirs%co2
        case ('land_mask')
          read(param_values,'(a)') this%dirs%land_mask
        case ('fert')
          read(param_values,'(a)') this%dirs%fert
        case ('irri')
          read(param_values,'(a)') this%dirs%irri
        case ('params')
          read(param_values,'(a)') this%dirs%params

        case default
          write(*,*) 'Error in case dirs!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
          stop
      end select

    case ('pft_mapping')
      select case (trim(param_name))
        case ('map')
          if (mod(nfields(param_values),2)==0) then
            write(*,*) &
            'Unpaired percentage cover pairs in landcover-pft mapping:'
            write(*,*) trim(tag),' ',trim(param_values)
            stop
          endif
          read(param_values,*) i,(percent(j),pft(j),j=1,(nfields(param_values)-1)/2)
          this%npft_mapping = this%npft_mapping + 1
          do j=1,(nfields(param_values)-1)/2
            this%pft_mapping(this%npft_mapping)%percent(j) = percent(j)
            this%pft_mapping(this%npft_mapping)%pft(j) = pft(j)
          enddo
          this%pft_mapping(this%npft_mapping)%n = (nfields(param_values)-1)/2
          this%pft_mapping(this%npft_mapping)%i = i
        case default
          write(*,*) 'Error in case pft-mapping!!!'
          write(*,*) 'Need to define read method for parameter ''',trim(param_name),'''!!!'
          write(*,*) 'This is done in the ''input_file'' module in the subroutine ''set_input_parameters''!!!'
      end select
    case default
      write(*,*) 'Stuck in a case statement in ''input_file'''
      write(*,*) trim(tag),' ',trim(param_name),' ',trim(param_values)
      stop
  end select

  end subroutine set_input_parameters
  
  
  
  
  
!**********************************************************************!
  function nfields(text)
!**********************************************************************!
  character(len=*) :: text
  integer :: nfields,nwords,pos,i 
!----------------------------------------------------------------------!

  pos = 1 
  nfields = 0
  
  do 
    i = verify(text(pos:),' ')   !-- Find next non-blank.
    if (i == 0) exit             !-- No word found. 
    nfields = nfields + 1        !-- Found something. 
    pos = pos + i - 1            !-- Move to start of the word. 
    i = scan(text(pos:), ' ')    !-- Find next blank. 
    if (i == 0) exit             !-- No blank found. 
    pos = pos + i - 1            !-- Move to the blank. 
  enddo

  end function nfields
  
  
  
  
  
!**********************************************************************!
  subroutine read_tag_parameter(param_name,param_values,tag,param_or_closing_tag)
!**********************************************************************!
  character(len=str_len) :: line
  character(len=str_len) :: param_name,tag
  character(len=str_len) :: param_values
  integer :: param_or_closing_tag,iostatus
  integer :: ii
!----------------------------------------------------------------------!

  do
    read(11,'(a)',iostat=iostatus) line
    line = adjustl(line)
    if (iostatus<0) then
      write(*,*) 'End of file reached in ''read_tag_parameter''!!!'
      stop
    endif
  
    if (len_trim(line) == str_len) then
      write(*,*) 'Line in input file is too long!!!',str_len
      write(*,*) line
      stop
    endif
  
    ii = index(line,'!')
    if (ii > 0)  line = line(1:ii-1)
  
! Leave the loop if line appears to be a parameter value or a close tab
! ie. first character is not a '!' or a ' '.
    if (.not.((line(1:1) == '!').or.(line(1:1) == ' '))) exit
  enddo
  
  ii = index(line,'::')
  if (ii > 0) then
! Set parameter value  
    param_name = trim(line(1:ii-1))
    param_or_closing_tag = 1
    param_values = adjustl(trim(line(ii+2:len(line))))
  elseif (line(1:2) == '</') then
! Closing tag
    ii = index(line,'>')
    param_name = line(3:ii-1)
! Check closing tag matches opening tag
    if (param_name /= tag) then
      write(*,*) 'Closing tag doesn''t match opening tag.'
      write(*,*) trim(tag),' is not equal to ',trim(param_name)
      stop
    endif
    param_or_closing_tag = 0
  else
    write(*,*) 'Error in reading input file tag !',trim(tag),'!'
    write(*,*) 'Lines must be either 1) blanks 2) first character '
    write(*,*) 'must be ''!'' 3) line should contain parameter '
    write(*,*) 'name, ''::'', and parameter value.'
    write(*,*) trim(line)
    stop
  endif
    
  end subroutine read_tag_parameter  
  
  
  
  
  
!**********************************************************************!
  function check_tag_exists(tag,nTags,tag_list) result(check)
!**********************************************************************!
    character(len=str_len) :: tag
    character(len=str_len),dimension(max_pftps) :: tag_list
    integer :: nTags,leq,i
    logical :: check
!----------------------------------------------------------------------!

    check = .false.
    do i=1,nTags
      if (tag == tag_list(i)) check = .true.
    enddo
  end function check_tag_exists
  
  
  
!**********************************************************************!
  subroutine find_next_open_tag(tag,iostatus)
!**********************************************************************!
    character(len=str_len) :: line
    character(len=str_len) :: tag
    integer :: iostatus,ii
!----------------------------------------------------------------------!
   
    do
      read(11,'(a)',iostat=iostatus) line
      if (iostatus /= 0) exit
      if (iostatus<0) then
        write(*,*) 'End of file reached in ''find_next_open_tag''!!!'
        stop
      endif

      if (len_trim(line) == str_len) then
        write(*,*) 'Line in input file is too long!!!'
        write(*,*) line
        stop
      endif
      line = adjustl(line)
      ii = index(line,'!')
      if (ii > 0)  line = line(1:ii-1)

      if (line(1:1) == '<') exit
    enddo
    
    if (iostatus == 0) then
      if (index(line,'>') == 0) then
        write(*,*) 'No end of tag character ''>''!!!'
        write(*,*) line
        stop
      endif
      tag = line(2:index(line,'>')-1)
    endif
        
  end subroutine find_next_open_tag

end module input_file
