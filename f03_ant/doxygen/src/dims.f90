module dims

integer, parameter :: max_cohorts      = 1000
integer, parameter :: max_pfts         = 20    ! pfts
integer, parameter :: max_pftps        = 100   ! pft parameters
integer, parameter :: max_years        = 1000
integer, parameter :: max_age          = 300

integer, parameter :: max_sites        = 100000

integer, parameter :: lai_comp_length  = 10
integer, parameter :: max_lai_comps    = 400

integer, parameter :: stem_comp_length = 30
integer, parameter :: max_stem_comps   = 20

integer, parameter :: root_comp_length = 30
integer, parameter :: max_root_comps   = 20

integer, parameter :: suma_comp_length = 30
integer, parameter :: max_suma_comps   = 20

integer, parameter :: max_outputs      = 75
integer, parameter :: str_len          = 1000

integer, parameter :: max_par_loops    = 11

end module dims
