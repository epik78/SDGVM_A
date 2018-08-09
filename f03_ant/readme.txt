* data.f        : major changes
hydrology.f   : added 'sl' argument sublimation
metdos.f      : no change
sdgvm0.f
weathergenerator.f : added mean preserving interpolation Rymes & Myers 2001, and SWR
doly.f
light.f
nppcalc.f
sdgvm1.f
func.f : added 'douts' parameter to set lenght of string
luna.f
parameter_adjustment.f : no change
soil.f : added w_scalar and t_scalar
growth.f
luna_reorg.f
phenology.f
sunshade.f


*******************************
sdgvm

altered 'co2' from yearly to daily now co2(years,12,31)

*******************************
DATA.f90

subroutine readco2
added the ability to read daily co2 values

subroutine landuse


**************************************
