!**********************************************************************!
!                                                                      !
!         screen_output_parameters :: screen_output_parameters         ! 
!                    ------------------------                          !
!                                                                      !
! module screen_output_parameters                                      !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Defines ScreenOutputParameters structures and declares
!! the sop variable.sop is assigned values in read_param function.
!! @details ScreenOutputParameters%site_out holds integer on how many 
!! grid cells will show on screen to show progress during execution.
!! If set to 1 then every grid cell will appear on screen.
!! ScreenOutputParameters%year_out holds how the period in years that
!! a grid cell will show on screen during execution.If set to 10 then
!! every 10 years the progress for each grid cell will show on screen.
!! ScreenOutputParameters%l_regional is logical and set to true for 
!! regional runs and false for country runs
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
module screen_output_parameters
!**********************************************************************!
use real_precision

type ScreenOutputParameters
  integer :: site_out
  integer :: year_out
  logical :: l_regional
end type

type (ScreenOutputParameters) sop

end module screen_output_parameters

