# ======================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ======================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
#FORTRAN_COMPILER = gfortran
FC = gfortran -O3 -fbounds-check -fdefault-real-8 -fbacktrace -fcheck=all
#FC = gfortran -O3 -fcheck=all

s = src
m = mod
o = obj
b = bin
e = f90

all: \
 $(o)/a_dims.o \
 $(o)/a_real_precision.o \
 $(o)/a_screen_output.o \
 $(o)/file_io_methods.o \
 $(o)/a_tuning_parameters.o \
 $(o)/light_methods.o \
 $(o)/a_misc_values.o \
 $(o)/sunshade.o \
 $(o)/b_output_filenames.o \
 $(o)/a_site_parameters.o \
 $(o)/a_pft_parameters.o \
 $(o)/b_input_file.o \
 $(o)/func.o \
 $(o)/weather_generator.o \
 $(o)/a_system_state.o \
 $(o)/productivity_luna_methods.o \
 $(o)/luna_methods.o \
 $(o)/state_methods.o \
 $(o)/open_files.o \
 $(o)/data.o \
 $(o)/veg_dynamics.o \
 $(o)/crops.o \
 $(o)/b_output_methods.o \
 $(o)/b_input_methods.o \
 $(o)/hydrology_methods.o \
 $(o)/phenology_methods.o \
 $(o)/productivity_methods.o \
 $(o)/daily_step.o \
 $(o)/soil_methods.o \
 $(o)/sdgvm.o

	$(FC) -o $b/sdgvm.exe \
 $(o)/a_dims.o \
 $(o)/a_real_precision.o \
 $(o)/a_screen_output.o \
 $(o)/file_io_methods.o \
 $(o)/a_tuning_parameters.o \
 $(o)/light_methods.o \
 $(o)/a_misc_values.o \
 $(o)/sunshade.o \
 $(o)/b_output_filenames.o \
 $(o)/a_site_parameters.o \
 $(o)/a_pft_parameters.o \
 $(o)/b_input_file.o \
 $(o)/func.o \
 $(o)/weather_generator.o \
 $(o)/a_system_state.o \
 $(o)/productivity_luna_methods.o \
 $(o)/luna_methods.o \
 $(o)/state_methods.o \
 $(o)/open_files.o \
 $(o)/data.o \
 $(o)/veg_dynamics.o \
 $(o)/crops.o \
 $(o)/b_output_methods.o \
 $(o)/b_input_methods.o \
 $(o)/hydrology_methods.o \
 $(o)/phenology_methods.o \
 $(o)/productivity_methods.o \
 $(o)/daily_step.o \
 $(o)/soil_methods.o \
 $(o)/sdgvm.o -Jmod

$(o)/a_dims.o: $(s)/a_dims.$(e)
	$(FC) -o $(o)/a_dims.o -c $(s)/a_dims.$(e) -J$(m)

$(o)/a_real_precision.o: $(s)/a_real_precision.$(e) 
	$(FC) -o $(o)/a_real_precision.o -c $(s)/a_real_precision.$(e) -J$(m)

$(o)/a_screen_output.o: $(s)/a_screen_output.$(e) $(s)/a_real_precision.$(e)
	$(FC) -o $(o)/a_screen_output.o -c $(s)/a_screen_output.$(e) -J$(m)

$(o)/file_io_methods.o: $(s)/file_io_methods.$(e) $(s)/a_dims.$(e) 
	$(FC) -o $(o)/file_io_methods.o -c $(s)/file_io_methods.$(e) -J$(m)

$(o)/a_tuning_parameters.o: $(s)/a_tuning_parameters.$(e) $(s)/a_real_precision.$(e)
	$(FC) -o $(o)/a_tuning_parameters.o -c $(s)/a_tuning_parameters.$(e) -J$(m)

$(o)/light_methods.o: $(s)/light_methods.$(e) $(s)/a_real_precision.$(e)
	$(FC) -o $(o)/light_methods.o -c $(s)/light_methods.$(e) -J$(m)

$(o)/a_misc_values.o: $(s)/a_misc_values.$(e) $(s)/a_real_precision.$(e)
	$(FC) -o $(o)/a_misc_values.o -c $(s)/a_misc_values.$(e) -J$(m)

$(o)/sunshade.o: $(s)/sunshade.$(e) $(s)/a_real_precision.$(e)
	$(FC) -o $(o)/sunshade.o -c $(s)/sunshade.$(e) -J$(m)

$(o)/b_output_filenames.o: $(s)/b_output_filenames.$(e) $(s)/a_dims.$(e)
	$(FC) -o $(o)/b_output_filenames.o -c $(s)/b_output_filenames.$(e) -J$(m)

$(o)/a_site_parameters.o: $(s)/a_site_parameters.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e)
	$(FC) -o $(o)/a_site_parameters.o -c $(s)/a_site_parameters.$(e) -J$(m)

$(o)/a_pft_parameters.o: $(s)/a_pft_parameters.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e)
	$(FC) -o $(o)/a_pft_parameters.o -c $(s)/a_pft_parameters.$(e) -J$(m)

$(o)/b_input_file.o: $(s)/b_input_file.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e)
	$(FC) -o $(o)/b_input_file.o -c $(s)/b_input_file.$(e) -J$(m)

$(o)/func.o: $(s)/func.$(e) $(s)/a_real_precision.$(e) $(s)/a_pft_parameters.$(e) $(s)/a_dims.$(e)
	$(FC) -o $(o)/func.o -c $(s)/func.$(e) -J$(m)

$(o)/weather_generator.o: $(s)/weather_generator.$(e) $(s)/a_real_precision.$(e) $(s)/b_input_file.$(e) $(s)/func.$(e)
	$(FC) -o $(o)/weather_generator.o -c $(s)/weather_generator.$(e) -J$(m)

$(o)/a_system_state.o: $(s)/a_system_state.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e)
	$(FC) -o $(o)/a_system_state.o -c $(s)/a_system_state.$(e) -J$(m)

$(o)/productivity_luna_methods.o: $(s)/productivity_luna_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e) $(s)/a_system_state.$(e) $(s)/a_site_parameters.$(e) $(s)/a_pft_parameters.$(e)
	$(FC) -o $(o)/productivity_luna_methods.o -c $(s)/productivity_luna_methods.$(e) -J$(m)

$(o)/luna_methods.o: $(s)/luna_methods.$(e) $(s)/a_real_precision.$(e) $(s)/productivity_luna_methods.$(e)
	$(FC) -o $(o)/luna_methods.o -c $(s)/luna_methods.$(e) -J$(m)

$(o)/state_methods.o: $(s)/state_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e) $(s)/b_input_file.$(e) $(s)/a_system_state.$(e) $(s)/a_site_parameters.$(e) $(s)/a_pft_parameters.$(e) $(s)/a_misc_values.$(e) $(s)/a_tuning_parameters.$(e) $(s)/file_io_methods.$(e) $(s)/b_output_filenames.$(e)
	$(FC) -o $(o)/state_methods.o -c $(s)/state_methods.$(e) -J$(m)

$(o)/open_files.o: $(s)/open_files.$(e) $(s)/a_real_precision.$(e) $(s)/b_input_file.$(e) $(s)/file_io_methods.$(e) $(s)/func.$(e) $(s)/a_pft_parameters.$(e)
	$(FC) -o $(o)/open_files.o -c $(s)/open_files.$(e) -J$(m)

$(o)/data.o: $(s)/data.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e) $(s)/func.$(e) $(s)/b_input_file.$(e) $(s)/a_site_parameters.$(e) $(s)/a_system_state.$(e) 
	$(FC) -o $(o)/data.o -c $(s)/data.$(e) -J$(m)

$(o)/veg_dynamics.o: $(s)/veg_dynamics.$(e) $(s)/file_io_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e)  $(s)/a_system_state.$(e) $(s)/a_site_parameters.$(e) $(s)/state_methods.$(e) $(s)/func.$(e) $(s)/a_tuning_parameters.$(e)
	$(FC) -o $(o)/veg_dynamics.o -c $(s)/veg_dynamics.$(e) -J$(m)

$(o)/crops.o: $(s)/crops.$(e) $(s)/a_real_precision.$(e) $(s)/b_input_file.$(e) $(s)/a_pft_parameters.$(e) $(s)/a_site_parameters.$(e) $(s)/func.$(e) $(s)/light_methods.$(e) $(s)/a_system_state.$(e) $(s)/b_input_file.$(e)
	$(FC) -o $(o)/crops.o -c $(s)/crops.$(e) -J$(m)

$(o)/b_output_methods.o: $(s)/b_output_methods.$(e) $(s)/a_dims.$(e) $(s)/a_real_precision.$(e) $(s)/file_io_methods.$(e) $(s)/b_input_file.$(e) $(s)/data.$(e) $(s)/weather_generator.$(e) $(s)/veg_dynamics.$(e) $(s)/func.$(e) $(s)/weather_generator.$(e)
	$(FC) -o $(o)/b_output_methods.o -c $(s)/b_output_methods.$(e) -J$(m)

$(o)/b_input_methods.o: $(s)/b_input_methods.$(e) $(s)/file_io_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_pft_parameters.$(e) $(s)/a_site_parameters.$(e) $(s)/b_output_methods.$(e) $(s)/func.$(e) $(s)/a_screen_output.$(e) $(s)/a_tuning_parameters.$(e) $(s)/open_files.$(e)
	$(FC) -o $(o)/b_input_methods.o -c $(s)/b_input_methods.$(e) -J$(m)

$(o)/hydrology_methods.o: $(s)/hydrology_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_system_state.$(e) $(s)/a_site_parameters.$(e) $(s)/a_tuning_parameters.$(e) $(s)/a_pft_parameters.$(e) $(s)/b_input_file.$(e)
	$(FC) -o $(o)/hydrology_methods.o -c $(s)/hydrology_methods.$(e) -J$(m)

$(o)/phenology_methods.o: $(s)/phenology_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e) $(s)/a_system_state.$(e) $(s)/a_pft_parameters.$(e) $(s)/a_site_parameters.$(e) $(s)/a_tuning_parameters.$(e) $(s)/a_misc_values.$(e) $(s)/light_methods.$(e)
	$(FC) -o $(o)/phenology_methods.o -c $(s)/phenology_methods.$(e) -J$(m)

$(o)/productivity_methods.o: $(s)/productivity_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e) $(s)/sunshade.$(e) $(s)/a_system_state.$(e) $(s)/a_pft_parameters.$(e) $(s)/a_site_parameters.$(e) $(s)/a_tuning_parameters.$(e) $(s)/productivity_luna_methods.$(e) $(s)/luna_methods.$(e) $(s)/func.$(e) $(s)/light_methods.$(e)
	$(FC) -o $(o)/productivity_methods.o -c $(s)/productivity_methods.$(e) -J$(m)

$(o)/daily_step.o: $(s)/daily_step.$(e) $(s)/a_real_precision.$(e) $(s)/phenology_methods.$(e) $(s)/productivity_methods.$(e) $(s)/a_misc_values.$(e) $(s)/a_site_parameters.$(e) $(s)/a_tuning_parameters.$(e) $(s)/a_system_state.$(e)
	$(FC) -o $(o)/daily_step.o -c $(s)/daily_step.$(e) -J$(m)

$(o)/soil_methods.o: $(s)/soil_methods.$(e) $(s)/a_real_precision.$(e) $(s)/state_methods.$(e) $(s)/a_misc_values.$(e)  $(s)/a_tuning_parameters.$(e) $(s)/file_io_methods.$(e)
	$(FC) -o $(o)/soil_methods.o -c $(s)/soil_methods.$(e) -J$(m)

$(o)/sdgvm.o: $(s)/sdgvm.$(e) $(s)/file_io_methods.$(e) $(s)/a_real_precision.$(e) $(s)/a_dims.$(e) $(s)/b_input_file.$(e) $(s)/a_system_state.$(e) $(s)/a_pft_parameters.$(e) $(s)/state_methods.$(e) $(s)/b_output_methods.$(e) $(s)/b_input_methods.$(e) $(s)/phenology_methods.$(e) $(s)/daily_step.$(e) $(s)/hydrology_methods.$(e) $(s)/light_methods.$(e) $(s)/soil_methods.$(e) $(s)/a_misc_values.$(e) $(s)/crops.$(e)
	$(FC) -o $(o)/sdgvm.o -c $(s)/sdgvm.$(e) -J$(m)

clean:
	rm -f $b/sdgvm.exe $b/sdgvm $m/*.mod $o/*.o $s/*.obj

