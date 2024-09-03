
cmake_minimum_required(VERSION 3.15)

project(lbmex
		LANGUAGES Fortran C)

enable_language(CUDA OPTIONAL)

add_library(slbm  MODULE sim.F90 sim_slbm.F90 kinds.f90)
add_library(lw    MODULE sim.F90 sim_lw.F90 kinds.f90)
add_library(lw4   MODULE sim.F90 sim_lw4.F90 kinds.f90)
add_library(lw6   MODULE sim.F90 sim_lw6.F90 kinds.f90)
#add_library(dugks MODULE sim.f90 sim_dugks.f90)

add_library(fvm sim.F90 kinds.f90 sim_fvm.F90)
#target_link_options(fvm PUBLIC -acc)

add_executable(fvm_main
	fvm_main.F90
	sim_fvm.F90 sim.F90 kinds.f90 taylor_green.f90
	print_progress.c)


add_library(lwcuf MODULE sim.F90 kinds.f90 liblwcuf.F90)
target_compile_options(lwcuf PRIVATE -O3 -cuda -gpu=cc7.5,keep -stdpar=gpu -Minfo=all,stdpar -static-nvidia)
target_link_options(lwcuf PUBLIC -cuda -acc)

add_library(lw4cuf MODULE sim.F90 kinds.f90 sim_lw4cuf.F90)
target_compile_options(lw4cuf PRIVATE -O3 -cuda -gpu=cc7.5,keep -stdpar=gpu -Minfo=all,stdpar -static-nvidia)
target_link_options(lw4cuf PUBLIC -cuda -acc)


add_executable(standard_lbm 
	standard_lbm.F90 
	kinds.f90 
	sim.F90
	taylor_green.f90)