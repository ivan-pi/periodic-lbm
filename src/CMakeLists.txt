set(SRC
	precision.F90
	interp.f90
	fvm_bardow.F90
	periodic_dugks.F90
	periodic_lbm.f90
	collision_bgk.F90
	collision_bgk_improved.f90
	collision_trt.F90
	collision_regularized.F90
	output/vtk.F90
	output/gnuplot.F90
	benchmarks/barotropic_vortex_case.F90
	benchmarks/shear_layer.f90
	benchmarks/taylor_green.f90
	)

# Process C++ files independently
set(CXX_SRC
	output/vtk_polydata.cc)

add_library(${PROJECT_NAME}_cxx ${CXX_SRC})

add_library(${PROJECT_NAME} ${SRC})
target_link_libraries(${PROJECT_NAME}
	PUBLIC ${PROJECT_NAME}_cxx)

if(OpenMP_Fortran_FOUND)
	target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran)
endif()