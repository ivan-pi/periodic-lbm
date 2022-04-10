cmake_minimum_required(VERSION 3.14.0)

project(plbm 
		LANGUAGES Fortran C CXX
		DESCRIPTION "Collection of periodic LBM subroutines")


find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
  message(STATUS "Compiling with OpenMP")
endif()


find_package(PkgConfig REQUIRED)
if(PkgConfig_FOUND)
  message(STATUS "Using PkgConfig")
#  pkg_check_modules(MKL REQUIRED IMPORTED_TARGET GLOBAL mkl-dynamic-lp64-iomp)
endif()


find_package(Boost 1.65)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
	message(STATUS "Boost was NOT found")
endif()

add_subdirectory(src)
add_subdirectory(app)