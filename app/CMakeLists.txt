
macro(ADDEXEC name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} "${PROJECT_NAME}")
	set_property(TARGET ${name}
				 PROPERTY LINKER_LANGUAGE Fortran)
endmacro()

# We need to include the directory that contains
# module files. See the link below for an explanation:
#
#   https://stackoverflow.com/questions/43917443/fortran-module-files-not-found-by-cmake
#
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../src)

ADDEXEC(main_vortex)
ADDEXEC(main_taylor_green)

include_directories(/home/ivan/local/librsb_gnu/include)

ADDEXEC(main_rsb)
target_link_directories(main_rsb PUBLIC /home/ivan/local/librsb_gnu/lib)
target_link_libraries(main_rsb rsb)