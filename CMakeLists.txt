cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)
project(reentry_calc Fortran)

set (EXEC a.out)
# set (SRCDIR src)
message(STATUS "CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")

set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -fbacktrace -D_debug")

add_compile_options(
    -O2 -fPIC -cpp 
    )

add_library(reentry_calc STATIC
	${CMAKE_SOURCE_DIR}/src/mod_vector.f90
	${CMAKE_SOURCE_DIR}/src/calcutils.f90
    ${CMAKE_SOURCE_DIR}/src/reentry_calc.f90
)

add_executable(${EXEC}
    ${CMAKE_SOURCE_DIR}/src/debug.f90
)
target_link_libraries(${EXEC} reentry_calc)
# add_executable( ${EXEC}
# 	${CMAKE_SOURCE_DIR}/src/mod_vector.f90
# 	${CMAKE_SOURCE_DIR}/src/calcutils.f90
#     ${CMAKE_SOURCE_DIR}/src/reentry_calc.f90
#     ${CMAKE_SOURCE_DIR}/src/debug.f90
#     )

