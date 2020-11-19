cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)
project(reentry_calc Fortran)

set (EXEC a.out)
# set (SRCDIR src)
message(STATUS "CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")

set(CMAKE_Fortran_FLAGS_DEBUG "-g -cpp -D_debug")

add_compile_options(
    -O2
    )

add_executable( ${EXEC}
    #${CMAKE_SOURCE_DIR}/src/mod_vector.f90
    #${CMAKE_SOURCE_DIR}/src/calcutils.f90
    ${CMAKE_SOURCE_DIR}/src/reentry_calc.f90
    )

