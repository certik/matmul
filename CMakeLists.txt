cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_SOURCE_DIR}/cmake/UserOverride.cmake)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "Build type (Debug, Release)")
endif()

project(matmul)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)

# Make sure that CMAKE_BUILD_TYPE is either Debug or Release:
if (NOT CMAKE_BUILD_TYPE MATCHES "Debug|Release")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

if (APPLE)
    set(DEFAULT_MATMUL_BLAS "Accelerate")
else()
    set(DEFAULT_MATMUL_BLAS "Fortran")
endif()
set(MATMUL_BLAS ${DEFAULT_MATMUL_BLAS}
    CACHE STRING "The BLAS library that Matmul should use")
if (NOT MATMUL_BLAS MATCHES "Accelerate|OpenBLAS|Fortran")
    message(FATAL_ERROR "MATMUL_BLAS must be one of: OpenBLAS, Accelerate, Fortran (current value: '${MATMUL_BLAS}')")
endif ()
if (MATMUL_BLAS STREQUAL "Accelerate")
    find_package(OMP)
    if(NOT OMP_FOUND)
        find_package(OpenMP REQUIRED COMPONENTS Fortran)
    endif()
elseif (MATMUL_BLAS STREQUAL "OpenBLAS")
    find_package(OPENBLAS REQUIRED)
    find_package(OMP)
    if(NOT OMP_FOUND)
        find_package(OpenMP REQUIRED COMPONENTS Fortran)
    endif()
else()
    # pass
endif()

enable_testing()

set(SRC
    m.f90
    )
if (MATMUL_BLAS STREQUAL "Accelerate")
    list(APPEND SRC
        linalg_accelerate.c
        linalg_c.f90
    )
elseif (MATMUL_BLAS STREQUAL "OpenBLAS")
    list(APPEND SRC
        linalg_openblas.c
        linalg_c.f90
    )
else()
    list(APPEND SRC
        linalg_f.f90
        omp_dummy.f90
    )
endif()
add_library(matmullib ${SRC})
if (MATMUL_BLAS STREQUAL "Accelerate")
    target_link_options(matmullib PUBLIC -framework accelerate)
    target_link_libraries(matmullib PUBLIC "$<IF:$<BOOL:${OMP_FOUND}>,p::omp,OpenMP::OpenMP_Fortran>")
elseif (MATMUL_BLAS STREQUAL "OpenBLAS")
    target_link_libraries(matmullib p::openblas
        "$<IF:$<BOOL:${OMP_FOUND}>,p::omp,OpenMP::OpenMP_Fortran>")
endif()

add_executable(matmul main.f90)
target_link_libraries(matmul matmullib)

message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message("MATMUL_BLAS: ${MATMUL_BLAS}")
