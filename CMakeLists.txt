cmake_minimum_required(VERSION 3.5)
project(k_aug C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-g -Wall -ansi -pedantic-errors -Wall -W")
#set( PROJECT_LINK_LIBS amplsolver.a )

set(k_aug_version_mayor 1)
set(k_aug_version_minor 0)

set(USE_PARDISO 0)
set(USE_MC30 0)
set(PRINT_VERBOSE 0)

configure_file("${PROJECT_SOURCE_DIR}/src/common/config_kaug.h.in"
        "${PROJECT_BINARY_DIR}/config_kaug.h")

include_directories(${PROJECT_BINARY_DIR} thirdparty/asl/solvers
        thirdparty/mumps/MUMPS/include thirdparty/mumps/MUMPS/libseq)

link_directories(thirdparty/asl/solvers thirdparty/metis/metis/build/Linux-x86_64/libmetis thirdparty/scotch/scotch/lib
        thirdparty/mumps/MUMPS/lib thirdparty/openblas/OpenBLAS thirdparty/mumps/MUMPS/PORD/lib thirdparty/mumps/MUMPS/libseq
        )

if (USE_MC30)
    link_directories(thirdparty/hsl/mc30/mc30/src)
    message("INCLUDE LINK DIRECTORIES MC30 ${USE_MC30}")
else(USE_MC30)
    link_directories(thirdparty/hsl/mc19/mc19/src)
endif(USE_MC30)

set(SRC_FILES ${PROJECT_SOURCE_DIR}/src/k_aug/main_.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/assemble_rhs_dcdp.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/assemble_rhs_rh.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/csr_driver.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/dot_driver.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/find_inequalities.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/get_jac_asl_aug.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/get_hess_asl_aug.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/inertia_strategy.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/mu_adjust_primal.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/sigma_compute.c
        ${PROJECT_SOURCE_DIR}/src/k_aug/suffix_decl_hand.c
        ${PROJECT_SOURCE_DIR}/src/interfaces/hsl/mc19_driver.c
        ${PROJECT_SOURCE_DIR}/src/interfaces/mumps/mumps_driver.c
        ${PROJECT_SOURCE_DIR}/src/matrix/dsyev_driver.c
        ${PROJECT_SOURCE_DIR}/src/matrix/dpotri_driver.c
        )


if (USE_MC30)
    file(GLOB SOURCES "src/k_aug/*.c" "src/matrix/*.c" "src/interfaces/mumps/*.c" "src/HSL/*.c")
    message("INCLUDE LINK DIRECTORIES MC30 ${USE_MC30}")
else(USE_MC30)
    #file(GLOB SOURCES "src/k_aug/*.c" "src/k_aug/dot_driver/dot_driver.c" "src/matrix/*.c" "src/interfaces/mumps/*.c" "src/interfaces/hsl/*.c")
    #file(SOURCES SRC_FILES)
endif(USE_MC30)

add_executable(k_aug ${SRC_FILES})
#add_executable(dot_sens "src/k_aug/dot_driver/dot_driver.c")
add_executable(dot_sens src/k_aug/dot_driver/dot_driver.c)
#add_library(amplsolver STATIC IMPORTED)
#set_target_properties(amplsolver PROPERTIES IMPORTED_LOCATION thirdparty/asl/amplsolver.a)


if (USE_MC30)
    target_link_libraries(k_aug amplsolver m dl openblas gfortran dmumps mumps_common pthread esmumps scotch scotcherr
            metis pord mpiseq mc30 gfortran)
    message("LINKING MC30 ${USE_MC30}")
else(USE_MC30)
    target_link_libraries(k_aug amplsolver m dl openblas gfortran dmumps mumps_common pthread esmumps scotch scotcherr
            metis pord mpiseq mc19 gfortran)

endif(USE_MC30)

target_link_libraries(dot_sens amplsolver m dl openblas gfortran)

set_target_properties( k_aug
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

#message("${USE_MC30}")

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
foreach (var ${SRC_FILES})
    message("YE " ${SRC_FILES})
endforeach()