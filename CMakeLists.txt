cmake_minimum_required(VERSION 3.5)
project(k_aug C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-g -Wall -ansi -pedantic-errors -Wall -W")
message(STATUS "\nWelcome to k_aug.\n\n I hope you have already compiled the dependencies.\n\n")
message(STATUS "\nIf you want to have all the debug files generated each time, set PRINT_VERBOSE to 0 in this file.\n")
#set( PROJECT_LINK_LIBS amplsolver.a )
execute_process(COMMAND uname -s OUTPUT_VARIABLE OS)

#MESSAGE( STATUS "OS: " ${OS} )
execute_process(COMMAND uname -m OUTPUT_VARIABLE AR)
#MESSAGE( STATUS "AR: " ${AR} )
#set(newvar ${OS}-${AR})
#message(STATUS "newvar :" ${newvar})
string(REGEX REPLACE "\n$" "" OS2 "${OS}")
string(REGEX REPLACE "\n$" "" AR2 "${AR}")
set(newvar2 ${OS2}-${AR2})
message(STATUS "Your current OS and ARCH:" ${newvar2})
message(STATUS "Now k_aug assumes that you have all the coin-or libraries, i.e. coinhsl, coinmetis and coinasl")
message(STATUS "It is assumed that these libraries are located at /usr/local/lib if not please change the HINTS the CMakeLists.txt")

if (${OS} MATCHES "CYGWIN")
	message(STATUS "\n\nCYGWIN DETECTED:\n Make sure LAPACK is in the PATH after compilation")
  set(USE_CYGWIN 1)
endif (${OS} MATCHES "CYGWIN")

set(k_aug_version_mayor 1)
set(k_aug_version_minor 0)

set(USE_PARDISO 0)
set(USE_MC30 0)
set(PRINT_VERBOSE 1)

configure_file("${PROJECT_SOURCE_DIR}/src/common/config_kaug.h.in"
        "${PROJECT_BINARY_DIR}/config_kaug.h")

include_directories(${PROJECT_BINARY_DIR})

# coinor stuff
# PLEASE MAKE SURE THIS ARE IN THE CORRECT DIRECTORY
include_directories(/usr/local/include/coin-or/asl)
find_library(COINASL coinasl HINTS /usr/local/lib)
find_library(COINHSL coinhsl HINTS /usr/local/lib)
find_library(COINMETIS coinmetis HINTS /usr/local/lib)
# openblas/blas stuff
find_library(BLAS NAMES blas openblas HINTS /usr/lib /usr/local/lib)
find_library(LAPACK NAMES lapack openblas HINTS /usr/lib /usr/local/lib)


if (USE_MC30)
    link_directories(thirdparty/hsl/mc30/mc30/src)
    message("INCLUDE LINK DIRECTORIES MC30 ${USE_MC30}")
else(USE_MC30)
	#link_directories(thirdparty/hsl/mc19/mc19/src)
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
        ${PROJECT_SOURCE_DIR}/src/interfaces/hsl/ma57_driver.c
        ${PROJECT_SOURCE_DIR}/src/matrix/dsyev_driver.c
        ${PROJECT_SOURCE_DIR}/src/matrix/dpotri_driver.c
        )


add_executable(k_aug ${SRC_FILES})
add_executable(dot_sens src/k_aug/dot_driver/dot_driver.c)
if (USE_MC30)
    #target_link_libraries(k_aug amplsolver m dl gfortran dmumps mumps_common pthread esmumps scotch scotcherr metis pord mpiseq mc30 gfortran openblas ma57)
    target_link_libraries(k_aug mc30 ${COINMETIS} ${COINASL} ${COINHSL} ${BLAS} ${LAPACK} gfortran m)
    message("LINKING MC30 ${USE_MC30}")
else(USE_MC30)
    target_link_libraries(k_aug ${COINMETIS} ${COINASL} ${COINHSL} ${BLAS} ${LAPACK} gfortran m)
endif(USE_MC30)


target_link_libraries(dot_sens ${COINASL} ${BLAS})


set_target_properties( k_aug
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

message(STATUS "\n\nVery cool.\n\n Care for some Jams?\n\n \t\t https://www.youtube.com/watch?v=6GEI3PpXEAo \n\n Please run make.")



