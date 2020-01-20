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
if (${OS} MATCHES "CYGWIN")
	message(STATUS "\n\nCYGWIN DETECTED:\n Make sure openblas is in the PATH after compilation")
endif (${OS} MATCHES "CYGWIN")

set(k_aug_version_mayor 1)
set(k_aug_version_minor 0)

set(USE_PARDISO 0)
set(USE_MC30 0)
set(PRINT_VERBOSE 1)

configure_file("${PROJECT_SOURCE_DIR}/src/common/config_kaug.h.in"
        "${PROJECT_BINARY_DIR}/config_kaug.h")

include_directories(${PROJECT_BINARY_DIR} thirdparty/asl/solvers)

link_directories(thirdparty/asl/solvers thirdparty/metis/metis/build/${newvar2}/libmetis /home/dav0/in_dev_/test_pynumero/miniconda/pkgs/libopenblas-0.3.6-h5a2b251_1/lib /home/dav0/in_dev_/test_pynumero/miniconda/pkgs/libgfortran-ng-7.3.0-hdf63c60_0/lib)


## Please modify this.
link_directories(/home/dav0/PycharmProjects/ipopt_van/build/lib)  ## THIS IS YOUR OWN ipopt build lib directory
## Make sure you have libcoinhsl.so and libcoinmetis.so!!



link_directories(thirdparty/asl/solvers thirdparty/metis/metis/build/${newvar2}/libmetis thirdparty/scotch/scotch/lib
        #thirdparty/mumps/MUMPS/lib thirdparty/openblas/OpenBLAS thirdparty/mumps/MUMPS/PORD/lib thirdparty/mumps/MUMPS/libseq thirdparty/hsl/ma57/build/src
        thirdparty/mumps/MUMPS/lib thirdparty/openblas/OpenBLAS thirdparty/mumps/MUMPS/PORD/lib thirdparty/mumps/MUMPS/libseq /home/dav0/build_ipopt/lib
        )
#link_directories(thirdparty/asl/solvers thirdparty/metis/metis/build/CYGWIN_NT-6.1-x86_64/libmetis thirdparty/scotch/scotch/lib
#        thirdparty/mumps/MUMPS/lib thirdparty/openblas/OpenBLAS thirdparty/mumps/MUMPS/PORD/lib thirdparty/mumps/MUMPS/libseq
#        )
#link_directories("thirdparty/metis/metis/build/${newvar2}/libmetis")

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
        ${PROJECT_SOURCE_DIR}/src/interfaces/hsl/ma57_driver.c
        ${PROJECT_SOURCE_DIR}/src/matrix/dsyev_driver.c
        ${PROJECT_SOURCE_DIR}/src/matrix/dpotri_driver.c
        )


    #if (USE_MC30)
    #file(GLOB SOURCES "src/k_aug/*.c" "src/matrix/*.c" "src/interfaces/mumps/*.c" "src/HSL/*.c")
    #message("INCLUDE LINK DIRECTORIES MC30 ${USE_MC30}")
    #else(USE_MC30)
    #file(GLOB SOURCES "src/k_aug/*.c" "src/k_aug/dot_driver/dot_driver.c" "src/matrix/*.c" "src/interfaces/mumps/*.c" "src/interfaces/hsl/*.c")
    #file(SOURCES SRC_FILES)
    #endif(USE_MC30)

add_executable(k_aug ${SRC_FILES})
#add_executable(dot_sens "src/k_aug/dot_driver/dot_driver.c")
add_executable(dot_sens src/k_aug/dot_driver/dot_driver.c)

if (USE_MC30)
    #target_link_libraries(k_aug amplsolver m dl gfortran dmumps mumps_common pthread esmumps scotch scotcherr metis pord mpiseq mc30 gfortran openblas ma57)
    target_link_libraries(k_aug amplsolver m dl gfortran dmumps mumps_common pthread esmumps scotch scotcherr metis pord mpiseq mc30 gfortran openblas coinhsl)
    message("LINKING MC30 ${USE_MC30}")
else(USE_MC30)
	# CYGWIN: are you kidding me?!
    #target_link_libraries(k_aug amplsolver m dl gfortran dmumps mumps_common pthread esmumps scotch scotcherr metis pord mpiseq mc19 gfortran openblas ma57)
    target_link_libraries(k_aug amplsolver m dl gfortran pthread mc19 gfortran coinhsl gfortran coinmetis)
endif(USE_MC30
	)

target_link_libraries(k_aug blas lapack)
target_link_libraries(dot_sens amplsolver m dl gfortran blas lapack)


set_target_properties( k_aug
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

message(STATUS "\n\nPlease run make.")



