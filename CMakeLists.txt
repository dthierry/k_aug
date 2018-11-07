cmake_minimum_required(VERSION 3.5)
project(k_aug C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-g -Wall -ansi -pedantic-errors -Wall -W")
#set( PROJECT_LINK_LIBS amplsolver.a )
execute_process(COMMAND uname -s OUTPUT_VARIABLE OS)

MESSAGE( STATUS "OS: " ${OS} )
execute_process(COMMAND uname -m OUTPUT_VARIABLE AR)
MESSAGE( STATUS "AR: " ${AR} )
set(newvar ${OS}-${AR})
message(STATUS "newvar :" ${newvar})
string(REGEX REPLACE "\n$" "" OS2 "${OS}")
string(REGEX REPLACE "\n$" "" AR2 "${AR}")
set(newvar2 ${OS2}-${AR2})
message(STATUS "newvar2 :" ${newvar2})
set(k_aug_version_mayor 1)
set(k_aug_version_minor 0)

set(USE_PARDISO 0)
set(USE_MC30 0)
set(PRINT_VERBOSE 0)

configure_file("${PROJECT_SOURCE_DIR}/src/common/config_kaug.h.in"
        "${PROJECT_BINARY_DIR}/config_kaug.h")

include_directories(${PROJECT_BINARY_DIR} thirdparty/asl/solvers
        thirdparty/mumps/MUMPS/include thirdparty/mumps/MUMPS/libseq)

#link_directories(thirdparty/asl/solvers thirdparty/metis/metis/build/${newvar2}/libmetis thirdparty/scotch/scotch/lib
#        thirdparty/mumps/MUMPS/lib thirdparty/openblas/OpenBLAS thirdparty/mumps/MUMPS/PORD/lib thirdparty/mumps/MUMPS/libseq
#        )
link_directories(thirdparty/asl/solvers thirdparty/metis/metis/build/CYGWIN_NT-6.1-x86_64/libmetis thirdparty/scotch/scotch/lib
        thirdparty/mumps/MUMPS/lib thirdparty/openblas/OpenBLAS thirdparty/mumps/MUMPS/PORD/lib thirdparty/mumps/MUMPS/libseq
        )



#link_directories("thirdparty/metis/metis/build/${newvar2}/libmetis")
if (USE_MC30)
    link_directories(thirdparty/hsl/mc30/mc30/src)
    message("INCLUDE LINK DIRECTORIES MC30 ${USE_MC30}")
else(USE_MC30)
    link_directories(thirdparty/hsl/mc19/mc19/src)
endif(USE_MC30)


if (USE_MC30)
    file(GLOB SOURCES "src/k_aug/*.c" "src/matrix/*.c" "src/interfaces/mumps/*.c" "src/HSL/*.c")
    message("INCLUDE LINK DIRECTORIES MC30 ${USE_MC30}")
else(USE_MC30)
    file(GLOB SOURCES "src/k_aug/*.c" "src/matrix/*.c" "src/interfaces/mumps/*.c" "src/interfaces/hsl/*.c")
endif(USE_MC30)

get_property(test_LINK_DIRECTORIES DIRECTORY PROPERTY LINK_DIRECTORIES)
message(STATUS " :" ${test_LINK_DIRECTORIES})
add_executable(k_aug ${SOURCES})


if (USE_MC30)
    target_link_libraries(k_aug amplsolver m dl openblas gfortran dmumps mumps_common pthread esmumps scotch scotcherr
            metis pord mpiseq mc30 gfortran)
    message("LINKING MC30 ${USE_MC30}")
else(USE_MC30)
    target_link_libraries(k_aug amplsolver m dl openblas gfortran dmumps mumps_common pthread esmumps scotch scotcherr
            metis pord mpiseq mc19 gfortran)

endif(USE_MC30)


set_target_properties( k_aug
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

message("${USE_MC30}")

