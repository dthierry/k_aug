cmake_minimum_required(VERSION 3.10)
project(k_aug C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-g -Wall -ansi -pedantic-errors -Wall -W")
#set( PROJECT_LINK_LIBS amplsolver.a )

include_directories(thirdparty/asl/solvers thirdparty/mumps/MUMPS/include thirdparty/mumps/MUMPS/libseq)

link_directories(thirdparty/asl/solvers thirdparty/metis/metis/build/Linux-x86_64/libmetis thirdparty/scotch/scotch/lib
        thirdparty/mumps/MUMPS/lib thirdparty/openblas/OpenBLAS thirdparty/hsl/mc30/mc30/src
        thirdparty/mumps/MUMPS/PORD/lib thirdparty/mumps/MUMPS/libseq)


file(GLOB SOURCES "src/k_aug/*.c" "src/matrix/*.c" "src/interfaces/mumps/*.c" "src/HSL/*.c")
add_executable(k_aug
        ${SOURCES})
target_link_libraries(k_aug amplsolver m dl openblas mc30 gfortran dmumps mumps_common pthread esmumps scotch scotcherr
        metis pord mpiseq)