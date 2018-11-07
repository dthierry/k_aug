#!/bin/bash

# Set the following to the latest MUMPS version.
#  THERE MUST BE NO SPACE BEFORE AND AFTER THE EQUAL (=) OPERATOR.
mumps_ver=5.1.2

set -e

wgetcmd=wget
wgetcount=`which wget 2>/dev/null | wc -l`
if test ! $wgetcount = 1; then
  echo "Utility wget not found in your PATH."
  if test `uname` = Darwin; then
    wgetcmd=ftp
    echo "Using ftp command instead."
  elif test `uname` = FreeBSD; then
    wgetcmd=fetch
    echo "Using fetch command instead."
  else
    exit -1
  fi
fi

echo " "
echo "Running script for downloading the source code for MUMPS"
echo " "

rm -f MUMPS*.tgz

echo "Downloading the source code from ..."
$wgetcmd http://mumps.enseeiht.fr/MUMPS_${mumps_ver}.tar.gz

echo "Uncompressing the tarball..."
gunzip -f MUMPS_${mumps_ver}.tar.gz

echo "Unpacking the source code..."
tar xf MUMPS_${mumps_ver}.tar

echo "Deleting the tar file..."
rm MUMPS_${mumps_ver}.tar

rm -rf MUMPS
mv MUMPS_${mumps_ver} MUMPS


cd MUMPS
cp Make.inc/Makefile.debian.SEQ ./Makefile.inc
os=$(uname -s)
ar=$(uname -m)

if [[ "$(uname -s)" = CYGWIN* ]]; then
    echo "WINDOWS"
    sed -i '1s|^|current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))\n|' Makefile.inc 
    sed -i '1s|^|mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))\n|' Makefile.inc 
    sed -i 's|LSCOTCHDIR = /usr/lib|LSCOTCHDIR = $(current_dir)/../../scotch/scotch/lib|g' Makefile.inc
    sed -i 's|ISCOTCH   = -I/usr/include/scotch|ISCOTCH = -I$(current_dir)/../../scotch/scotch/include|g' Makefile.inc
    sed -i 's|LMETISDIR = /usr/lib|LMETISDIR = $(current_dir)/../../metis/metis/build//libmetis|g' Makefile.inc
    sed -i 's|IMETIS    = -I/usr/include/metis|IMETIS = -I$(current_dir)/../../metis/metis/include|g' Makefile.inc
    sed -i 's|LAPACK = -llapack|LAPACK = $(current_dir)/../../openblas/OpenBLAS/libopenblas.dll.a|g' Makefile.inc
    sed -i 's|LIBBLAS = -lblas|LIBBLAS = $(current_dir)/../../openblas/OpenBLAS/libopenblas.dll.a|g' Makefile.inc
elif [[ "$(uname -s)" = Linux* ]]; then
    echo "LINUX"
    sed -i '1s|^|current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))\n|' Makefile.inc 
    sed -i '1s|^|mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))\n|' Makefile.inc 
    sed -i 's|LSCOTCHDIR = /usr/lib|LSCOTCHDIR = $(current_dir)/../../scotch/scotch/lib|g' Makefile.inc
    sed -i 's|ISCOTCH   = -I/usr/include/scotch|ISCOTCH = -I$(current_dir)/../../scotch/scotch/include|g' Makefile.inc
    sed -i 's|LMETISDIR = /usr/lib|LMETISDIR = $(current_dir)/../../metis/metis/build/Linux-x86_64/libmetis|g' Makefile.inc
    sed -i 's|IMETIS    = -I/usr/include/metis|IMETIS = -I$(current_dir)/../../metis/metis/include|g' Makefile.inc
    sed -i 's|LAPACK = -llapack|LAPACK = $(current_dir)/../../openblas/OpenBLAS/libopenblas.a|g' Makefile.inc
    sed -i 's|LIBBLAS = -lblas|LIBBLAS = $(current_dir)/../../openblas/OpenBLAS/libopenblas.a|g' Makefile.inc
elif [[ "$(uname -s)" = Darwin* ]]; then
    echo "DARWIN"
    sed -i '1s|^|current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))\n|' Makefile.inc 
    sed -i '1s|^|mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))\n|' Makefile.inc 
    sed -i 's|LSCOTCHDIR = /usr/lib|LSCOTCHDIR = $(current_dir)/../../scotch/scotch/lib|g' Makefile.inc
    sed -i 's|ISCOTCH   = -I/usr/include/scotch|ISCOTCH = -I$(current_dir)/../../scotch/scotch/include|g' Makefile.inc
    sed -i 's|LMETISDIR = /usr/lib|LMETISDIR = $(current_dir)/../../metis/metis/build/Linux-x86_64/libmetis|g' Makefile.inc
    sed -i 's|IMETIS    = -I/usr/include/metis|IMETIS = -I$(current_dir)/../../metis/metis/include|g' Makefile.inc
    sed -i 's|LAPACK = -llapack|LAPACK = $(current_dir)/../../openblas/OpenBLAS/libopenblas.a|g' Makefile.inc
    sed -i 's|LIBBLAS = -lblas|LIBBLAS = $(current_dir)/../../openblas/OpenBLAS/libopenblas.a|g' Makefile.inc

else
    echo "NO SYSTEM"
fi

make
echo "
█▀▄▀█ █░░█ █▀▄▀█ █▀▀█ █▀▀ ░ ░   █▀▀█ █▀▀ █▀▀█ █▀▀▄ █░░█
█░▀░█ █░░█ █░▀░█ █░░█ ▀▀█ ▄ ▄   █▄▄▀ █▀▀ █▄▄█ █░░█ █▄▄█
▀░░░▀ ░▀▀▀ ▀░░░▀ █▀▀▀ ▀▀▀ ░ █   ▀░▀▀ ▀▀▀ ▀░░▀ ▀▀▀░ ▄▄▄█
"

