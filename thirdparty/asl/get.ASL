#!/bin/bash

set -e

#check for wget command
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
echo "Running script for downloading the source code for ASL"
echo " "

rm -f solvers.*
rm -rf solvers
$wgetcmd http://ampl.com/netlib/ampl/solvers.tgz
tar -zxvf solvers.tgz
echo "Making .."
cd ./solvers/
echo "#include <stddef.h>" > temp
echo "#include <sys/types.h>" >> temp
cat asl.h >> temp
mv temp asl.h

if [[ "$(uname -s)" = CYGWIN* ]]; then
    sed -i 's|a.out|a.exe|g' makefile.u
    echo "WINDOWS"
elif [[ "$(uname -s)" = Linux* ]]; then
    echo "LINUX"
else
    echo "NO SYSTEM"
fi



make -f makefile.u
cp amplsolver.a libamplsolver.a
cd ..
rm solvers.tgz

echo "
█▀▀█ █▀▀ █░░   █▀▀█ █▀▀ █▀▀█ █▀▀▄ █░░█
█▄▄█ ▀▀█ █░░   █▄▄▀ █▀▀ █▄▄█ █░░█ █▄▄█
▀░░▀ ▀▀▀ ▀▀▀   ▀░▀▀ ▀▀▀ ▀░░▀ ▀▀▀░ ▄▄▄█
"


