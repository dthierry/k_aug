#!/bin/bash
echo "If os=mac pls edit this file"
# activate the following line if os=mac, deactivate wget line
# curl "http://ampl.com/netlib/ampl/solvers.tgz" -o solvers.tgz

rm solvers.*
rm -rf solvers
wget http://ampl.com/netlib/ampl/solvers.tgz
tar -zxvf solvers.tgz
echo "Making .."
cd ./solvers/
make -f makefile.u


echo "
    ___o .--.
   /___| |OO|
  /'   |_|  |_
       (_    _)
       | |   \
       | |oo_/sjw \n\n done "

