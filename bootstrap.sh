#!/bin/sh

#installing grpc

cd /vagrant
sudo apt-get install build-essential autoconf libtool
sudo apt-get install libgflags-dev libgtest-dev
sudo apt-get install clang libc++-dev
sudo apt-get install cmake
git clone -b $(curl -L http://grpc.io/release) \\
https://github.com/grpc/grpc
cd grpc
git submodule update --init
make && sudo make install
cd third_party/protobuf
sudo apt-get install unzip
./autogen.sh
./configure
make && sudo make install
cd -
cd examples/cpp/helloworld/
sudo apt-get install pkg-config
make

#installing hypre

cd /vagrant
sudo apt-get install mpich
git clone https://github.com/LLNL/hypre
cd hypre/src
./configure
make