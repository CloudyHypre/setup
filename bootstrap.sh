#!/bin/sh

#installing grpc

cd /vagrant

sudo apt-get install git
#grpc
sudo apt-get install build-essential autoconf libtool
sudo apt-get install libgflags-dev libgtest-dev
sudo apt-get install clang libc++-dev
sudo apt-get install cmake
sudo apt-get install pkg-config
#protobuf
sudo apt-get install unzip
#hypre
sudo apt-get install mpich

#get grpc
git clone -b $(curl -L http://grpc.io/release) https://github.com/grpc/grpc
cd grpc
git checkout master
git submodule update --init

#change protobuf version
cd third_party/protobuf
git checkout v3.0.x

#install protobuf
./autogen.sh
./configure
make
sudo make install

#install grpc
cd -
make
sudo make install

#install example
cd examples/cpp/helloworld/
make

#install hypre
cd /vagrant
git clone https://github.com/LLNL/hypre
cd hypre/src
./configure
make

cd /vagrant/tryit/src/cpp/
make
cd /vagrant
