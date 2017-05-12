#!/bin/sh

#installing grpc

cd /vagrant

echo "installing git"
sudo apt-get install git

#grpc
sudo apt-get install -y build-essential autoconf libtool libgflags-dev libgtest-dev clang libc++-dev cmake pkg-config

#protobuf
sudo apt-get install -y autoconf automake libtool curl make g++ zip unzip

#hypre
sudo apt-get install -y mpich

rm -r grpc
rm -r hypre

#get grpc
git clone -b $(curl -L http://grpc.io/release) https://github.com/grpc/grpc
cd grpc
git checkout master
git submodule update --init

#change protobuf version
cd third_party/protobuf
git checkout v3.0.x

#install protobuf https://github.com/google/protobuf/blob/master/src/README.md
./autogen.sh
#./configure --prefix=/usr
./configure
make
make check
sudo make install
sudo ldconfig #/usr/local/lib

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
make clean
make

cd /vagrant
