#!/bin/bash

protoc --proto_path=src/protos --cpp_out=src/cpp/ src/protos/hypre.proto
protoc --proto_path=src/protos --grpc_out=src/cpp/ --plugin=protoc-gen-grpc=`which grpc_cpp_plugin` src/protos/hypre.proto
