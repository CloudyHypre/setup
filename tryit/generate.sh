#!/bin/bash

protoc --proto_path=src/protos --cpp_out=src/cpp/ src/protos/hypre.proto
