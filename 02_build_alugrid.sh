#!/bin/bash
pushd ALUGrid*
./configure --prefix=$(pwd)/../alugrid CXXFLAGS="-fPIC -DNDEBUG -O2"
make
make install
popd

