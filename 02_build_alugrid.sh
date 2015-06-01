#!/bin/bash
pushd ALUGrid*
./configure --prefix=$(pwd)/../alugrid CXXFLAGS="-fPIC -DNDEBUG"
make
make install
popd

