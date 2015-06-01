#!/bin/bash

./dune-common-2.3.1/bin/dunecontrol --opts=clang.opts --module=dune-burgers all

# rebuild dune-burgers with --enable-shared

touch dune-burgers/src/dune_burgers.cc
./dune-common-2.3.1/bin/dunecontrol --opts=dune-burgers-clang.opts --only=dune-burgers all

