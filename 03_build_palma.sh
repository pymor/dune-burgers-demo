#!/bin/bash

./dune-common/bin/dunecontrol --opts=palma.opts --module=dune-burgers all

# rebuild dune-burgers with --enable-shared

touch dune-burgers/src/dune_burgers.cc
./dune-common/bin/dunecontrol --opts=dune-burgers-palma.opts --only=dune-burgers all

