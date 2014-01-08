#!/bin/bash

# find no. of threads
cpus=`grep processor /proc/cpuinfo | wc -l`

# start maximum number of jobs in parallel
cat "$1" | xargs -L 1 -P $cpus ./superMC.e
