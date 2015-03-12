#!/bin/bash
rm -rf bin build
JADE_MPI_size=4 ./go.sh new |tee  run.log
mv run.log bin
