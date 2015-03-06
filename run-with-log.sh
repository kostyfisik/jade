#!/bin/bash
JADE_MPI_size=12 ./go.sh |tee  run.log
mv run.log bin
