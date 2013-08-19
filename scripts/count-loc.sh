#!/bin/bash
cd ../src
find . -name '*.h' -o -name '*.cc'| xargs wc -l

