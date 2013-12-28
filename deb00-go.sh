#!/bin/bash
rsync -ave ssh ../jade++/ deb00:/home/tig/jade++ && ssh deb00 "cd ~/jade++ && ./go.sh"