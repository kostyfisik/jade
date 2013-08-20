#!/bin/bash
hg -v push ssh://onza-dev@onzafdtd.org//home/nfs-shared/onza-dev/jade++
echo Updating dir on onzafdtd.org
ssh onza-dev@onzafdtd.org "cd jade++ && hg update"
