#!/usr/bin/env python
import os
from glob import glob
WL_work = 3.75
print "Filter is on for WL =", WL_work
for dirname, dirnames, filenames in os.walk('.'):
    # print path to all filenames.
    for filename in filenames:
        isGood = False
        if '-spectra.dat' in filename:
            print filename
            sign = filename[0:-12]
            with open(filename) as f:
                f.readline()
                f.readline()
                array = [[float(x) for x in line.split()] for line in f]
            test_line = array[0]
            diff = WL_work
            for line in array:
                new_diff = abs(line[0]-WL_work)
                if new_diff < diff:
                    diff = new_diff
                    test_line = line
            work_index = array.index(test_line)
            Qsca = [array[work_index-2][2],
                    array[work_index-1][2],
                    array[work_index][2],
                    array[work_index+1][2],
                    array[work_index+2][2]]
            Qsca_sort = sorted(Qsca)
            if Qsca == Qsca_sort:
                isGood = True
                # print Qsca
                # print filename
            if Qsca == Qsca_sort.reverse():
                isGood = True
                print filename
            if 2.0*(abs(Qsca[0]-Qsca[1]) + abs(Qsca[3]-Qsca[4])) > abs(Qsca[1]-Qsca[2]) + abs(Qsca[3]-Qsca[2]):
                isGood = True
            # print Qsca
            # print Qsca_sort
            for f in glob ('*'+sign+'*'):
                if isGood != True:
                    os.remove(f)
            # sign = filename[0:-12]
                # for line in f:
                # data = [int(s) for s in line.split()]
                # print data
                # if data[0] > 3.75:
                #     print data
                
#            print sign
#         print os.path.join(dirname, filename)
# print "Echo"
