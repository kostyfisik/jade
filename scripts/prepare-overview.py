#!/usr/bin/env python
import os
from glob import glob
WL_work = 3.75
os.system("rm R*-owerview.dat")            
for dirname, dirnames, filenames in os.walk('.'):
    # print path to all filenames.
    for filename in filenames:
        isGood = False
        if '-spectra.dat' in filename:
            # print filename
            sign = filename[0:-12]
            r =  sign[7:13]
            width = sign[23:28]
            rcs = sign[37:44]
            layers_num = 'n00'
            for num in [4,8,16,32,64]:
                if sign.find('-n'+str(num)+'-') > -1:
                    layers_num= 'R'+str(r)+'-n'+str(num)+'-owerview'
                    if not os.path.isfile(layers_num+'.dat'):
                        with open(layers_num+'.dat', "w") as myfile:
                            myfile.write('# '+layers_num+'\n')
                    with open(layers_num+'.dat', "a") as myfile:
                        myfile.write(width + '\t' + rcs + '\n')
                    # print width + '\t' + rcs
            # with open(filename) as f:
            #     f.readline()
            #     f.readline()
            #     array = [[float(x) for x in line.split()] for line in f]
            # test_line = array[0]
            # diff = WL_work
            # for line in array:
            #     new_diff = abs(line[0]-WL_work)
            #     if new_diff < diff:
            #         diff = new_diff
            #         test_line = line
            # work_index = array.index(test_line)
            # Qsca = [array[work_index-2][2],
            #         array[work_index-1][2],
            #         array[work_index][2],
            #         array[work_index+1][2],
            #         array[work_index+2][2]]
            # Qsca_sort = sorted(Qsca)
            # if Qsca == Qsca_sort:
            #     isGood = True
            #     # print Qsca
            #     # print filename
            # if Qsca == Qsca_sort.reverse():
            #     isGood = True
            #     print filename
            # if 2.0*(abs(Qsca[0]-Qsca[1]) + abs(Qsca[3]-Qsca[4])) > abs(Qsca[1]-Qsca[2]) + abs(Qsca[3]-Qsca[2]):
            #     isGood = True
            # # print Qsca
            # # print Qsca_sort
            # for f in glob ('*'+sign+'*'):
            #     if isGood != True:
            #         os.remove(f)
            # sign = filename[0:-12]
                # for line in f:
                # data = [int(s) for s in line.split()]
                # print data
                # if data[0] > 3.75:
                #     print data
                
#            print sign
#         print os.path.join(dirname, filename)
# print "Echo"

# target_size = 6
# for dirpath, dirs, files in os.walk('.'):
#     for file in files: 
#         path = os.path.join(dirpath, file)
#         if (os.stat(path).st_size == target_size) \
#                or (os.stat(path).st_size == target_size-1):
#             os.remove(path)
