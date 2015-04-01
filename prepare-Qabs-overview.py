#!/usr/bin/env python
import os
import string
from glob import glob
import numpy as np
##################################################################################################
##################################################################################################
##################################################################################################
#filename sould contain spectra columns: wl, Qext, Qsca, Qabs, Qbk
def find_max(filename, column=3): #default column=3 equals for Qabs
    wl_max = 0
    Q_max = 0
#    print filename
    if '-spectra.dat' in filename and 'bulk-Total' in filename:
        with open(filename) as f:
             f.readline()
             f.readline()
             array = []
             for line in f:
                 nums = line.split()
                 if len(nums)==5:
                     row = []
                     for x in nums:
                         if x.count('.') < 2:
                             row.append(float(x))
                         else:
                             row.append(0.0)
                     if row[column] > Q_max:
                         wl_max = row[0]
                         Q_max = row[column]
    return wl_max, Q_max        
##################################################################################################
##################################################################################################
##################################################################################################


#WL_work = 3.75
label = []
p = np.empty([0])
#os.system("rm -f R*-overview.dat")            
for dirname, dirnames, filenames in os.walk('.'):
    # print path to all filenames.
    isFirst = True
    for filename in filenames:
        if '-spectra.dat' in filename and 'bulk-Total' in filename:
        #if 'fails0-spectra.dat' in filename:
            # print filename
            sign = filename[0:-12]
            params = sign.split('-')
            print (sign)
            row = np.empty([0])
            for element in params:
                digits = element.translate(None,string.letters+'%')
                if len(digits) > 0:
                    #print element
                    value = float(digits)
                    row = np.append(row, value)
                    # if isFirst:
                    #     if (element[-2:] == 'nm'):
                    #         element = element[:-2]
                    #     name = element.translate(None,'.').translate(None,string.digits)
                    #     label.append(name)
                        
                        #print ("%s : %f"%(name, value))
            print(row)
            wl_max, Q_max = find_max(filename)
            #s = row[-4]+row[-2]+row[-3]
            #s = row[-3]
            #print row
            row = np.append(row,Q_max)
            row = np.append(row,wl_max)
            #row = np.append(row, row[1]*)
            p = np.append(p,row)
            #print p
            isFirst = False
p.shape = (len(p)/len(row),len(row))
label.append('TotalR')
label.append('Qabs')
label.append('epsRe')
label.append('epsIm')
label.append('fails')
label.append('Qabs_max')
label.append('WL_Qabs_max')

np.savetxt("overview-Qabs.txt",  p[p[:, 0].argsort()], fmt = "%.5f")
#print p
##################################################################################################
##################################################################################################
##################################################################################################

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass
mark = markers[3:3+len(label)]
##################################################################################################
#                    Qabs plot
##################################################################################################
if 'Qabs' in label:
    fig, ax1 = plt.subplots()
    ax1.grid(b=True, which='major', color='0.8', linestyle='-')
    ax1.set_axisbelow(True)
    ax1.plot(p[:,0], p[:,1], 'o', linewidth=2.0, ms=4, label=label[1], markeredgecolor='none')
    ax1.plot(p[:,0], p[:,-2], '^', linewidth=2.0, ms=4, label=label[-2], markeredgecolor='none')
    #fig.suptitle(label[1], fontsize=20)
    ax1.set_xlabel('%s, nm'%label[0], fontsize=16)
    ax1.set_ylabel(label[1], fontsize=16)
    # fig = plt.plot(p[:,0], p[:,1],'g^', linewidth=2.0)
    # fig.suptitle(label[1])
    # ax1.legend(loc='upper left')
    ax1.set_ylim([0,8])
    ax1.legend(loc='upper right')
    ax2 = ax1.twinx()
    ax2.plot(p[:,0], p[:,-1], linestyle='None', marker='x', ms=4, linewidth=3.0, label=label[-1], color='k', fillstyle='none')
    ax2.set_ylabel('%s, nm'%label[-1], fontsize=16)
    ax2.legend(loc=4)
    plt.savefig("overview-%s.svg"%label[1])
    #print p
##################################################################################################
#                    Layers width plot
##################################################################################################
    fig, ax1 = plt.subplots()
    mark = markers[3:3+len(label)]
    for plot_num in range(2,len(label)-2):
        plt.plot(p[:,0], p[:,plot_num], linestyle='None', ms=4, marker=mark[plot_num], linewidth=2.0, label=label[plot_num], markeredgecolor='none')
    #fig.suptitle(label[1], fontsize=20)
    # from matplotlib.font_manager import FontProperties
    # fontP = FontProperties()
    # fontP.set_size('small')
    plt.legend(loc='upper left')
    #plt.legend()
    plt.xlabel('%s, nm'%label[0], fontsize=16)
    plt.ylabel('Width, nm', fontsize=16)
    ax1.grid(b=True, which='major', color='0.8', linestyle='-')
    ax1.set_axisbelow(True)
    ax2 = ax1.twinx()
    ax2.plot(p[:,0], p[:,-1], linestyle='None', marker='x', ms=4, linewidth=3.0, label=label[-1], color='k', fillstyle='none')
    ax2.set_ylabel('%s, nm'%label[-1], fontsize=16)
    ax2.legend(loc=9)
    # fig = plt.plot(p[:,0], p[:,1],'g^', linewidth=2.0)
    # fig.suptitle(label[1])
    plt.savefig("overview-%s.svg"%label[2])
##################################################################################################
#                    Error and WL plot
##################################################################################################
    fig = plt.figure()
    for plot_num in range(1,3):
        plt.plot(p[:,0], p[:,-plot_num], linestyle='None', marker=mark[-plot_num], linewidth=2.0, label=label[-plot_num])
    #fig.suptitle(label[1], fontsize=20)
    # from matplotlib.font_manager import FontProperties
    # fontP = FontProperties()
    # fontP.set_size('small')
    plt.legend(loc='upper left')
    #plt.legend()
    plt.xlabel('%s,nm'%label[0], fontsize=16)
    plt.ylabel('%s,nm'%label[-1], fontsize=16)
    # fig = plt.plot(p[:,0], p[:,1],'g^', linewidth=2.0)
    # fig.suptitle(label[1])
    plt.savefig("overview-00-%s.svg"%label[-2])
else:
    fig = plt.figure()
    x=1
    y=2
    plt.plot(p[:,x], p[:,y], 'o', linewidth=2.0)
    #fig.suptitle(label[1], fontsize=20)
    plt.xlabel('%s'%label[x], fontsize=16)
    plt.ylabel(label[y], fontsize=16)
    # fig = plt.plot(p[:,0], p[:,1],'g^', linewidth=2.0)
    # fig.suptitle(label[1])
    plt.savefig("overview-%s.svg"%label[y])

# subarray = p[:,0]
# print subarray
#print(sort)

            # print (params)
            # print label
            # break
            # for num in [1,2,3,4,8,16,32,64]:
            #     if sign.find('-n'+str(num)+'-') > -1:
            #         layers_num= 'R'+str(r)+'-n'+str(num)+'-overview'
            #         if not os.path.isfile(layers_num+'.dat'):
            #             with open(layers_num+'.dat', "w") as myfile:
            #                 myfile.write('# '+layers_num+'\n')
            #         with open(layers_num+'.dat', "a") as myfile:
            #             myfile.write(width + '\t' + rcs + '\n')
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
