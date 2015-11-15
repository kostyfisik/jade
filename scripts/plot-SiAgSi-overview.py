#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
def load_data(fname):
    data = np.loadtxt(fname+".dat")
    space = [np.nan]*len(data[0,:])
    #space = [-100]*len(data[0,:])
    min_value = np.min(data, axis=0)
    max_value = np.max(data, axis=0)
    range_value = max_value-min_value
    max_step = range_value/20.0
    data_spaced = data[0:1,:]
    for i in xrange(2, len(data[:,0])+1):
        diff = np.absolute(data[i-1:i,:] - data[i-2:i-1,:])
        need_space = False
        for j in xrange(2,len(max_step)):
            if j > 4:
                continue
            if max_step[j]<diff[0,j]:
                need_space = True
        if need_space:
            data_spaced = np.concatenate((data_spaced,[space]))
        data_spaced = np.concatenate((data_spaced,data[i-1:i,:]))
    return data, data_spaced

fname = "overview-Qsca"
data, data_spaced = load_data(fname)
for i in xrange(1, len(data[:,1])-2):
    if data[i-2,1]<=data[i,1] and data[i+2,1]<=data[i,1]:
        print(data[i,:])

fname2 = "sweep-ch"
data2, data_spaced2 = load_data(fname2)

max1 = (2*1+1)/(2*np.power(2*np.pi*data[:,0]/500,2))
max2 = (2*2+1)/(2*np.power(2*np.pi*data[:,0]/500,2))



############################# Plotting ######################
import numpy.ma as ma
vals = ma.array(data_spaced)
mvals = ma.masked_where(np.nan in data_spaced, vals)

fig, axs = plt.subplots(3,figsize=(4,6), sharex=True)#, sharey=True)
NACS=1
Qsca=0
Design=2
for ax in axs:
    ax.locator_params(axis='y',nbins=4)
    # for label in ['left', 'right', 'top', 'bottom']:
    #     ax.spines[label].set_position(('outward',-1.3))
    #ax.tick_params(axis='x', pad=30)

plotwidth=2.0
cax = axs[NACS].plot(data_spaced2[:,0], data_spaced2[:,1], linewidth=plotwidth,
                     solid_joinstyle='round', solid_capstyle='round', color='black'
                     , label=r"$\tilde{a}_1$"
)
cax = axs[NACS].plot(data_spaced2[:,0], data_spaced2[:,2], linewidth=plotwidth/1.5,
                     solid_joinstyle='round', solid_capstyle='round', color='red'
                     , label=r"$\tilde{b}_1$"
)
cax = axs[NACS].plot(data_spaced2[:,0], data_spaced2[:,3], linewidth=plotwidth,
                     solid_joinstyle='round', solid_capstyle='round', color='green'
                     , label=r"$\tilde{a}_2$"
)
cax = axs[NACS].plot(data_spaced2[:,0], data_spaced2[:,4], linewidth=plotwidth/1.5,
                     solid_joinstyle='round', solid_capstyle='round', color='blue'
                     , label=r"$\tilde{b}_2$"
)
axs[NACS].axhline(y=0.25, ls='--', dashes=[2,2], color='gray')
lg=axs[NACS].legend(loc='center left',prop={'size':11})
#lg=axs[Qsca].legend(loc='upper right',prop={'size':8})
#lg.get_frame().set_linewidth(0.0)
axs[NACS].annotate('0.25', xy=(27, 0.25), fontsize=9, color='gray',
                horizontalalignment='left', verticalalignment='bottom')

lg.draw_frame(False)



cax = axs[Qsca].plot(data_spaced[:,0], data_spaced[:,1], linewidth=plotwidth,
                     solid_joinstyle='round', solid_capstyle='round', color='black',
                     label=r"Si/Ag/Si")
#Analyic
cax = axs[Qsca].plot(data[:,0], max1, '--',linewidth=plotwidth/2.0,
                     solid_joinstyle='round', solid_capstyle='round', color='red'
#                     , label="max(n=1)"
                     )
dashes = [5, 2] # points on, off, ...
cax[0].set_dashes(dashes)
cax = axs[Qsca].plot(data[:,0], max2, '--',linewidth=plotwidth/2.0, color='blue'
#                     , label="max(n=2)"
                     )
dashes = [2, 2] # points on, off, ...
cax[0].set_dashes(dashes)
lg=axs[Qsca].legend(loc='upper left',prop={'size':10})
axs[Qsca].text(55, 1.2, r'max($n=1$)', fontsize=10, color='red')
axs[Qsca].text(55, 5.9, r'max($n=2$)', fontsize=10, color='blue')
#lg=axs[Qsca].legend(loc='upper right',prop={'size':8})
lg.draw_frame(False)
axs[Qsca].arrow(36, 4.5, 0, 0.6, head_width=1.2, head_length=0.3, fc='k', ec='k')
axs[Qsca].arrow(62.6, 3.35, 0, 0.6, head_width=1.2, head_length=0.3, fc='k', ec='k')
axs[Qsca].arrow(81.4, 4.3, 0, -0.6, head_width=1.2, head_length=0.3, fc='k', ec='k')


cax = axs[Design].plot(data_spaced[:,0], data_spaced[:,4], linewidth=plotwidth,
                       solid_joinstyle='round', solid_capstyle='round', color='black'
                       , label="outer shell"
                       )
cax = axs[Design].plot(data_spaced[:,0], data_spaced[:,3], linewidth=plotwidth,
                       solid_joinstyle='round', solid_capstyle='round', color='green'
                       , label="inner shell"
                       )
cax = axs[Design].plot(data_spaced[:,0], data_spaced[:,2], linewidth=plotwidth/1.5,
                       solid_joinstyle='round', solid_capstyle='round', color='red'
                       , label="core"
                       )
lg=axs[Design].legend(loc='upper left',prop={'size':10})
lg.draw_frame(False)

axs[NACS].set_ylabel(r'$\tilde{a}_n ,\ \tilde{b}_n$', labelpad=-0.9)
axs[NACS].set_ylim(0, 0.29)

axs[Qsca].set_ylabel(r'$Q_{abs}$', labelpad=8.8)
axs[Qsca].set_ylim(0, 7)
axs[Design].set_ylabel('Width, nm', labelpad=2)
axs[Design].set_ylim(0, 75)
axs[Design].set_xlabel(r'$R$, nm', labelpad=2)
plt.xlim(0,  89)
#plt.xlim(0,  160)
axs[NACS].annotate('(b)', xy=(0.99, 0.985), xycoords='axes fraction', fontsize=10,
                horizontalalignment='right', verticalalignment='top')
axs[Qsca].annotate('(a)', xy=(0.99, 0.985), xycoords='axes fraction', fontsize=10,
                horizontalalignment='right', verticalalignment='top')
axs[Design].annotate('(c)', xy=(0.99, 0.985), xycoords='axes fraction', fontsize=10,
                horizontalalignment='right', verticalalignment='top')
axs[Design].locator_params(axis='x',nbins=5)

fig.subplots_adjust(hspace=.05)

plt.savefig(fname+".pdf",pad_inches=0.02, bbox_inches='tight')
#plt.draw()

#plt.show()

plt.clf()
plt.close()

# cax = axs[0,0].imshow(Eabs_data, interpolation = 'nearest', cmap = cm.jet,
#                       origin = 'lower'
#                       #, vmin = min_tick, vmax = max_tick
#                       , extent = (min(scale_x), max(scale_x), min(scale_z), max(scale_z))
#                       #,norm = LogNorm()
#                       )


