import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data_path_5 = "./Data/mean_std_ceytometry_SMAF.csv"
SMAF = np.loadtxt(open(data_path_5),delimiter=",",skiprows=0)
#print(SMAF.shape)


data_path_3 = "./Data/DeepAE distictc measurements_ceytometry.csv"
data_measurement = np.loadtxt(open(data_path_3),delimiter=",",skiprows=0)


PCC_3 = data_measurement[:,(9,6,3,0)]
EM_3 = data_measurement[:,(10,7,4,1)]
MAE_3 = data_measurement[:,(11,8,5,2)]

ind = np.arange(4)  # the x locations for the groups
width = 0.15  # the width of the bars
color_list = plt.cm.Set3(np.linspace(0, 1, 12))

from matplotlib import rcParams
rcParams.update({'font.size': 14, 'font.family': 'STIXGeneral'})

from matplotlib import ticker
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((0,1))

fig, axes = plt.subplots(nrows=1,ncols=3, figsize=(15,4))


rects1 = axes[0].bar(ind - width/2 - 2*width, SMAF[(2,10,18,26),0], 
                       width, yerr=SMAF[(3,11,19,27),0], color=color_list[0], label='SVD')
rects2 = axes[0].bar(ind - width/2 - width, SMAF[(4,12,20,28),0],  
                       width, yerr=SMAF[(5,13,21,29),0], color=color_list[7], label='k-SVD')
rects3 = axes[0].bar(ind - width/2, SMAF[(6,14,22,30),0],  
                       width, yerr=SMAF[(7,15,23,31),0], color=color_list[2], label='sNMF')
rects4 = axes[0].bar(ind + width/2, SMAF[(0,8,16,24),0],  
                       width, yerr=SMAF[(1,9,17,25),0], color=color_list[3], label='CS-SMAF')
rects5 = axes[0].bar(ind + width/2 + width,   np.mean(PCC_3[0:5, :], axis=0),  
                       width, yerr=np.std(PCC_3[0:5, :], axis=0), color=color_list[4], label='DeepAE')
axes[0].set_ylabel('PCC')
axes[0].set_xticks(ind)
axes[0].set_ylim(0.4,1)
axes[0].set_xticklabels(('10', '25', '50', '100'))
#axes[0,0].set_title('GSE45234')


rects1 = axes[1].bar(ind - width/2 - 2*width, SMAF[(2,10,18,26),1], 
                       width, yerr=SMAF[(3,11,19,27),0], color=color_list[0], label='SVD')
rects2 = axes[1].bar(ind - width/2 - width, SMAF[(4,12,20,28),1],  
                       width, yerr=SMAF[(5,13,21,29),0], color=color_list[7], label='k-SVD')
rects3 = axes[1].bar(ind - width/2, SMAF[(6,14,22,30),1],  
                       width, yerr=SMAF[(7,15,23,31),0], color=color_list[2], label='sNMF')
rects4 = axes[1].bar(ind + width/2, SMAF[(0,8,16,24),1],  
                       width, yerr=SMAF[(1,9,17,25),1], color=color_list[3], label='CS-SMAF')
rects5 = axes[1].bar(ind + width/2 + width, np.mean(EM_3[0:5, :], axis=0),  
                       width, yerr=np.std(EM_3[0:5, :], axis=0), color=color_list[4], label='DeepAE')

axes[1].set_ylabel('EM')
axes[1].set_xticks(ind)
axes[1].set_xticklabels(('10', '25', '50', '100'))
axes[1].set_title('Mass cytometry data')
axes[1].yaxis.set_major_formatter(formatter) 
axes[1].set_xlabel('Measurements')


rects1 = axes[2].bar(ind - width/2 - 2*width, SMAF[(2,10,18,26),2], 
                       width, yerr=SMAF[(3,11,19,27),0], color=color_list[0], label='SVD')
rects2 = axes[2].bar(ind - width/2 - width, SMAF[(4,12,20,28),2],  
                       width, yerr=SMAF[(5,13,21,29),2], color=color_list[7], label='k-SVD')
rects3 = axes[2].bar(ind - width/2, SMAF[(6,14,22,30),2],  
                       width, yerr=SMAF[(7,15,23,31),0], color=color_list[2], label='sNMF')
rects4 = axes[2].bar(ind + width/2, SMAF[(0,8,16,24),2],  
                       width, yerr=SMAF[(1,9,17,25),2], color=color_list[3], label='CS-SMAF')
rects5 = axes[2].bar(ind + width/2 + width, np.mean(MAE_3[0:5, :], axis=0),  
                       width, yerr=np.std(MAE_3[0:5, :], axis=0), color=color_list[4], label='DeepAE')

axes[2].set_ylabel('MAE')
axes[2].set_xticks(ind)
axes[2].set_xticklabels(('10', '25', '50', '100'))
#axes[0,2].set_title('GSE45234')
axes[2].yaxis.set_major_formatter(formatter) 
axes[2].legend(loc=1)

plt.show()
fig.savefig('./Data/ceytomtry.pdf', dpi=200, bbox_inches='tight')