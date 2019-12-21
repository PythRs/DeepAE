import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data_path_5 = "./Data/mean_std_SMAF.csv"
SMAF = np.loadtxt(open(data_path_5), delimiter=",", skiprows=0)
# print(SMAF.shape)
PCC_SMAF = SMAF[:, (0, 3, 6, 9, 12, 15, 18, 21, 24)]
EM_SMAF = SMAF[:, (1, 4, 7, 10, 13, 16, 19, 22, 25)]
MAE_SMAF = SMAF[:, (2, 5, 8, 11, 14, 17, 20, 23, 26)]

data_path_3 = "./Data/DeepAE distinct measurements.csv"
data_measurement = np.loadtxt(open(data_path_3), delimiter=",", skiprows=0)

from matplotlib import rcParams

rcParams.update({'font.size': 14, 'font.family': 'STIXGeneral'})
PCC_3 = data_measurement[:, (9, 6, 3, 0)]
EM_3 = data_measurement[:, (10, 7, 4, 1)]
MAE_3 = data_measurement[:, (11, 8, 5, 2)]

ind = np.arange(4)  # the x locations for the groups
width = 0.15  # the width of the bars
color_list = plt.cm.Set3(np.linspace(0, 1, 12))

from matplotlib import rcParams

rcParams.update({'font.size': 14})

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 12))

rects1 = axes[0, 0].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 2],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 2], color=color_list[0], label='SVD')
rects2 = axes[0, 0].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 2],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 2], color=color_list[7], label='k-SVD')
rects3 = axes[0, 0].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 2],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 2], color=color_list[2], label='sNMF')
rects4 = axes[0, 0].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 2],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 2], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 0].bar(ind + width / 2 + width, np.mean(PCC_3[0:5, :], axis=0),
                        width, yerr=np.std(PCC_3[0:5, :], axis=0), color=color_list[4])
axes[0, 0].set_ylabel('PCC')
axes[0, 0].set_xticks(ind)
axes[0, 0].set_ylim(0.2, 1)
axes[0, 0].set_xticklabels(('10', '25', '50', '100'))
axes[0, 0].set_title('GSE60361')

rects1 = axes[0, 1].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 0],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 0], color=color_list[0], label='SVD')
rects2 = axes[0, 1].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 0],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 0], color=color_list[7], label='k-SVD')
rects3 = axes[0, 1].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 0],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 0], color=color_list[2], label='sNMF')
rects4 = axes[0, 1].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 0],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 0], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 1].bar(ind + width / 2 + width, np.mean(PCC_3[5:10, :], axis=0),
                        width, yerr=np.std(PCC_3[5:10, :], axis=0), color=color_list[4])

# axes[0,1].set_ylabel('PCC')
axes[0, 1].set_xticks(ind)
axes[0, 1].set_xticklabels(('10', '25', '50', '100'))
axes[0, 1].set_title('GSE102475')
axes[0, 1].set_ylim(0.5, 1)

rects1 = axes[0, 2].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 1],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 1], color=color_list[0], label='SVD')
rects2 = axes[0, 2].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 1],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 1], color=color_list[7], label='k-SVD')
rects3 = axes[0, 2].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 1],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 1], color=color_list[2], label='sNMF')
rects4 = axes[0, 2].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 1],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 1], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 2].bar(ind + width / 2 + width, np.mean(PCC_3[10:15, :], axis=0),
                        width, yerr=np.std(PCC_3[10:15, :], axis=0), color=color_list[4])

# axes[0,2].set_ylabel('PCC')
axes[0, 2].set_xticks(ind)
axes[0, 2].set_xticklabels(('10', '25', '50', '100'))
axes[0, 2].set_title('GSE52529')
axes[0, 2].set_ylim(0.5, 1)

rects1 = axes[1, 0].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 7],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 7], color=color_list[0], label='SVD')
rects2 = axes[1, 0].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 7],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 7], color=color_list[7], label='k-SVD')
rects3 = axes[1, 0].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 7],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 7], color=color_list[2], label='sNMF')
rects4 = axes[1, 0].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 7],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 7], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 0].bar(ind + width / 2 + width, np.mean(PCC_3[15:20, :], axis=0),
                        width, yerr=np.std(PCC_3[15:20, :], axis=0), color=color_list[4])
axes[1, 0].set_ylabel('PCC')
axes[1, 0].set_xticks(ind)
axes[1, 0].set_xticklabels(('10', '25', '50', '100'))
axes[1, 0].set_title('GSE69405')
axes[1, 0].set_ylim(0.5, 1)

rects1 = axes[1, 1].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 5],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 5], color=color_list[0], label='SVD')
rects2 = axes[1, 1].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 5],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 5], color=color_list[7], label='k-SVD')
rects3 = axes[1, 1].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 5],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 5], color=color_list[2], label='sNMF')
rects4 = axes[1, 1].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 5],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 5], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 1].bar(ind + width / 2 + width, np.mean(PCC_3[20:25, :], axis=0),
                        width, yerr=np.std(PCC_3[20:25, :], axis=0), color=color_list[4])
# axes[1,1].set_ylabel('PCC')
axes[1, 1].set_xticks(ind)
axes[1, 1].set_xticklabels(('10', '25', '50', '100'))
axes[1, 1].set_title('GSE65525')
axes[1, 1].set_ylim(0, 1)

rects1 = axes[1, 2].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 6],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 6], color=color_list[0], label='SVD')
rects2 = axes[1, 2].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 6],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 6], color=color_list[7], label='k-SVD')
rects3 = axes[1, 2].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 6],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 6], color=color_list[2], label='sNMF')
rects4 = axes[1, 2].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 6],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 6], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 2].bar(ind + width / 2 + width, np.mean(PCC_3[25:30, :], axis=0),
                        width, yerr=np.std(PCC_3[25:30, :], axis=0), color=color_list[4])
# axes[1,2].set_ylabel('PCC')
axes[1, 2].set_xticks(ind)
axes[1, 2].set_xticklabels(('10', '25', '50', '100'))
axes[1, 2].set_title('GSE78779')
axes[1, 2].set_ylim(0.5, 1)

rects1 = axes[2, 0].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 8],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 8], color=color_list[0], label='SVD')
rects2 = axes[2, 0].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 8],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 8], color=color_list[7], label='k-SVD')
rects3 = axes[2, 0].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 8],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 8], color=color_list[2], label='sNMF')
rects4 = axes[2, 0].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 8],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 8], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 0].bar(ind + width / 2 + width, np.mean(PCC_3[30:35, :], axis=0),
                        width, yerr=np.std(PCC_3[30:35, :], axis=0), color=color_list[4])
axes[2, 0].set_ylabel('PCC')
axes[2, 0].set_xticks(ind)
axes[2, 0].set_xticklabels(('10', '25', '50', '100'))
axes[2, 0].set_title('GSE84133')
axes[2, 0].set_ylim(0, 1)

rects1 = axes[2, 1].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 3],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 3], color=color_list[0], label='SVD')
rects2 = axes[2, 1].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 3],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 3], color=color_list[7], label='k-SVD')
rects3 = axes[2, 1].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 3],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 3], color=color_list[2], label='sNMF')
rects4 = axes[2, 1].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 3],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 3], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 1].bar(ind + width / 2 + width, np.mean(PCC_3[35:40, :], axis=0),
                        width, yerr=np.std(PCC_3[35:40, :], axis=0), color=color_list[4])
# axes[2,1].set_ylabel('PCC')
axes[2, 1].set_xticks(ind)
axes[2, 1].set_xticklabels(('10', '25', '50', '100'))
axes[2, 1].set_title('GSE62270')
axes[2, 1].set_xlabel('Measurements')
axes[2, 1].set_ylim(0.5, 1)

rects1 = axes[2, 2].bar(ind - width / 2 - 2 * width, PCC_SMAF[(2, 10, 18, 26), 4],
                        width, yerr=PCC_SMAF[(3, 11, 19, 27), 4], color=color_list[0], label='SVD')
rects2 = axes[2, 2].bar(ind - width / 2 - width, PCC_SMAF[(4, 12, 20, 28), 4],
                        width, yerr=PCC_SMAF[(5, 13, 21, 29), 4], color=color_list[7], label='k-SVD')
rects3 = axes[2, 2].bar(ind - width / 2, PCC_SMAF[(6, 14, 22, 30), 4],
                        width, yerr=PCC_SMAF[(7, 15, 23, 31), 4], color=color_list[2], label='sNMF')
rects4 = axes[2, 2].bar(ind + width / 2, PCC_SMAF[(0, 8, 16, 24), 4],
                        width, yerr=PCC_SMAF[(1, 9, 17, 25), 4], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 2].bar(ind + width / 2 + width, np.mean(PCC_3[40:45, :], axis=0),
                        width, yerr=np.std(PCC_3[40:45, :], axis=0), color=color_list[4], label='DeepAE')
# axes[2,2].set_ylabel('PCC')
axes[2, 2].set_xticks(ind)
axes[2, 2].set_xticklabels(('10', '25', '50', '100'))
axes[2, 2].set_title('GSE48968')
axes[2, 2].legend(loc=4)
axes[2, 2].set_ylim(0.5, 1)

plt.show()
fig.savefig('./Data/PCC.pdf', dpi=600, bbox_inches='tight')

ind = np.arange(4)  # the x locations for the groups
width = 0.15  # the width of the bars
color_list = plt.cm.Set3(np.linspace(0, 1, 12))

from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 1))

from matplotlib import rcParams

rcParams.update({'font.size': 14, 'font.family': 'STIXGeneral'})

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 12))

rects1 = axes[0, 0].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 2],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 2], color=color_list[0], label='SVD')
rects2 = axes[0, 0].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 2],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 2], color=color_list[7], label='k-SVD')
rects3 = axes[0, 0].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 2],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 2], color=color_list[2], label='sNMF')
rects4 = axes[0, 0].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 2],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 2], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 0].bar(ind + width / 2 + width, np.mean(EM_3[0:5, :], axis=0),
                        width, yerr=np.std(EM_3[0:5, :], axis=0), color=color_list[4])
axes[0, 0].set_ylabel('EM')
axes[0, 0].set_xticks(ind)
axes[0, 0].yaxis.set_major_formatter(formatter)
axes[0, 0].set_xticklabels(('10', '25', '50', '100'))
axes[0, 0].set_title('GSE60361')

rects1 = axes[0, 1].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 0],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 0], color=color_list[0], label='SVD')
rects2 = axes[0, 1].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 0],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 0], color=color_list[7], label='k-SVD')
rects3 = axes[0, 1].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 0],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 0], color=color_list[2], label='sNMF')
rects4 = axes[0, 1].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 0],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 0], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 1].bar(ind + width / 2 + width, np.mean(EM_3[5:10, :], axis=0),
                        width, yerr=np.std(EM_3[5:10, :], axis=0), color=color_list[4])

# axes[0,1].set_ylabel('PCC')
axes[0, 1].set_xticks(ind)
axes[0, 1].set_xticklabels(('10', '25', '50', '100'))
axes[0, 1].set_title('GSE102475')
axes[0, 1].yaxis.set_major_formatter(formatter)

rects1 = axes[0, 2].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 1],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 1], color=color_list[0], label='SVD')
rects2 = axes[0, 2].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 1],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 1], color=color_list[7], label='k-SVD')
rects3 = axes[0, 2].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 1],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 1], color=color_list[2], label='sNMF')
rects4 = axes[0, 2].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 1],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 1], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 2].bar(ind + width / 2 + width, np.mean(EM_3[10:15, :], axis=0),
                        width, yerr=np.std(EM_3[10:15, :], axis=0), color=color_list[4])
axes[0, 2].set_xticks(ind)
axes[0, 2].set_xticklabels(('10', '25', '50', '100'))
axes[0, 2].set_title('GSE52529')
axes[0, 2].yaxis.set_major_formatter(formatter)

rects1 = axes[1, 0].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 7],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 7], color=color_list[0], label='SVD')
rects2 = axes[1, 0].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 7],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 7], color=color_list[7], label='k-SVD')
rects3 = axes[1, 0].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 7],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 7], color=color_list[2], label='sNMF')
rects4 = axes[1, 0].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 7],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 7], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 0].bar(ind + width / 2 + width, np.mean(EM_3[15:20, :], axis=0),
                        width, yerr=np.std(EM_3[15:20, :], axis=0), color=color_list[4])
axes[1, 0].set_ylabel('EM')
axes[1, 0].set_xticks(ind)
axes[1, 0].set_xticklabels(('10', '25', '50', '100'))
axes[1, 0].set_title('GSE69405')
axes[1, 0].yaxis.set_major_formatter(formatter)

rects1 = axes[1, 1].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 5],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 5], color=color_list[0], label='SVD')
rects2 = axes[1, 1].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 5],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 5], color=color_list[7], label='k-SVD')
rects3 = axes[1, 1].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 5],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 5], color=color_list[2], label='sNMF')
rects4 = axes[1, 1].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 5],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 5], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 1].bar(ind + width / 2 + width, np.mean(EM_3[20:25, :], axis=0),
                        width, yerr=np.std(EM_3[20:25, :], axis=0), color=color_list[4])
# axes[1,1].set_ylabel('PCC')
axes[1, 1].set_xticks(ind)
axes[1, 1].set_xticklabels(('10', '25', '50', '100'))
axes[1, 1].set_title('GSE65525')
axes[1, 1].yaxis.set_major_formatter(formatter)

rects1 = axes[1, 2].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 6],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 6], color=color_list[0], label='SVD')
rects2 = axes[1, 2].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 6],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 6], color=color_list[7], label='k-SVD')
rects3 = axes[1, 2].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 6],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 6], color=color_list[2], label='sNMF')
rects4 = axes[1, 2].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 6],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 6], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 2].bar(ind + width / 2 + width, np.mean(EM_3[25:30, :], axis=0),
                        width, yerr=np.std(EM_3[25:30, :], axis=0), color=color_list[4])
# axes[1,2].set_ylabel('PCC')
axes[1, 2].set_xticks(ind)
axes[1, 2].set_xticklabels(('10', '25', '50', '100'))
axes[1, 2].set_title('GSE78779')
axes[1, 2].yaxis.set_major_formatter(formatter)

rects1 = axes[2, 0].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 8],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 8], color=color_list[0], label='SVD')
rects2 = axes[2, 0].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 8],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 8], color=color_list[7], label='k-SVD')
rects3 = axes[2, 0].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 8],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 8], color=color_list[2], label='sNMF')
rects4 = axes[2, 0].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 8],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 8], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 0].bar(ind + width / 2 + width, np.mean(EM_3[30:35, :], axis=0),
                        width, yerr=np.std(EM_3[30:35, :], axis=0), color=color_list[4])
axes[2, 0].set_ylabel('EM')
axes[2, 0].set_xticks(ind)
axes[2, 0].set_xticklabels(('10', '25', '50', '100'))
axes[2, 0].set_title('GSE84133')
axes[2, 0].yaxis.set_major_formatter(formatter)

rects1 = axes[2, 1].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 3],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 3], color=color_list[0], label='SVD')
rects2 = axes[2, 1].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 3],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 3], color=color_list[7], label='k-SVD')
rects3 = axes[2, 1].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 3],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 3], color=color_list[2], label='sNMF')
rects4 = axes[2, 1].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 3],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 3], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 1].bar(ind + width / 2 + width, np.mean(EM_3[35:40, :], axis=0),
                        width, yerr=np.std(EM_3[35:40, :], axis=0), color=color_list[4])
# axes[2,1].set_ylabel('PCC')
axes[2, 1].set_xticks(ind)
axes[2, 1].set_xticklabels(('10', '25', '50', '100'))
axes[2, 1].set_title('GSE62270')
axes[2, 1].set_xlabel('Measurements')
axes[2, 1].yaxis.set_major_formatter(formatter)

rects1 = axes[2, 2].bar(ind - width / 2 - 2 * width, EM_SMAF[(2, 10, 18, 26), 4],
                        width, yerr=EM_SMAF[(3, 11, 19, 27), 4], color=color_list[0], label='SVD')
rects2 = axes[2, 2].bar(ind - width / 2 - width, EM_SMAF[(4, 12, 20, 28), 4],
                        width, yerr=EM_SMAF[(5, 13, 21, 29), 4], color=color_list[7], label='k-SVD')
rects3 = axes[2, 2].bar(ind - width / 2, EM_SMAF[(6, 14, 22, 30), 4],
                        width, yerr=EM_SMAF[(7, 15, 23, 31), 4], color=color_list[2], label='sNMF')
rects4 = axes[2, 2].bar(ind + width / 2, EM_SMAF[(0, 8, 16, 24), 4],
                        width, yerr=EM_SMAF[(1, 9, 17, 25), 4], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 2].bar(ind + width / 2 + width, np.mean(EM_3[40:45, :], axis=0),
                        width, yerr=np.std(EM_3[40:45, :], axis=0), color=color_list[4], label='DeepAE')
# axes[2,2].set_ylabel('PCC')
axes[2, 2].set_xticks(ind)
axes[2, 2].set_xticklabels(('10', '25', '50', '100'))
axes[2, 2].set_title('GSE48968')
axes[2, 2].legend(loc=4)
axes[2, 2].yaxis.set_major_formatter(formatter)

plt.show()
fig.savefig('./Data/EM.pdf', dpi=600, bbox_inches='tight')

ind = np.arange(4)  # the x locations for the groups
width = 0.15  # the width of the bars
color_list = plt.cm.Set3(np.linspace(0, 1, 12))

from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 1))

from matplotlib import rcParams

rcParams.update({'font.size': 14, 'font.family': 'STIXGeneral'})

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 12))

rects1 = axes[0, 0].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 2],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 2], color=color_list[0], label='SVD')
rects2 = axes[0, 0].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 2],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 2], color=color_list[7], label='k-SVD')
rects3 = axes[0, 0].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 2],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 2], color=color_list[2], label='sNMF')
rects4 = axes[0, 0].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 2],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 2], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 0].bar(ind + width / 2 + width, np.mean(MAE_3[0:5, :], axis=0),
                        width, yerr=np.std(MAE_3[0:5, :], axis=0), color=color_list[4])
axes[0, 0].set_ylabel('MAE')
axes[0, 0].set_xticks(ind)
axes[0, 0].yaxis.set_major_formatter(formatter)
axes[0, 0].set_xticklabels(('10', '25', '50', '100'))
axes[0, 0].set_title('GSE60361')

rects1 = axes[0, 1].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 0],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 0], color=color_list[0], label='SVD')
rects2 = axes[0, 1].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 0],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 0], color=color_list[7], label='k-SVD')
rects3 = axes[0, 1].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 0],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 0], color=color_list[2], label='sNMF')
rects4 = axes[0, 1].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 0],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 0], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 1].bar(ind + width / 2 + width, np.mean(MAE_3[5:10, :], axis=0),
                        width, yerr=np.std(MAE_3[5:10, :], axis=0), color=color_list[4])

# axes[0,1].set_ylabel('PCC')
axes[0, 1].set_xticks(ind)
axes[0, 1].set_xticklabels(('10', '25', '50', '100'))
axes[0, 1].set_title('GSE102475')
# axes[0,1].set_ylim(0.5,1)
axes[0, 1].yaxis.set_major_formatter(formatter)

rects1 = axes[0, 2].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 1],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 1], color=color_list[0], label='SVD')
rects2 = axes[0, 2].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 1],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 1], color=color_list[7], label='k-SVD')
rects3 = axes[0, 2].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 1],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 1], color=color_list[2], label='sNMF')
rects4 = axes[0, 2].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 1],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 1], color=color_list[3], label='CS-SMAF')
rects5 = axes[0, 2].bar(ind + width / 2 + width, np.mean(MAE_3[10:15, :], axis=0),
                        width, yerr=np.std(MAE_3[10:15, :], axis=0), color=color_list[4])

axes[0, 2].set_xticks(ind)
axes[0, 2].set_xticklabels(('10', '25', '50', '100'))
axes[0, 2].set_title('GSE52529')
axes[0, 2].yaxis.set_major_formatter(formatter)

rects1 = axes[1, 0].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 7],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 7], color=color_list[0], label='SVD')
rects2 = axes[1, 0].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 7],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 7], color=color_list[7], label='k-SVD')
rects3 = axes[1, 0].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 7],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 7], color=color_list[2], label='sNMF')
rects4 = axes[1, 0].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 7],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 7], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 0].bar(ind + width / 2 + width, np.mean(MAE_3[15:20, :], axis=0),
                        width, yerr=np.std(MAE_3[15:20, :], axis=0), color=color_list[4])
axes[1, 0].set_ylabel('MAE')
axes[1, 0].set_xticks(ind)
axes[1, 0].set_xticklabels(('10', '25', '50', '100'))
axes[1, 0].set_title('GSE69405')
axes[1, 0].yaxis.set_major_formatter(formatter)

rects1 = axes[1, 1].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 5],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 5], color=color_list[0], label='SVD')
rects2 = axes[1, 1].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 5],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 5], color=color_list[7], label='k-SVD')
rects3 = axes[1, 1].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 5],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 5], color=color_list[2], label='sNMF')
rects4 = axes[1, 1].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 5],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 5], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 1].bar(ind + width / 2 + width, np.mean(MAE_3[20:25, :], axis=0),
                        width, yerr=np.std(MAE_3[20:25, :], axis=0), color=color_list[4])
# axes[1,1].set_ylabel('PCC')
axes[1, 1].set_xticks(ind)
axes[1, 1].set_xticklabels(('10', '25', '50', '100'))
axes[1, 1].set_title('GSE65525')
axes[1, 1].yaxis.set_major_formatter(formatter)

rects1 = axes[1, 2].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 6],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 6], color=color_list[0], label='SVD')
rects2 = axes[1, 2].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 6],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 6], color=color_list[7], label='k-SVD')
rects3 = axes[1, 2].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 6],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 6], color=color_list[2], label='sNMF')
rects4 = axes[1, 2].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 6],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 6], color=color_list[3], label='CS-SMAF')
rects5 = axes[1, 2].bar(ind + width / 2 + width, np.mean(MAE_3[25:30, :], axis=0),
                        width, yerr=np.std(MAE_3[25:30, :], axis=0), color=color_list[4])
# axes[1,2].set_ylabel('PCC')
axes[1, 2].set_xticks(ind)
axes[1, 2].set_xticklabels(('10', '25', '50', '100'))
axes[1, 2].set_title('GSE78779')
axes[1, 2].yaxis.set_major_formatter(formatter)

rects1 = axes[2, 0].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 8],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 8], color=color_list[0], label='SVD')
rects2 = axes[2, 0].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 8],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 8], color=color_list[7], label='k-SVD')
rects3 = axes[2, 0].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 8],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 8], color=color_list[2], label='sNMF')
rects4 = axes[2, 0].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 8],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 8], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 0].bar(ind + width / 2 + width, np.mean(MAE_3[30:35, :], axis=0),
                        width, yerr=np.std(MAE_3[30:35, :], axis=0), color=color_list[4])
axes[2, 0].set_ylabel('MAE')
axes[2, 0].set_xticks(ind)
axes[2, 0].set_xticklabels(('10', '25', '50', '100'))
axes[2, 0].set_title('GSE84133')
axes[2, 0].yaxis.set_major_formatter(formatter)

rects1 = axes[2, 1].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 3],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 3], color=color_list[0], label='SVD')
rects2 = axes[2, 1].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 3],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 3], color=color_list[7], label='k-SVD')
rects3 = axes[2, 1].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 3],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 3], color=color_list[2], label='sNMF')
rects4 = axes[2, 1].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 3],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 3], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 1].bar(ind + width / 2 + width, np.mean(MAE_3[35:40, :], axis=0),
                        width, yerr=np.std(MAE_3[35:40, :], axis=0), color=color_list[4])
# axes[2,1].set_ylabel('PCC')
axes[2, 1].set_xticks(ind)
axes[2, 1].set_xticklabels(('10', '25', '50', '100'))
axes[2, 1].set_title('GSE62270')
axes[2, 1].set_xlabel('Measurements')
axes[2, 1].yaxis.set_major_formatter(formatter)

rects1 = axes[2, 2].bar(ind - width / 2 - 2 * width, MAE_SMAF[(2, 10, 18, 26), 4],
                        width, yerr=MAE_SMAF[(3, 11, 19, 27), 4], color=color_list[0], label='SVD')
rects2 = axes[2, 2].bar(ind - width / 2 - width, MAE_SMAF[(4, 12, 20, 28), 4],
                        width, yerr=MAE_SMAF[(5, 13, 21, 29), 4], color=color_list[7], label='k-SVD')
rects3 = axes[2, 2].bar(ind - width / 2, MAE_SMAF[(6, 14, 22, 30), 4],
                        width, yerr=MAE_SMAF[(7, 15, 23, 31), 4], color=color_list[2], label='sNMF')
rects4 = axes[2, 2].bar(ind + width / 2, MAE_SMAF[(0, 8, 16, 24), 4],
                        width, yerr=MAE_SMAF[(1, 9, 17, 25), 4], color=color_list[3], label='CS-SMAF')
rects5 = axes[2, 2].bar(ind + width / 2 + width, np.mean(MAE_3[40:45, :], axis=0),
                        width, yerr=np.std(MAE_3[40:45, :], axis=0), color=color_list[4], label='DeepAE')
# axes[2,2].set_ylabel('PCC')
axes[2, 2].set_xticks(ind)
axes[2, 2].set_xticklabels(('10', '25', '50', '100'))
axes[2, 2].set_title('GSE48968')
axes[2, 2].legend(loc=4)
axes[2, 2].yaxis.set_major_formatter(formatter)

plt.show()
fig.savefig('./Data/MAE.pdf', dpi=600, bbox_inches='tight')

m10 = np.mean(PCC_SMAF[0:8, :], axis=1)
m25 = np.mean(PCC_SMAF[8:16, :], axis=1)
m50 = np.mean(PCC_SMAF[16:24, :], axis=1)
m100 = np.mean(PCC_SMAF[24:32, :], axis=1)
from matplotlib import rcParams

rcParams.update({'font.size': 14, 'font.family': 'STIXGeneral'})

fig, axes = plt.subplots(figsize=(5, 5))
x = (1, 2, 3, 4)

y = (m10[0], m25[0], m50[0], m100[0])
yerrs = (m10[1], m25[1], m50[1], m100[1])
plt.errorbar(x, y, yerr=yerrs, label='SVD')

y = (m10[2], m25[2], m50[2], m100[2])
yerrs = (m10[3], m25[3], m50[3], m100[3])
plt.errorbar(x, y, yerr=yerrs, label='k-SVD')

y = (m10[4], m25[4], m50[4], m100[4])
yerrs = m10[5], m25[5], m50[5], m100[5]
plt.errorbar(x, y, yerr=yerrs, label='sNMF')

y = (m10[6], m25[6], m50[6], m100[6])
yerrs = (m10[7], m25[7], m50[7], m100[7])
plt.errorbar(x, y, yerr=yerrs, label='CS-SMAF')

y = (np.mean(PCC_3[:, 0]), np.mean(PCC_3[:, 1]), np.mean(PCC_3[:, 2]), np.mean(PCC_3[:, 3]))
yerrs = (np.std(PCC_3[:, 0]), np.std(PCC_3[:, 1]), np.std(PCC_3[:, 2]), np.std(PCC_3[:, 3]))
plt.errorbar(x, y, yerr=yerrs, label='DeepAE')

plt.legend(loc=4)
axes.set_xticklabels(('0', '10', '25', '50', '100'))
axes.set_ylabel('PCC')
axes.set_xlabel('Measurements')
plt.ylim(0.5, 1)
plt.show()
fig.savefig('./Data/meas.pdf', dpi=600, bbox_inches='tight')
