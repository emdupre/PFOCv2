import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def fill_corr_matrices(files):

    # Generate empty matrices to fill
    corr = np.zeros((len(files),24,24,5))
    zcorr = np.zeros((len(files),24,24,5))
    avgcorr = np.zeros((24,24,5))

    # Loop through subjects, generating correlation matrices for each condition
    subj = 0
    for f in files:
        data = np.loadtxt(f)
        x, y = data.shape
        seeds = y//8
        data = data.reshape(x,seeds,8)

        for c in range(data.shape[0]):
            corr[subj,:,:,c] = np.corrcoef(data[c])
            zcorr[subj,:,:,c] = np.arctanh(corr[subj,:,:,c].round(3))
        subj += 1

    avgcorr = zcorr.mean(axis = 0)
    return avgcorr


def plot_corr_matrices(avgcorr,group):
    # Average matrices and project correlation matrices for each condition 
    transcorr = np.tanh(avgcorr)
    condNames = ['Control', 'Null', 'Past', 'Future', 'Other']
    for cond in range(5):
        sns.set(style = "white") # set seaborn style
        f, ax = plt.subplots(figsize = (11, 9)) # generate empty figure
        sns.heatmap(transcorr[:,:,cond], cmap = "RdBu_r", center = 0, square = True,
                    xticklabels=['PFCdp','IPL','STS','MPFC','pHG','PCC','FEF', 'IPS', 'SPL7a', 'MT+', 'SPL7p','PrCv'],
                    yticklabels=['PFCdp','IPL','STS','MPFC','pHG','PCC','FEF', 'IPS', 'SPL7a', 'MT+', 'SPL7p','PrCv'])
        plt.axvline(6, color='k', lw = 1.5)
        plt.axhline(6, color='k', lw = 1.5)
        plt.savefig("corrMatrix_%s_%s.png" % (group, condNames[cond]))


if __name__ == '__main__':

    # Select files 
    Youngfiles = glob.glob('YeoCogNets_nbhd1_fMRI_grp1_subj*_voxeldata.txt')
    Olderfiles = glob.glob('YeoCogNets_nbhd1_fMRI_grp2_subj*_voxeldata.txt')

    # Calculate corr matrices
    YoungAvgCorr = fill_corr_matrices(Youngfiles)
    OldAvgCorr = fill_corr_matrices(Olderfiles)
    DiffAvgCorr = YoungAvgCorr - OldAvgCorr

    # Plot all matrices
    plot_corr_matrices(YoungAvgCorr, 'Young')
    plot_corr_matrices(OldAvgCorr, 'Older')
    plot_corr_matrices(DiffAvgCorr, 'YoungvOld')


# DN region labels: ['PFCdp','IPL','STS','MPFC','pHG','PCC']
# DAN region labels: ['FEF', 'IPS', 'SPL7a', 'MT+', 'SPL7p','PrCv']


