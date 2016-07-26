import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def fill_corr_matrices(files):

    # Generate empty matrices to fill
    corr = np.zeros((len(files),12,12,5))
    zcorr = np.zeros((len(files),12,12,5))
    avgcorr = np.zeros((12,12,5))

    # Loop through subjects, generating correlation matrices for each condition
    subj = 0
    for f in files:
        df = pd.read_csv(f, header=None, delim_whitespace=True)
        seed1 = df[list(range(0,8))]
        seed2 = df[list(range(8,16))]
        seed3 = df[list(range(16,24))]
        seed4 = df[list(range(24,32))]
        seed5 = df[list(range(32,40))]
        seed6 = df[list(range(40,48))]
        seed7 = df[list(range(48,56))]
        seed8 = df[list(range(56,64))]
        seed9 = df[list(range(64,72))]
        seed10 = df[list(range(72,80))]
        seed11 = df[list(range(80,88))]
        seed12 = df[list(range(88,96))]

        seeds = [seed1,seed2,seed3,seed4,seed5,seed6,seed7,seed8,seed9,seed10,seed11,seed12]

        for s in seeds:
            s.columns = ['TR0','TR1','TR2','TR3','TR4','TR5','TR6','TR7']

        for cond in range(5): # loop through the five conditions
            for s in range(len(seeds)): # loop through each seed
                for c in range(len(seeds)): # loop through each seed again
                    corr[subj,s,c,cond] = seeds[s].corrwith(seeds[c], axis = 1)[cond] # correlate and assign
            zcorr[subj,:,:,cond] = np.arctanh(corr[subj,:,:,cond].round(3)) # pull each condition and r-to-z transform
        subj += 1

    avgzcorr = zcorr.mean(axis = 0) # average but keep z-transformed in case comparing groups
    return avgzcorr


def plot_corr_matrices(avgzcorr,group):
    # Convert matrices back to r-values and project
    avgrcorr = np.tanh(avgzcorr) # transform to r-values for interpretation
    condNames = ['Control', 'Null', 'Past', 'Future', 'Other']
    for cond in range(5):
        sns.set(style = "white") # set seaborn style
        f, ax = plt.subplots(figsize = (11, 9)) # generate empty figure
        sns.heatmap(avgrcorr[:,:,cond], cmap = "RdBu_r", center = 0, square = True,
                    xticklabels=['PFCdp','IPL','STS','MPFC','pHG','PCC','FEF', 'IPS', 'SPL7a', 'MT+', 'SPL7p','PrCv'],
                    yticklabels=['PFCdp','IPL','STS','MPFC','pHG','PCC','FEF', 'IPS', 'SPL7a', 'MT+', 'SPL7p','PrCv'])
        # DN region labels: ['PFCdp','IPL','STS','MPFC','pHG','PCC']
        # DAN region labels: ['FEF', 'IPS', 'SPL7a', 'MT+', 'SPL7p','PrCv']
        plt.axvline(6, color='k', lw = 1.5) # vertical line to delineate DN and DAN
        plt.axhline(6, color='k', lw = 1.5) # hortizontal line to delineate DN and DAN
        plt.savefig("%s_%s_corrMatrix.png" % (group, condNames[cond]))


if __name__ == '__main__':

    # Select files 
    Youngfiles = glob.glob('YeoNet_fmri_grp1_subj*_voxeldata.txt')
    Olderfiles = glob.glob('YeoNet_fmri_grp2_subj*_voxeldata.txt')

    # Calculate corr matrices
    YoungAvgCorr = fill_corr_matrices(Youngfiles)
    OldAvgCorr = fill_corr_matrices(Olderfiles)
    DiffAvgCorr = YoungAvgCorr - OldAvgCorr

    # Plot all matrices
    plot_corr_matrices(YoungAvgzCorr, 'Young')
    plot_corr_matrices(OldAvgzCorr, 'Older')
    plot_corr_matrices(DiffAvgzCorr, 'Diff')




