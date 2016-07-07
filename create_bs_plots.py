import os
import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    files = glob.glob('PFOC_3dQwarpYO_LV*bs_plot.mat')
    sns.set_style("white")

    for f in files:
        fname = os.path.splitext(f)[0]
        data_array = h5py.File(f,'r')
        bsmean = np.array(data_array.get('bs_data_mean'))

        plt.figure(figsize=(12, 9))
        plt.figure(facecolor="white")
        ax = plt.subplot(111)
        ax.spines["top"].set_visible(False)    
        ax.spines["right"].set_visible(False)    
        ax.get_xaxis().tick_bottom()    
        ax.get_yaxis().tick_left()
        ax.xaxis.set_tick_params(pad=15)
        ax.yaxis.set_tick_params(pad=15)
        ax.axhline(y=0, color="k")
        plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on") 

        plt.xlabel('TRs',fontsize=16)
        plt.ylabel('Mean Brain Score',fontsize=16)
        plt.xlim(-0.1,7.1)
        plt.ylim(round(np.amin(bsmean)-1), round(np.amax(bsmean)+1))
        plt.tick_params(labelsize=16)
        plt.plot(bsmean.T[2], color="green",  linewidth=2, linestyle="-", marker="D", label="Past")
        plt.plot(bsmean.T[3], color="red",  linewidth=2, linestyle="-", marker="D", label="Future")
        plt.plot(bsmean.T[4], color="blue",  linewidth=2, linestyle="-", marker="D", label="Other")
        plt.plot(bsmean.T[0], color="black",  linewidth=2, linestyle="-", marker="D", label="Control")
        plt.plot(bsmean.T[1], color="orange",  linewidth=2, linestyle="-", marker="D", label="Null")
        lgd = plt.legend(loc='center right',numpoints=1,fontsize=16,handleheight=2,frameon=False,bbox_to_anchor=(1.3,0.5))
        plt.savefig(fname+'.png',bbox_extra_artists=(lgd,),bbox_inches='tight')