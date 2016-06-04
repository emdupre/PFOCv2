# Note: If this is your first time running R from within 
# python, you'll need to execute the following commands:
# !conda install -c r rpy2
# from rpy2.robjects.packages import importr 
# utils = importr('utils')
# utils.install_packages('ggplot2', repos='http://cran.rstudio.com/')
# utils.install_packages('wesanderson', repos='http://cran.rstudio.com/')

import os
import glob
import numpy as np
import scipy.io as sio
import h5py
import pandas as pd
import rpy2.robjects as robj
import rpy2.robjects.pandas2ri
from rpy2.robjects.packages import importr

###############################
### Setting default parameters:
alpha = 0.05
GroupLevels = ['Young','Old']
###############################

def set_group_lvls(nRepeat,GroupLevels,nGroup):
    Group = []
    for i in range(0,nGroup):
        Group.extend([GroupLevels[i]]*nRepeat)
    Group = pd.DataFrame(Group,columns=['Group'])
    return Group


def event_related(data_array):
    Estimate = np.array(data_array.get('boot_result')['orig_usc'])
    UL = np.array(data_array.get('boot_result')['ulusc'])
    LL = np.array(data_array.get('boot_result')['llusc'])
    Significance = np.array(data_array.get('perm_result')['s_prob'])
    return Estimate, UL, LL, Significance


def extract_hdf5(data_array,ftype):
    if ftype == 'block':
        Estimate = np.array(data_array.get('result')['boot_result']['orig_usc'])
        UL = np.array(data_array.get('result')['boot_result']['ulusc'])
        LL = np.array(data_array.get('result')['boot_result']['llusc'])
        Significance = np.array(data_array.get('result')['perm_result']['sprob'])
        nGroup = np.array(data_array.get('result')['num_subj_lst']).shape[1]
        nRepeat = Estimate.shape[1]/nGroup
    elif ftype == 'event':
        Estimate, UL, LL, Significance = event_related(data_array)
        nGroup = np.array(data_array.get('subj_group')).shape[1]
        nRepeat = Estimate.shape[1]/nGroup

    mask=Significance[0]<alpha
    sigLV = np.where(mask)
    sigEstimate=Estimate[sigLV].T
    sigUL=UL[sigLV].T
    sigLL=LL[sigLV].T

    Group = set_group_lvls(nRepeat,GroupLevels,nGroup)
    cond_array = data_array['cond_name']
    Condition = []
    for i in range(0, cond_array.shape[0]):
        ascii_int = data_array[cond_array[i][0]]
        Condition.append(''.join(chr(i) for i in ascii_int[:]))
    Condition = pd.DataFrame(Condition, columns = ['Condition'])
    Condition = pd.concat([Condition] * nGroup, ignore_index=True)

    colnames = (['Estimate_LV'+str(i) for i in sigLV[0]+1] + 
    ['UL_LV'+str(i) for i in sigLV[0]+1] + 
    ['LL_LV'+str(i) for i in sigLV[0]+1])

    df = pd.DataFrame(np.hstack((sigEstimate,sigUL,sigLL)))
    df.columns = colnames
    df = pd.concat([df,Condition,Group],axis = 1)

    return df


def extract_unicode(data_array,ftype):
    if ftype == 'block':
        Estimate = data_array.get('result')['boot_result'][0,0]['orig_usc']
        UL = data_array.get('result')['boot_result'][0,0]['ulusc']
        LL = data_array.get('result')['boot_result'][0,0]['llusc']
        Significance = data_array.get('result')['perm_result'][0,0]['sprob']
        nGroup = data_array.get('result')['num_subj_lst'][0,0].shape[1]
        nRepeat = Estimate[0,0].shape[1]/nGroup
    elif ftype == 'event':
        Estimate, UL, LL, Significance = event_related(data_array)
        nGroup = data_array.get('subj_group').shape[1]
        nRepeat = Estimate[0,0].shape[1]/nGroup

    mask=Significance[0,0]<alpha
    sigLV = np.where(mask)[0]
    sigEstimate=pd.DataFrame(Estimate[0,0])[sigLV]
    sigUL=pd.DataFrame(UL[0,0])[sigLV]
    sigLL=pd.DataFrame(LL[0,0])[sigLV]

    Group = set_group_lvls(nRepeat,GroupLevels,nGroup)
    cond_array = data_array['cond_name']
    Condition = []
    for i in range(0, cond_array.shape[1]):
        Condition.extend(cond_array[0][i])
    Condition = pd.DataFrame(Condition, columns = ['Condition'])
    Condition = pd.concat([Condition] * nGroup, ignore_index=True)

    colnames = (['Estimate_LV'+str(i) for i in sigLV+1] + 
    ['UL_LV'+str(i) for i in sigLV+1] + 
    ['LL_LV'+str(i) for i in sigLV+1])

    df = pd.concat([sigEstimate,sigUL,sigLL],axis=1)
    df.columns = colnames
    df = pd.concat([df,Condition,Group],axis = 1)

    return df


def plot_w_ggplot2(f,df):
    plot_func = robj.r("""
        library(ggplot2)
        library(wesanderson)

        function(fname,pandasDF){

            nsigLVs = (ncol(pandasDF)-2)/3
            if (nlevels(pandasDF$Group)==1){
              for(i in 1:nsigLVs){
                ggsave(filename=paste(file_path_sans_ext(fname),colnames(pandasDF[i]),".png",sep=""),
                        plot=ggplot(pandasDF, aes(x=Condition, y=pandasDF[i])) +
                        geom_bar(width=.75,position=position_dodge(), stat="identity",
                                 size=.2, fill="#899DA4") +
                        geom_errorbar(aes(ymin=pandasDF[i+(2*nsigLVs)],
                                          ymax=pandasDF[i+nsigLVs]),
                                      width=.1,
                                      position=position_dodge(.75),
                                      colour="black") +
                        theme_minimal(base_size = 28, base_family = "Arial") +
                        theme(axis.text.y = element_blank()) +
                        theme(axis.title.y = element_blank()) +
                        theme(axis.title.x = element_text(margin = margin(t= 22))))
                }
            } else if (nlevels(pandasDF$Group)>1){
                for(i in 1:nsigLVs){
                    ggsave(filename=paste(file_path_sans_ext(fname),colnames(pandasDF[i]),".png",sep=""),
                            plot=ggplot(pandasDF, aes(x=Condition, y=pandasDF[i], fill=Group)) +
                            geom_bar(width=.75,position=position_dodge(), stat="identity",
                                     size=.2) +
                            geom_errorbar(aes(ymin=pandasDF[i+(2*nsigLVs)],
                                              ymax=pandasDF[i+nsigLVs]),
                                          width=.1,
                                          position=position_dodge(.75),
                                          colour="black") +
                            theme_minimal(base_size = 28, base_family = "Arial") +
                            theme(axis.text.y = element_blank()) +
                            theme(axis.title.y = element_blank()) +
                            theme(axis.title.x = element_text(margin = margin(t= 22))) +
                            scale_fill_manual(values=wes_palette("Royal1")))
                }
            }
        }
        """)
    robj.pandas2ri.activate()
    df_R = robj.conversion.py2ri(df)
    plot_func(f,df_R)


if __name__ == '__main__':
    files = glob.glob('*result.mat')

    for f in files:
        if f.find('_BfMRIresult.mat') >= 0: ftype = 'block'
        elif f.find('_fMRIresult.mat') >=0: ftype = 'event'
        else: print('Check file type, or give up all hope.')

        try:
            data_array = sio.loadmat(f)
            df = extract_unicode(data_array,ftype)
        except NotImplementedError:
            data_array = h5py.File(f,'r')
            df = extract_hdf5(data_array,ftype)

        plot_w_ggplot2(f,df)
