#!/usr/bin/env python

'''
# Script to extract and organize onset times
# from PLS datamat batch text files. 
# Can also write for AFNI GLM with TRs converted
# to timing in seconds.
'''

# first attempt!
# fname = "1511044_buttonPress_timing.txt"
# shift_time = -2
# onsetTimes = []
# with open(fname, 'r+') as f:
#     for line in f:
#         onsetTimes.append([float(i) for i in line.strip().split("\t")])
#     for run in range(count(onsetTimes)):
#         for onset in range(len(onsetTimes[run])):
#             onsetTimes[run][onset] = round(onsetTimes[run][onset] + shift_time,3)
#     f.seek(0)
#     f.writelines("\t".join(str(onset) for onset in run) + '\n' for run in onsetTimes)

import fileinput
from operator import itemgetter

def split_dmat_txt(fname):
    '''Split PLS datamat text files into lists,
    retaining run information.'''
    runsOnly = []
    fparts = []
    with open(fname, 'r+') as f:
        fparts.extend(f.read().split("data_files"))
    for part in fparts[1:len(fparts)+1]:
        runsOnly.append(part.strip().split("\n"))
        for sess in runsOnly:
            match  = ['run', 'event_onsets']
            sess[:] = [l for l in sess if any(m in l for m in match)]
    for sess in range(len(runsOnly)):
        for s in range(len(runsOnly[sess])):
            runsOnly[sess][s] = runsOnly[sess][s].split("\t")
    for sess in runsOnly:
        for cond in sess:
            cond[:] = [c for c in cond if "event_onsets" not in c]
    return runsOnly


def shift_onsets(shift_TR, runsOnly):
    '''Add specified number of TRs to every onset in list.'''
    for sess in runsOnly:
        for cond in sess[1:]:
            cond[:] = ([float(i) for i in cond])
    for sess in range(len(runsOnly)):
        for cond in range(1,len(runsOnly[sess])):
            for onset in range(len(runsOnly[sess][cond])):
                runsOnly[sess][cond][onset] = round(runsOnly[sess][cond][onset] + 
                                                    shift_TR,3)
    return runsOnly


def sort_and_write(TR_in_S, runsOnly, fname):
    '''Sort runs by chronological order, convert onsets 
    to seconds, and pull conditions to text files for AFNI.'''
    runsOnly.sort(key=itemgetter(0))
    for sess in range(len(runsOnly)):
        for cond in range(1,len(runsOnly[sess])):
            with open("%s_condition_%s.txt" % (fname.split('.')[0], cond), "a+") as newfile:
                newfile.write("\t".join(str(onset*TR_in_S) for onset in runsOnly[sess][cond])
                                 + '\n')
    return runsOnly


def write_back_to_dmat(runsOnly, fname):
    '''Write shifted onsets back to dmat text file.'''
    for sess in range(len(runsOnly)):
        for cond in range(1,len(runsOnly[sess])):
            runsOnly[sess][cond].insert(0,'event_onsets')
        runsOnly[sess][0].insert(0,'data_files')
        for cond in range(len(runsOnly[sess])):
            runsOnly[sess][cond] = '\t'.join([str(onset) for onset in runsOnly[sess][cond]])
    for sess in range(len(runsOnly)):
        for cond in range(1,len(runsOnly[sess])):
            runsOnly[sess][cond] = runsOnly[sess][cond].replace('.0','')
    fileloc = []
    with open(fname, 'r+') as f:
        text = f.read()
        for sess in runsOnly:
            fileloc.append(text.find(sess[0]))
        for sess in range(len(runsOnly)):
            f.seek(fileloc[sess])
            f.write('\n'.join(runsOnly[sess])+ '\n')

if __name__ == '__main__':
    import glob
    files_to_shift = glob.glob('PFOC1*.txt')

    # to shift existing dmats.
    # for fil in files_to_shift:
    #     runsOnly = split_dmat_txt(fil)
    #     runsOnly = shift_onsets(-2, runsOnly)
    #     write_back_to_dmat(runsOnly, fil)

    # to write existing dmats to text files for AFNI.
    for fil in files_to_shift:
        runsOnly = split_dmat_txt(fil)
        runsOnly = shift_onsets(0, runsOnly)
        runsOnly = sort_and_write(2, runsOnly, fil)

