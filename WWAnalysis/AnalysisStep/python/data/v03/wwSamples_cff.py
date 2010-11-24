import FWCore.ParameterSet.Config as cms
from glob import glob

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

wwSamples = {
      #id    #name      #scale to 1pb^-1 #is data  #skim event files
    'id01': ['wz'    ,  8.29E-06,    False,     []],
    'id02': ['zz'    ,  2.81E-06,    False,     []],
    'id03': ['dyee'  ,  8.33E-04,    False,     []],
    'id04': ['dymm'  ,  8.33E-04,    False,     []],
    'id05': ['dytt'  ,  8.35E-04,    False,     []],
    'id06': ['wm'    ,  2.27E-03,    False,     []],
    'id07': ['we'    ,  2.72E-03,    False,     []],
    'id08': ['wt'    ,  2.62E-03,    False,     []],
    'id09': ['ttbar' ,  1.66E-05,    False,     []],
    'id10': ['ts'    ,  3.01E-06,    False,     []],
    'id11': ['tt'    ,  4.32E-05,    False,     []],
    'id12': ['tw'    ,  6.94E-06,    False,     []],
    'id13': ['ww'    ,  4.10E-05,    False,     []],
    'id16': ['wwNoPU',  4.10E-05,    False,     []],
}


#loop over the looper
for key in wwSamples:
    wwSamples[key][3] += [ '%s'%x for x in glob(prefix+key+'/*.root') ]

