import FWCore.ParameterSet.Config as cms
from glob import glob

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

sep172010BSamples = {
      #id    #name      #scale to 1pb^-1 #is data  #skim event files
    'id61': ['MuSep17'        , 1, True, [] ],
    'id62': ['EGSep17'        , 1, True, [] ],
    'id63': ['Mu2010B'        , 1, True, [] ],
    'id64': ['Electron2010B'  , 1, True, [] ],
}


for key in sep172010BSamples:
    sep172010BSamples[key][3] += [ '%s'%x for x in glob(prefix+key+'/*.root') ]























