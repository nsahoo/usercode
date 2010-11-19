import FWCore.ParameterSet.Config as cms
from FWCore.PythonUtilities.LumiList import LumiList
import math
from glob import glob

process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.yieldProducer_cfi")
process.yieldParams = process.FWLiteParams.clone()

process.yieldParams.selectionParams.wwmumu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu.mZ   = cms.double(-1.0) 
process.yieldParams.selectionParams.wwelmu.pMet = cms.double(20.0) 

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

looper = []
if False:
    #Loop over MC names, ids and scales
    looper = [ 
 #histo  #folder  #scale  #print events #veto mu triggered
        ['wz'    , 'id01', 2.77E-04, False, False],
        ['zz'    , 'id02', 9.39E-05, False, False],
        ['dyee'  , 'id03', 2.79E-02, False, False],
        ['dymm'  , 'id04', 2.79E-02, False, False],
        ['dytt'  , 'id05', 2.79E-02, False, False],
        ['wm'    , 'id06', 7.58E-02, False, False],
        ['we'    , 'id07', 9.10E-02, False, False],
        ['wt'    , 'id08', 8.74E-02, False, False],
        ['ttbar' , 'id09', 5.55E-04, False, False],
        ['ts'    , 'id10', 1.01E-04, False, False],
        ['tt'    , 'id11', 1.45E-03, False, False],
        ['tw'    , 'id12', 2.32E-04, False, False],
        ['ww'    , 'id13', 1.37E-03, False, False],
        ['wwNoPU', 'id16', 1.37E-03, False, False],
    ]
else:
    lumis = LumiList(filename = '/nfs/bluearc/group/jsons/Skim_Mu_And_El_byGiovanni/Final.json')
    process.yieldParams.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.yieldParams.lumisToProcess = lumis.getCMSSWString().split(',')

    #Loop over Data names, ids and scales
    looper = [
        ['muSep17', 'id61', 1.0, True, False],
        ['egSep17', 'id62', 1.0, True, True],
        ['mu2010B', 'id63', 1.0, True, False],
        ['eg2010B', 'id64', 1.0, True, True],
    ]
    #set the residual file to non-null to apply in data
    process.yieldParams.resFile = 'WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2L3Residual.txt'

    #add another set of plots to make with all data
#     process.yieldParams.inputParams.data = cms.PSet(
#         files = cms.vstring(  ),
#         scale = cms.double(1.),
#         printEvents = cms.bool(True)
#     )
#     process.yieldParams.inputParams.data.files += [ '%s'%x for x in glob(prefix+'id61/*.root') ]
#     process.yieldParams.inputParams.data.files += [ '%s'%x for x in glob(prefix+'id62/*.root') ]
#     process.yieldParams.inputParams.data.files += [ '%s'%x for x in glob(prefix+'id63/*.root') ]
#     process.yieldParams.inputParams.data.files += [ '%s'%x for x in glob(prefix+'id64/*.root') ]


#loop over the looper
for iname, iid, iscale, toPrint, muVeto in looper:
    #make new pset
    setattr(
        process.yieldParams.inputParams,
        iname,
        cms.PSet(
            files = cms.vstring(),
            scale = cms.double(iscale),
            printEvents = cms.bool(toPrint),
            vetoMuTriggered = cms.bool(muVeto)
        )
    )
    #glob the appropriate files to it
    getattr(process.yieldParams.inputParams,iname).files += [ '%s'%x for x in glob(prefix+iid+'/*.root') ]



