
import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.yieldProducer_cfi")
process.yieldParams = process.FWLiteParams.clone()

process.yieldParams.selectionParams.wwmumu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu.mZ   = cms.double(-1.0) 
process.yieldParams.selectionParams.wwelmu.pMet = cms.double(20.0) 

from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff         import wwSamples
isData = False

looper = wwSamples

if isData:
    lumis = LumiList(filename = '/nfs/bluearc/group/jsons/Skim_Mu_And_El_byGiovanni/Final.json')
    process.yieldParams.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.yieldParams.lumisToProcess = lumis.getCMSSWString().split(',')


#loop over the looper
# for iname, iid, iscale, toPrint, muVeto in looper:
for x in looper:
    #make new pset
    setattr(
        process.yieldParams.inputParams,
        looper[x][0],
        cms.PSet(
            files = cms.vstring(),
            scale = cms.double(looper[x][1]),
            printEvents = cms.bool(isData),
        )
    )
    getattr(process.yieldParams.inputParams,looper[x][0]).files += looper[x][3]



