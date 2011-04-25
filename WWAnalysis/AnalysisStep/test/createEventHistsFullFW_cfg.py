import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.yieldProducer_cfi import *
from WWAnalysis.AnalysisStep.cutPSets_cfi import *

process = cms.Process("Test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:RMMEFN'
#         'file:/data/mwlebour/WW_414_SKIM_V04_STEP2_V00/101160/ggToH160toWWto2L2Nu_1_1_jGU.root'
    )
)
# from glob import glob
# process.source.fileNames = [ 'file:%s'%x for x in glob('/nfs/bluearc/group/trees/hww/WW_414_SKIM_V00/RMMEFOLDER/*.root') ]

# import commands
# myDir = "/castor/cern.ch/user/m/mangano/hww/skims/WW_413_SKIM_V00/id101160.GGtoH160toWWto2L2Nu/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RMMEFN')
#     fileName = cms.string('hists.root')
)


process.eventHists = cms.EDAnalyzer("CreateEventHistsEDMWrapped",
    FWLiteParams.clone(),
    sampleName = cms.string("RMMENUM.RMMENAME")
#     sampleName = cms.string("id101160.ggToH160toWWto2L2Nu")
)
process.p = cms.Path(process.eventHists)


#0
swapOutPSetValues(process.eventHists.selectionParams.wwelel0,hReOpt160)
# process.eventHists.selectionParams.wwmumu0 = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmu0 = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuel0 = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmu0,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuel0,oppositeFlavor)

#LHL
# process.eventHists.selectionParams.wwmumuIDLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuIDLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelIDLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelIDLHL = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuIDLHL,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelIDLHL,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuISOLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuISOLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelISOLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelISOLHL = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuISOLHL,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelISOLHL,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuCONVLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuCONVLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelCONVLHL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelCONVLHL = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuCONVLHL,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelCONVLHL,oppositeFlavor)

process.eventHists.selectionParams.wwmumuIPLHL = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelmuIPLHL = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwmuelIPLHL = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelelIPLHL = process.eventHists.selectionParams.wwelel0.clone()
swapOutPSetValues(process.eventHists.selectionParams.wwelmuIPLHL,oppositeFlavor)
swapOutPSetValues(process.eventHists.selectionParams.wwmuelIPLHL,oppositeFlavor)

# LHT
# process.eventHists.selectionParams.wwmumuIDLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuIDLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelIDLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelIDLHT = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuIDLHT,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelIDLHT,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuISOLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuISOLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelISOLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelISOLHT = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuISOLHT,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelISOLHT,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuCONVLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuCONVLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelCONVLHT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelCONVLHT = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuCONVLHT,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelCONVLHT,oppositeFlavor)

process.eventHists.selectionParams.wwmumuIPLHT = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelmuIPLHT = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwmuelIPLHT = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelelIPLHT = process.eventHists.selectionParams.wwelel0.clone()
swapOutPSetValues(process.eventHists.selectionParams.wwelmuIPLHT,oppositeFlavor)
swapOutPSetValues(process.eventHists.selectionParams.wwmuelIPLHT,oppositeFlavor)

# CBL
# process.eventHists.selectionParams.wwmumuIDCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuIDCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelIDCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelIDCBL = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuIDCBL,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelIDCBL,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuISOCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuISOCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelISOCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelISOCBL = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuISOCBL,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelISOCBL,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuCONVCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuCONVCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelCONVCBL = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelCONVCBL = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuCONVCBL,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelCONVCBL,oppositeFlavor)

process.eventHists.selectionParams.wwmumuIPCBL = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelmuIPCBL = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwmuelIPCBL = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelelIPCBL = process.eventHists.selectionParams.wwelel0.clone()
swapOutPSetValues(process.eventHists.selectionParams.wwelmuIPCBL,oppositeFlavor)
swapOutPSetValues(process.eventHists.selectionParams.wwmuelIPCBL,oppositeFlavor)

# CBT
# process.eventHists.selectionParams.wwmumuIDCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuIDCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelIDCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelIDCBT = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuIDCBT,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelIDCBT,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuISOCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuISOCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelISOCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelISOCBT = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuISOCBT,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelISOCBT,oppositeFlavor)
# 
# process.eventHists.selectionParams.wwmumuCONVCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelmuCONVCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwmuelCONVCBT = process.eventHists.selectionParams.wwelel0.clone()
# process.eventHists.selectionParams.wwelelCONVCBT = process.eventHists.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.eventHists.selectionParams.wwelmuCONVCBT,oppositeFlavor)
# swapOutPSetValues(process.eventHists.selectionParams.wwmuelCONVCBT,oppositeFlavor)

process.eventHists.selectionParams.wwmumuIPCBT = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelmuIPCBT = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwmuelIPCBT = process.eventHists.selectionParams.wwelel0.clone()
process.eventHists.selectionParams.wwelelIPCBT = process.eventHists.selectionParams.wwelel0.clone()
swapOutPSetValues(process.eventHists.selectionParams.wwelmuIPCBT,oppositeFlavor)
swapOutPSetValues(process.eventHists.selectionParams.wwmuelIPCBT,oppositeFlavor)

delattr(process.eventHists.selectionParams,"wwelel0")

