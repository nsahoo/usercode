import FWCore.ParameterSet.Config as cms
import re
import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]

process = cms.Process("STEP3")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'file:/nfs/bluearc/group/trees/hww/R414_S1_V06_S2_V04_S3_V00/101160.ggToH160toWWto2L2Nu/ggToH160toWWto2L2Nu_1_2_yUh.root',
    ),
)
process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("WWAnalysis.AnalysisStep.step3_cff")

if len(args) == 0: args = [ 'ggToH160toWWto2L2Nu', 101160, 0.003621062529384]
if len(args) != 3: raise RuntimeError, "step3.py dataset id json (for data) or step3.py dataset id scalefactor (for MC)"
## step3.py dataset id json   for data
## step3.py dataset id scalef for MC
dataset = ['MC','ggToH160toWWto2L2Nu']; id = 101160; 
scalef  = 0.003621062529384;
json    = None
mhiggs  = 0
from WWAnalysis.AnalysisStep.fourthScaleFactors_cff import *
fourthGenSF = 1
# from WWAnalysis.AnalysisStep.scaleFactors_cff import *
# if args[1] in dataSamples or args[1] in data42xSamples:
args[0] = args[0].replace("2011Av2","").replace("2011AMay10","").replace("2011Av4","").replace("2011A","")
if args[0] in [ 'SingleElectron', 'DoubleElectron', 'SingleMuon', 'DoubleMuon', 'MuEG']:
    dataset = [args[0]]; id = args[1]
    json    = args[2]
    scalef = 1
else:
    dataset = ['MC', args[0]]; id = args[1];
    scalef  = float(args[2])
    m = re.match("ggToH(\\d+)to.*", args[0])
    if m: mhiggs = int(m.group(1))
    if m: fourthGenSF = fourthGenScales[int(m.group(1))]
process.step3Tree.cut = process.step3Tree.cut.value().replace("DATASET", dataset[0])
process.step3Tree.variables.trigmatch = process.step3Tree.variables.trigmatch.value().replace("DATASET",dataset[0])
process.step3Tree.variables.trigguil  = process.step3Tree.variables.trigguil.value().replace("DATASET",dataset[0])
process.step3Tree.variables.trigbits  = process.step3Tree.variables.trigbits.value().replace("DATASET",dataset[0])
process.step3Tree.variables.dataset = str(id)

if dataset[0] == "MC":
#     process.step3Tree.eventWeight = cms.InputTag("mcWeight");
#     process.mcWeight.baseWeight = scalef
    process.step3Tree.variables.lumiWeight = "%.12f" % scalef
    if mhiggs > 0:
#         process.higgsPt.inputFilename = "WWAnalysis/Misc/Scales/scalefactor.mh%d.dat" % mhiggs
        process.higgsPt.inputFilename = "HiggsAnalysis/HiggsToWW2Leptons/data/kfactors_Std/kfactors_mh%(mass)d_ren%(mass)d_fac%(mass)d.dat" % {"mass":mhiggs}
        process.step3Tree.variables.fourWeight = "%.12f" % fourthGenSF
    else:
        process.step3Tree.variables.fourWeight = "1"
else:
    from FWCore.PythonUtilities.LumiList import LumiList
    import os    
    lumis = LumiList(filename = os.getenv('CMSSW_BASE')+'/src/WWAnalysis/Misc/Jsons/%s.json'%json)
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess = lumis.getCMSSWString().split(',')
    process.step3Tree.variables.lumiWeight = "1"
    process.step3Tree.variables.fourWeight = "1"

for X in "elel", "mumu", "elmu", "muel":
    tree = process.step3Tree.clone(src = cms.InputTag("ww%sIPMerge"% X));
    setattr(process, X+"PtWeight", process.ptWeight.clone(src = cms.InputTag("ww%sIPMerge"% X)))
    tree.variables.ptWeight = cms.InputTag(X+"PtWeight")
    setattr(process, X+"PuWeight", process.puWeight.clone(src = cms.InputTag("ww%sIPMerge"% X)))
    tree.variables.puWeight = cms.InputTag(X+"PuWeight")
    setattr(process, X+"Nvtx", process.nverticesModule.clone(probes = cms.InputTag("ww%sIPMerge"% X)))
    tree.variables.nvtx = cms.InputTag(X+"Nvtx")
    setattr(process, X+"Tree", tree)
    setattr(process, X+"Path", cms.Path(
        getattr(process,X+"Nvtx") + 
        getattr(process,X+"PtWeight") + 
        getattr(process,X+"PuWeight") + 
        tree
    ))

process.TFileService = cms.Service("TFileService",fileName = cms.string("tree.root"))
