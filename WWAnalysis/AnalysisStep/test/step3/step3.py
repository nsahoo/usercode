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

if len(args) == 0: args = [ 'ggToH160toWWto2L2Nu', 101160, 0.003621062529384]
if len(args) != 3: raise RuntimeError, "step3.py dataset id json (for data) or step3.py dataset id scalefactor (for MC)"
## step3.py dataset id json   for data
## step3.py dataset id scalef for MC
dataset = ['MC','ggToH160toWWto2L2Nu']; id = 101160; 
scalef  = 0.003621062529384;
json    = None
mhiggs  = 0
args[0] = args[0].replace("2011Av2","").replace("2011A","")
if args[0] in [ 'SingleElectron', 'DoubleElectron', 'SingleMuon', 'DoubleMuon', 'MuEG']:
    dataset = [args[0]]; id = args[1]
    json    = args[2]
else:
    dataset = ['MC', args[0]]; id = args[1];
    scalef  = float(args[2])
    m = re.match("ggToH(\\d+)to.*", args[0])
    if m: mhiggs = int(m.group(1))

nverticesModule = cms.EDProducer("VertexMultiplicityCounter",
    probes = cms.InputTag("REPLACE_ME"),
    objects = cms.InputTag("offlinePrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
)
Template = cms.EDFilter("ProbeTreeProducer",
    cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && "+
                     "leptEtaCut(2.4,2.5) && "+
                     "ptMax > 20 && "+
                     "((abs(pdgIdByPt(1)) == 13 && ptMin > 10) || "+
                     " (abs(pdgIdByPt(1)) == 11 && ptMin > 15)) && "
                     "nExtraLep(10) == 0 "
                    +" && passesIP"
                   +(" && triggerMatchingCut('%s')" % dataset[0])
    ),
    variables = cms.PSet(
        hypo = cms.string("hypo()"),
        mll  = cms.string("mll()"),
        ptll = cms.string("pTll()"),
        yll  = cms.string("yll()"),
        pt1  = cms.string("ptMax"),
        pt2  = cms.string("ptMin"),
        met  = cms.string("pfMet"),
        mmet = cms.string("minMet"),
        pmet = cms.string("projPfMet"),
        #pmmet = cms.string("projMet('MINMET')"), ## no, min of proj and proj of min are not the same
        pmmet = cms.string("min(projPfMet,projChargedMet)"),
        dphill = cms.string("dPhill()"),
        drll   = cms.string("dRll()"),
        dphilljet  = cms.string("dPhillLeadingJet()"),
        dphillcjet = cms.string("dPhillLeadingJet(2.5)"),
        dphillmet  = cms.string("dPhillMet('PFMET')"),
        dphillmmet = cms.string("dPhillMet('MINMET')"),
        dphilmet = cms.string("dPhilMet('PFMET')"),
        mtw1 = cms.string("mTByPt(0,'PFMET')"),
        mtw2 = cms.string("mTByPt(1,'PFMET')"),
        mth  = cms.string("mTHiggs('PFMET')"),
        gammaMRStar = cms.string("gammaMRStar"),
        njet  = cms.string("nJets(30)"),
        ncjet = cms.string("nCentralJets(30,2.5)"),
        nbjet = cms.string("bTaggedJetsOver(30,2.1)"),
        bestiso  = cms.string("allIsoByIso(0)/ptByIso(0)"),
        worstiso = cms.string("allIsoByIso(1)/ptByIso(1)"),
        worstJetLepPt = cms.string("max(matchedJetPt(0, 0.5)/pt(0), matchedJetPt(1, 0.5)/pt(1))"),
        dataset = cms.string(str(id))
    ),
    flags = cms.PSet(
        sameflav = cms.string("hypo == 3 || hypo == 6"),
        zveto    = cms.string("abs(mll-91.1876)>15 || hypo == 4 || hypo == 5"),
        bveto    = cms.string("bTaggedJetsUnder(30,2.1) == 0 && nSoftMu(3) == 0"),
        bveto_ip = cms.string("bTaggedJetsUnder(30,2.1) == 0"),
        bveto_mu = cms.string("nSoftMu(3) == 0"),
    ),
)

if dataset[0] == "MC":
    from WWAnalysis.AnalysisStep.pileupReweighting_cfi import reWeightVector
    Template.eventWeight = cms.InputTag("mcWeight");
    process.mcWeight     = cms.EDProducer("CombinedWeightProducer",
        baseWeight = cms.double(scalef),
        puWeight   = cms.vdouble(*reWeightVector[:]),
        puLabel    = cms.InputTag("addPileupInfo"),
    )
    if mhiggs > 0:
        process.higgsPt = cms.EDProducer("HWWKFactorProducer",
            genParticlesTag = cms.InputTag("onlyHiggsGen"),
            inputFilename = cms.untracked.string("WWAnalysis/Misc/Scales/scalefactor.mh%d.dat" % mhiggs),
            ProcessID = cms.untracked.int32(10010),
            Debug =cms.untracked.bool(False)
        )
        process.mcWeight.ptWeight = cms.InputTag("higgsPt")
        process.weightPath = cms.Path(process.higgsPt + process.mcWeight)
    else:
        process.weightPath = cms.Path(process.mcWeight)
else:
    from FWCore.PythonUtilities.LumiList import LumiList
    import os    
    lumis = LumiList(filename = os.getenv('CMSSW_BASE')+'/src/WWAnalysis/Misc/Jsons/%s.json'%json)
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess = lumis.getCMSSWString().split(',')

    Template.addRunLumiInfo = cms.addRunLumiInfo = cms.bool(True)

for X in "elel", "mumu", "elmu", "muel":
    setattr(process, X+"Nvtx", nverticesModule.clone(probes = cms.InputTag("ww%sCONVLHT"% X)))
    tree = Template.clone(src = cms.InputTag("ww%sCONVLHT"% X));
    tree.variables.nvtx = cms.InputTag(X+"Nvtx")
    setattr(process, X+"Tree", tree)
    setattr(process, X+"Path", cms.Path(getattr(process,X+"Nvtx") + tree))# + getattr(process,X+"TreeNoPU")))

process.TFileService = cms.Service("TFileService",fileName = cms.string("tree.root"))
