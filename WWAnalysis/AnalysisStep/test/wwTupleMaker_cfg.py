import FWCore.ParameterSet.Config as cms

process = cms.Process("EWKTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:data/REPLACEME.root'),
    #fileNames = cms.untracked.vstring('file:test.root'),
)
process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.hwwee = cms.EDAnalyzer("SkimTupleProducer", 
    src = cms.untracked.InputTag("hwwee"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("nlep"), quantity = cms.untracked.string("nLep()")),
        cms.PSet( tag = cms.untracked.string("neel"), quantity = cms.untracked.string("nExtraEl()")),
        cms.PSet( tag = cms.untracked.string("nel"), quantity = cms.untracked.string("nEl()")),
        cms.PSet( tag = cms.untracked.string("pdgid1"), quantity = cms.untracked.string("pdgid(0)")),
        cms.PSet( tag = cms.untracked.string("pdgid2"), quantity = cms.untracked.string("pdgid(1)")),
        cms.PSet( tag = cms.untracked.string("pt1"), quantity = cms.untracked.string("pt(0)")),
        cms.PSet( tag = cms.untracked.string("pt2"), quantity = cms.untracked.string("pt(1)")),
        cms.PSet( tag = cms.untracked.string("eta1"), quantity = cms.untracked.string("eta(0)")),
        cms.PSet( tag = cms.untracked.string("eta2"), quantity = cms.untracked.string("eta(1)")),
        cms.PSet( tag = cms.untracked.string("phi1"), quantity = cms.untracked.string("phi(0)")),
        cms.PSet( tag = cms.untracked.string("phi2"), quantity = cms.untracked.string("phi(1)")),
        cms.PSet( tag = cms.untracked.string("q1"), quantity = cms.untracked.string("q(0)")),
        cms.PSet( tag = cms.untracked.string("q2"), quantity = cms.untracked.string("q(1)")),
        cms.PSet( tag = cms.untracked.string("nj20"), quantity = cms.untracked.string("nJets(20)")),
        cms.PSet( tag = cms.untracked.string("ncj20"), quantity = cms.untracked.string("nCentralJets(20)")),
        cms.PSet( tag = cms.untracked.string("nj30"), quantity = cms.untracked.string("nJets(30)")),
        cms.PSet( tag = cms.untracked.string("ncj30"), quantity = cms.untracked.string("nCentralJets(30)")),
        cms.PSet( tag = cms.untracked.string("pf"), quantity = cms.untracked.string("pfMet()")),
        cms.PSet( tag = cms.untracked.string("tc"), quantity = cms.untracked.string("tcMet()")),
        cms.PSet( tag = cms.untracked.string("dphilpf1"), quantity = cms.untracked.string("dPhilPfMet(0)")),
        cms.PSet( tag = cms.untracked.string("dphiltc1"), quantity = cms.untracked.string("dPhilTcMet(0)")),
        cms.PSet( tag = cms.untracked.string("dphilpf2"), quantity = cms.untracked.string("dPhilPfMet(1)")),
        cms.PSet( tag = cms.untracked.string("dphiltc2"), quantity = cms.untracked.string("dPhilTcMet(1)")),
        cms.PSet( tag = cms.untracked.string("ppf"), quantity = cms.untracked.string("projPfMet(0,1)")),
        cms.PSet( tag = cms.untracked.string("ptc"), quantity = cms.untracked.string("projTcMet(0,1)")),
        cms.PSet( tag = cms.untracked.string("ppf1"), quantity = cms.untracked.string("projPfMet(0)")),
        cms.PSet( tag = cms.untracked.string("ptc1"), quantity = cms.untracked.string("projTcMet(0)")),
        cms.PSet( tag = cms.untracked.string("ppf2"), quantity = cms.untracked.string("projPfMet(1)")),
        cms.PSet( tag = cms.untracked.string("ptc2"), quantity = cms.untracked.string("projTcMet(1)")),
        cms.PSet( tag = cms.untracked.string("mll"), quantity = cms.untracked.string("mll()")),
        cms.PSet( tag = cms.untracked.string("dphill"), quantity = cms.untracked.string("dPhill()")),
        cms.PSet( tag = cms.untracked.string("ptll"), quantity = cms.untracked.string("pTll()")),
        cms.PSet( tag = cms.untracked.string("drll"), quantity = cms.untracked.string("dRll()")),
        cms.PSet( tag = cms.untracked.string("pfmt1"), quantity = cms.untracked.string("pfMT(0)")),
        cms.PSet( tag = cms.untracked.string("pfmt2"), quantity = cms.untracked.string("pfMT(1)")),
        cms.PSet( tag = cms.untracked.string("tcmt1"), quantity = cms.untracked.string("tcMT(0)")),
        cms.PSet( tag = cms.untracked.string("tcmt2"), quantity = cms.untracked.string("tcMT(1)")),
        cms.PSet( tag = cms.untracked.string("tcmll"), quantity = cms.untracked.string("tcMTll()")),
        cms.PSet( tag = cms.untracked.string("pfmll"), quantity = cms.untracked.string("pfMTll()")),
    ),
    weight = cms.untracked.double(RMMEWEIGHT)
    #weight = cms.untracked.double(1)
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("data/REPLACEME.tuple.root"))
#process.TFileService = cms.Service("TFileService",fileName = cms.string("test.tuple.root"))

process.p = cms.Path(process.hwwee)

