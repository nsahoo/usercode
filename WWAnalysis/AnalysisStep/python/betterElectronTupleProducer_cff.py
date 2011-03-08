import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.pileupReweighting_cfi import *
betterElectronTupleProducer = cms.EDAnalyzer("BetterElectronTupleProducer", 
    src = cms.untracked.InputTag("boostedElectrons"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("eid80"),              quantity = cms.untracked.string("electronID('eidVBTFRel80')")),
        cms.PSet( tag = cms.untracked.string("eid95"),              quantity = cms.untracked.string("electronID('eidVBTFRel95')")),
        cms.PSet( tag = cms.untracked.string("liklihoodID"),        quantity = cms.untracked.string("electronID('liklihoodID')")),
        cms.PSet( tag = cms.untracked.string("pt"),                 quantity = cms.untracked.string("pt")),
        cms.PSet( tag = cms.untracked.string("eta"),                quantity = cms.untracked.string("eta")),
        cms.PSet( tag = cms.untracked.string("phi"),                quantity = cms.untracked.string("phi")),
        cms.PSet( tag = cms.untracked.string("eb"),                 quantity = cms.untracked.string("isEB")),
        cms.PSet( tag = cms.untracked.string("ee"),                 quantity = cms.untracked.string("isEE")),
        cms.PSet( tag = cms.untracked.string("gap"),                quantity = cms.untracked.string("isEBEEGap")),
        cms.PSet( tag = cms.untracked.string("dxyPV"),              quantity = cms.untracked.string("userFloat('dxyPV')")),
        cms.PSet( tag = cms.untracked.string("dzPV"),               quantity = cms.untracked.string("userFloat('dzPV')")),
        cms.PSet( tag = cms.untracked.string("dcot"),               quantity = cms.untracked.string("userFloat('convValueMapProd:dcot')")),
        cms.PSet( tag = cms.untracked.string("dist"),               quantity = cms.untracked.string("userFloat('convValueMapProd:dist')")),
        cms.PSet( tag = cms.untracked.string("beta0003"),           quantity = cms.untracked.string("userFloat('beta0003')")),
        cms.PSet( tag = cms.untracked.string("beta0004"),           quantity = cms.untracked.string("userFloat('beta0004')")),
        cms.PSet( tag = cms.untracked.string("beta0703"),           quantity = cms.untracked.string("userFloat('beta0703')")),
        cms.PSet( tag = cms.untracked.string("beta0704"),           quantity = cms.untracked.string("userFloat('beta0704')")),
        cms.PSet( tag = cms.untracked.string("rho"),                quantity = cms.untracked.string("userFloat('rho')")),
        cms.PSet( tag = cms.untracked.string("sigma"),              quantity = cms.untracked.string("userFloat('sigma')")),
        cms.PSet( tag = cms.untracked.string("mom"),                quantity = cms.untracked.string("?genParticleRef().isNonnull()?genParticleRef().get().mother(0).mother(0).pdgId():9999")),
        cms.PSet( tag = cms.untracked.string("tk03"),               quantity = cms.untracked.string("dr03TkSumPt")),
        cms.PSet( tag = cms.untracked.string("ecal03"),             quantity = cms.untracked.string("dr03EcalRecHitSumEt")),
        cms.PSet( tag = cms.untracked.string("hcal03"),             quantity = cms.untracked.string("dr03HcalTowerSumEt")),
        cms.PSet( tag = cms.untracked.string("tk04"),               quantity = cms.untracked.string("dr04TkSumPt")),
        cms.PSet( tag = cms.untracked.string("ecal04"),             quantity = cms.untracked.string("dr04EcalRecHitSumEt")),
        cms.PSet( tag = cms.untracked.string("hcal04"),             quantity = cms.untracked.string("dr04HcalTowerSumEt")),
        cms.PSet( tag = cms.untracked.string("tkDefault"),          quantity = cms.untracked.string("userFloat('tkDefault')")),
        cms.PSet( tag = cms.untracked.string("ecalDefault"),        quantity = cms.untracked.string("userFloat('ecalDefault')")),
        cms.PSet( tag = cms.untracked.string("hcalDefault"),        quantity = cms.untracked.string("userFloat('hcalDefault')")),
    ), 
#     weight = cms.untracked.double(1),
    vtxWeights = cms.untracked.PSet(
        vtx05 = cms.untracked.vdouble( vtx05[:] ),
        vtx10 = cms.untracked.vdouble( vtx10[:] ),
        vtx15 = cms.untracked.vdouble( vtx15[:] ),
        vtx2010 = cms.untracked.vdouble( vtx2010[:] ),
    ),
    vtxLabel = cms.untracked.InputTag('offlinePrimaryVertices'),
    puLabel = cms.untracked.InputTag('addPileupInfo'),
    candCounts = cms.untracked.PSet(
#         nMu = cms.untracked.InputTag("goodMuons"),
    )

)  

