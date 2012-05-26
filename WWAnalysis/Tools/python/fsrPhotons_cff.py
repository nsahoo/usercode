import FWCore.ParameterSet.Config as cms

## PF Photons
fsrPhotons = cms.EDProducer("FSRPhotonProducer",
    srcCands = cms.InputTag("particleFlow"),
    ptThresh = cms.double(1.0), ## will tighten to 2 at analysis level
    extractMuonFSR = cms.bool(True),
)

## vetos as per muons 
fsrPhotonPFIsoChHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("fsrPhotons"),
    pfLabel     = cms.InputTag("pfNoPileUp"), 
    pfSelection = cms.string("charge != 0 && abs(pdgId) == 211"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0001), # self-veto 0.0001
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
fsrPhotonPFIsoNHad04 = fsrPhotonPFIsoChHad04.clone(
    pfSelection = "charge == 0 && abs(pdgId) == 130 && pt > 0.5", # neutral hadrons
    deltaRself  = 0.01, # larger veto cone for neutrals
)
fsrPhotonPFIsoPhoton04 = fsrPhotonPFIsoChHad04.clone(
    pfSelection = "charge == 0 && abs(pdgId) == 22 && pt > 0.5", # photons
    deltaRself  = 0.01, # larger veto cone for neutrals
)
# for deltaBeta corrections
fsrPhotonPFIsoChHadPU04 = fsrPhotonPFIsoChHad04.clone(
    pfLabel = 'pfPileUp'
)
# for rho corrections
import WWAnalysis.Tools.rhoValueMapProducer_cfi
rhoFsrPhoton = WWAnalysis.Tools.rhoValueMapProducer_cfi.rhoValueMapProducer.clone(
    leptonTag = "fsrPhotons",
    rhoTag = cms.untracked.InputTag("kt6PFJetsForIso","rho"),
)


import PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi 
boostedFsrPhotons = PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi.patPFParticles.clone(
    pfCandidateSource = 'fsrPhotons',
)
boostedFsrPhotons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("fsrPhotonPFIsoChHad04"),
    cms.InputTag("fsrPhotonPFIsoNHad04"),
    cms.InputTag("fsrPhotonPFIsoPhoton04"),
    cms.InputTag("fsrPhotonPFIsoChHadPU04"),
    cms.InputTag("rhoFsrPhoton"),
)

fsrPhotonSequence = cms.Sequence(
    fsrPhotons +
    fsrPhotonPFIsoChHad04 + 
    fsrPhotonPFIsoNHad04 + 
    fsrPhotonPFIsoPhoton04 + 
    fsrPhotonPFIsoChHadPU04 + 
    rhoFsrPhoton + 
    boostedFsrPhotons
)
