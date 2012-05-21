import FWCore.ParameterSet.Config as cms

muons4skim = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string("(isTrackerMuon||isGlobalMuon) && abs(eta) < 2.4"),
)
electrons4skim = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string("abs(eta) < 2.5"),
)
leptons4skim = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag( cms.InputTag("muons4skim"),  cms.InputTag("electrons4skim"), )
)
dileptons4skim = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('leptons4skim leptons4skim'),
    cut = cms.string('deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.01'),
    checkCharge = cms.bool(False)
)

skim2010 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('min(daughter(0).pt,daughter(1).pt) > 10 && ' +
                     'max(daughter(0).pt,daughter(1).pt) > 20'),
    filter = cms.bool(True),
)

skim40 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('mass > 40 && daughter(0).pdgId = -daughter(1).pdgId'), ## and SF/OS
    filter = cms.bool(True),
)
skim40NoOF  = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('mass > 40 && abs(daughter(0).pdgId) == abs(daughter(1).pdgId)'), ## and SF only
    filter = cms.bool(True),
)


reskim     = cms.Sequence(muons4skim + electrons4skim + leptons4skim + dileptons4skim + skim2010 + skim40)
reskimNoOS = cms.Sequence(muons4skim + electrons4skim + leptons4skim + dileptons4skim + skim2010 + skim40NoOF)

from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *
triggerResultsFilter.l1tResults = ''
triggerResultsFilter.throw = True
triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
triggerFilterDoubleElectron7TeV = triggerResultsFilter.clone(triggerConditions = [ 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*' ])
triggerFilterDoubleMuon7TeV     = triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu13_Mu8_v*', 'HLT_Mu17_Mu8_v*' ])
triggerFilter7TeV = triggerResultsFilter.clone(triggerConditions = [ 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 'HLT_Mu13_Mu8_v*', 'HLT_Mu17_Mu8_v*' ])
triggerFilter7TeV_77  = triggerResultsFilter.clone(triggerConditions = [ 'HLT_DoubleMu7_v*', 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*' ])
triggerFilter7TeV_138 = triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu13_Mu8_v*',  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*' ])
triggerFilter7TeV_178 = triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu17_Mu8_v*',  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*' ])

### Trigger emulation
triggerEmulator7TeV_178 = cms.EDFilter("GioTriggerEmulator",
    trigger = cms.InputTag("TriggerResults","","HLT"),
    muons   = cms.InputTag("boostedMuons"),
    doubleMu = cms.string("Mu17_Mu8"), # or "Mu13_Mu8" or "none" or "any" (=loosest unprescaled)
    doubleEl = cms.string("Ele17_Ele8_Tight"), # or "any" or "none"
    runForMC = cms.uint32(1), # run whose trigger table should be assumed in MC (180000 is good for Fall11 MC; use 1 or 178419 to test emulator on DoubleMu7 or Mu13/Mu7)
)