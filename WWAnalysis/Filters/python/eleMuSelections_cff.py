import FWCore.ParameterSet.Config as cms

# --- initial selection
eleSelection1 = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("gsfElectrons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && eta > -2.5 && eta < 2.5"),
)

muSelection1 = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && abs(eta) < 2.4"),
)

selection1Pairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("eleSelection1@+ muSelection1@-"),
    cut = cms.string("")
)

selection1PairsCounter = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("selection1Pairs"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


# --- IP selections
import WWAnalysis.Filters.eleIPSelection_cfi
eleIPSelection = WWAnalysis.Filters.eleIPSelection_cfi.eleIPSelection.clone(
  src = cms.InputTag("eleSelection1"),
)

import WWAnalysis.Filters.muIPSelection_cfi
muIPSelection = WWAnalysis.Filters.muIPSelection_cfi.muIPSelection.clone(
  src = cms.InputTag("muSelection1"),
)


# --- Iso selections
ELE_ISO_CUT=("(dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt  )/pt < 0.10 ");

eleIsoSelection = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ISO_CUT),
)

MUON_ISO_CUT=("(isolationR03().emEt +" +
            " isolationR03().hadEt +" +
            " isolationR03().sumPt)/pt < 0.15 ");

muIsoSelection = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muIPSelection","SelectedMuons"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ISO_CUT),
)



# --- ID selection
ELE_ID_CUT=("( (isEB && sigmaIetaIeta < 0.01 &&" +           
            " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
            " hadronicOverEm < 0.04) ||" +
            " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
            " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 &&" +
            " hadronicOverEm < 0.025) )"
            );

eleID = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ID_CUT),
)

MUON_ID_CUT=("(isGlobalMuon && isTrackerMuon &&" +
             " innerTrack.found >10 &&" +
             " globalTrack.normalizedChi2 <10 &&" +
             " standAloneMuon.hitPattern.numberOfValidMuonHits > 0)" 
             );

muID = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ID_CUT),
)


# --- conversion rejection  ---    
import WWAnalysis.Filters.convRejectionSelection_cfi      
eleConvRejection = WWAnalysis.Filters.convRejectionSelection_cfi.convRejectionSelection.clone(
  src = cms.InputTag("eleID"),
)


