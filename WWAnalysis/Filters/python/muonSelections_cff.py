import FWCore.ParameterSet.Config as cms

muPlusSelection1 = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && abs(eta) < 2.4 && charge==+1"),
)

muMinusSelection1 = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && abs(eta) < 2.4 && charge==-1"),
)



# --- IP selections
import WWAnalysis.Filters.muIPSelection_cfi
muPlusIPSelection = WWAnalysis.Filters.muIPSelection_cfi.muIPSelection.clone(
  src = cms.InputTag("muPlusSelection1"),
)

muMinusIPSelection = WWAnalysis.Filters.muIPSelection_cfi.muIPSelection.clone(
  src = cms.InputTag("muMinusSelection1"),
)




# --- Iso selections
MUON_ISO_CUT=("(isolationR03().emEt +" +
            " isolationR03().hadEt +" +
            " isolationR03().sumPt)/pt < 0.15 ");

muPlusIsoSelection = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muPlusIPSelection","SelectedMuons"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ISO_CUT),
)

muMinusIsoSelection = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muMinusIPSelection","SelectedMuons"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ISO_CUT),
)


# --- ID selection
MUON_ID_CUT=("(isGlobalMuon && isTrackerMuon &&" +
             " innerTrack.found >10 &&" +
             " globalTrack.normalizedChi2 <10 &&" +
             " standAloneMuon.hitPattern.numberOfValidMuonHits > 0)" 
             );


muPlusID = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muPlusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ID_CUT),
)

muMinusID = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muMinusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ID_CUT),
)
                           

