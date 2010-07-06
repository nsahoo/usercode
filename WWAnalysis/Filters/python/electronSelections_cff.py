import FWCore.ParameterSet.Config as cms

# --- initial selection
elePlusSelection1 = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("gsfElectrons"),
       filter = cms.bool(True),                              
#       cut = cms.string("pt> 20.0 && abs(eta) < 2.5 && charge==+1"),
       cut = cms.string("pt> 20.0 && eta > -2.5 && eta < 2.5 && charge==+1"),
)

eleMinusSelection1 = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("gsfElectrons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && eta > - 2.5 && eta < 2.5 && charge==-1"),
)




# --- IP selections
import WWAnalysis.Filters.eleIPSelection_cfi
elePlusIPSelection = WWAnalysis.Filters.eleIPSelection_cfi.eleIPSelection.clone(
  src = cms.InputTag("elePlusSelection1"),
)

eleMinusIPSelection = WWAnalysis.Filters.eleIPSelection_cfi.eleIPSelection.clone(
  src = cms.InputTag("eleMinusSelection1"),
)




# --- Iso selections
ELE_ISO_CUT=("(dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt  )/pt < 0.10 ");

elePlusIsoSelection = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("elePlusIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ISO_CUT),
)

eleMinusIsoSelection = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleMinusIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ISO_CUT),
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

elePlusID = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("elePlusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ID_CUT),
)

eleMinusID = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleMinusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ID_CUT),
)


# --- conversion rejection (FOR THE MOMENT IT IS JUST A DUMMY-SELECTOR/PLACEHOLDER) ---                          
convRejectionPlus = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("elePlusID"),
       filter = cms.bool(True),                              
       cut = cms.string("pt > 0.0"),
)

convRejectionMinus = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleMinusID"),
       filter = cms.bool(True),                              
       cut = cms.string("pt > 0.0"),
)



# --- diLepton selection
diEleSel1 = cms.EDProducer("CandViewShallowCloneCombiner",
    decay       = cms.string("convRejectionPlus convRejectionPlus"), 
    checkCharge = cms.bool(False),           
    cut         = cms.string("(mass > 12.0)"),
)

diEleFilter1 = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("diEleSel1"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


diEleSel2 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("diEleSel1"),
    filter = cms.bool(True),                              
    cut = cms.string("abs(mass-91.1876) > 15.0"),
)


# --- MET cuts
metSel1 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("tcMet"),
    filter = cms.bool(True),                              
    cut = cms.string("pt > 20.0"),
)

# --- soft muon veto
import WWAnalysis.Filters.softMuonVeto_cfi
eleSoftMuonVeto = WWAnalysis.Filters.softMuonVeto_cfi.softMuonVeto.clone(
  srcCompCands = cms.InputTag("diEleSel1"),
)
