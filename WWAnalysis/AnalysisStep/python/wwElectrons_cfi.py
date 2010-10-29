import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *

patEleIPSelector = cms.EDProducer("PatEleIPSelector",
  src = cms.InputTag("patElectronsWithTrigger"),
  vertexs = cms.InputTag("offlinePrimaryVertices"),
  d0Cut = cms.double(0.020),                              
)

#convRejSelection = cms.EDFilter("ConvRejectionSelectionPAT",
#  src            = cms.InputTag("patEleIPSelector"),
#  srcTracks      = cms.InputTag("generalTracks"),                               
#  filter         = cms.bool(False),
#  distMax        = cms.double(0.02),                              
#  dCotanThetaMax = cms.double(0.02),                              
#)


# --- Iso selections V1
ELE_ISO_CUT=("(( isEB && (dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt)/pt < 0.5) ||"+
             "( isEE && (dr03TkSumPt +" +
             " dr03EcalRecHitSumEt +" +
             " dr03HcalTowerSumEt)/pt < 0.5) ) ")

# --- Iso selections V0
ELE_ISO_CUT=("(dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt  )/pt < 0.10 ");






ELE_ID_CUT=("(" +
            "((isEB && sigmaIetaIeta < 0.01 &&" +           
            " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
            " hadronicOverEm < 0.04) ||" +
            " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
            " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 &&" +
            " hadronicOverEm < 0.025 )) &&" + 
            " abs(userFloat('convValueMapProd:dist')) > 0.02 &&" +
            " abs(userFloat('convValueMapProd:dcot')) > 0.02 &&" +
            " userInt('expectedHitsEle') == 0 )"
            )


wwElectrons = selectedPatElectrons.clone()
wwElectrons.src = "patEleIPSelector"

wwElectrons.filter = cms.bool(False)
wwElectrons.cut = ( "pt > 10 && " +
                    ELE_ID_CUT + " && " + 
                    ELE_ISO_CUT)

wwElectronSequence = cms.Sequence(
    patEleIPSelector*
    wwElectrons
)


