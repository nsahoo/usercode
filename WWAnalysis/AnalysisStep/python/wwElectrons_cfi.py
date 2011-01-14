import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *

boostedElectrons = cms.EDProducer("PatElectronBooster",
  electronTag = cms.untracked.InputTag("patElectronsWithTrigger"),
  vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
  electronIDSources = cms.PSet(
      myEidLoose = cms.InputTag("myEidLoose"),
      myEidMedium = cms.InputTag("myEidMedium"),
      myEidTight = cms.InputTag("myEidTight"),
  ),
)

from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentification_cfi import eidLoose as myEidLoose
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentification_cfi import eidMedium as myEidMedium
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentification_cfi import eidTight as myEidTight
myEidLoose.src = 'patElectronsWithTrigger'
myEidMedium.src = 'patElectronsWithTrigger'
myEidTight.src = 'patElectronsWithTrigger'

myEIDSequence = cms.Sequence(
    myEidLoose +
    myEidMedium +
    myEidTight
)


# --- Iso selections V1
ELE_ISO_CUT=("(( isEB && (dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt)/pt < 0.10) ||"+
             "( isEE && (dr03TkSumPt +" +
             " dr03EcalRecHitSumEt +" +
             " dr03HcalTowerSumEt)/pt < 0.10) ) ")

# --- Iso selections V0
# ELE_ISO_CUT=("(dr03TkSumPt +" +
#              " max(0,dr03EcalRecHitSumEt - 1.0) +" +
#              " dr03HcalTowerSumEt  )/pt < 0.10 ");

ELE_ID_CUT=("(" +
            "((isEB && sigmaIetaIeta < 0.01 &&" +           
            " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
            " hadronicOverEm < 0.04) ||" +
            " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
            " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 &&" +
            " hadronicOverEm < 0.025 )) &&" + 
            " (abs(userFloat('convValueMapProd:dist')) > 0.02   ||" +
            "  abs(userFloat('convValueMapProd:dcot')) > 0.02 ) &&" +
            " userInt('expectedHitsEle') == 0 && " + 
            " abs(userFloat('dxyPV')) < 0.02 && " +
            " abs(userFloat('dzPV'))  < 1.0  )  " )


wwElectrons = selectedPatElectrons.clone()
wwElectrons.src = "boostedElectrons"
wwElectrons.filter = cms.bool(False)
wwElectrons.cut = ( "pt > 10 && " +
                    ELE_ID_CUT + " && " + 
                    ELE_ISO_CUT)

wwElectronSequence = cms.Sequence(
    myEIDSequence *
    boostedElectrons *
    wwElectrons
)


