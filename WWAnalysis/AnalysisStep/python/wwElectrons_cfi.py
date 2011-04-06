import FWCore.ParameterSet.Config as cms

        
from PhysicsTools.PatAlgos.patSequences_cff import *

ELE_ID_CUT=("((isEB && sigmaIetaIeta < 0.01 &&" +           
            " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004) ||" +
            " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
            " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007)) ")
#            " hadronicOverEm < 0.04) ||" +
#            " hadronicOverEm < 0.025 )) ")

ELE_ISO_CUT=("(( isEB && (dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.1) ||"+
            "( isEE && (dr03TkSumPt +" +
             " dr03EcalRecHitSumEt +" +
             " dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.1) ) ")

ELE_NOCONV_CUT=(" (abs(userFloat('convValueMapProd:dist')) > 0.02   ||" +
                "  abs(userFloat('convValueMapProd:dcot')) > 0.02 ) &&"
                "  gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0")


ELE_IP_CUT=(" abs(userFloat('dxyPV')) < 0.02 && " +
            " abs(userFloat('dzPV'))  < 1.0 " )

wwEleMatch = selectedPatElectrons.clone()
wwEleMatch.src = "boostedElectrons"
wwEleMatch.filter = cms.bool(False)
wwEleMatch.cut = ( "pt > 10 && abs(eta)<2.5 ")
#wwEleMatch.cut = ( "pt > 10 && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==11 && abs(genParticleRef.get().mother().mother().pdgId()) ==24")



wwEleID = selectedPatElectrons.clone()
wwEleID.src = "wwEleMatch"
wwEleID.filter = cms.bool(False)
wwEleID.cut = ( ELE_ID_CUT )

wwEleISO = selectedPatElectrons.clone()
wwEleISO.src = "wwEleID"
wwEleISO.filter = cms.bool(False)
wwEleISO.cut = ( ELE_ISO_CUT )

wwEleNOCONV = selectedPatElectrons.clone()
wwEleNOCONV.src = "wwEleISO"
wwEleNOCONV.filter = cms.bool(False)
wwEleNOCONV.cut = ( ELE_NOCONV_CUT )

wwEleIP = selectedPatElectrons.clone()
wwEleIP.src = "wwEleNOCONV"
wwEleIP.filter = cms.bool(False)
wwEleIP.cut = ( ELE_IP_CUT )



wwElectronSequence = cms.Sequence(  
    wwEleMatch*wwEleID*wwEleISO*wwEleNOCONV*wwEleIP
)


