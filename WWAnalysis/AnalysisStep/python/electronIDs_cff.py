
import FWCore.ParameterSet.Config as cms

ELE_ID_CB_85_2011=("((isEB && sigmaIetaIeta < 0.01 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.039 && deltaPhiSuperClusterTrackAtVtx < 0.039 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.005 && deltaEtaSuperClusterTrackAtVtx < 0.005) || " + 
                   " ( (!isEB) && sigmaIetaIeta < 0.031  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.028 && deltaPhiSuperClusterTrackAtVtx < 0.028 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 )) ")

ELE_ID_CB_80_2011=("((isEB && sigmaIetaIeta < 0.01 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.027 && deltaPhiSuperClusterTrackAtVtx < 0.027 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.005 && deltaEtaSuperClusterTrackAtVtx < 0.005) || " + 
                   " ( (!isEB) && sigmaIetaIeta < 0.031  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.021 && deltaPhiSuperClusterTrackAtVtx < 0.021 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.006 && deltaEtaSuperClusterTrackAtVtx < 0.006 )) ")

ELE_ID_CB_80_2010=("((isEB && sigmaIetaIeta < 0.01 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
                   " hadronicOverEm < 0.04) ||" +
                   " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 " + 
                   " hadronicOverEm < 0.025 )) ")

ELE_ISO_CB_85_2011=("(( isEB && (dr03TkSumPt +" +
                 "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                 "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.053) ||"+
                 "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                 "  dr03EcalRecHitSumEt +" +
                 "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.042) ) ")

ELE_ISO_CB_80_2011=("(( isEB && (dr03TkSumPt +" +
                 "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                 "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.040) ||"+
                 "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                 "  dr03EcalRecHitSumEt +" +
                 "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.033) ) ")

ELE_ISO_CB_80_2010=("(( isEB && (dr03TkSumPt +" +
                 "  max(0,dr03EcalRecHitSumEt - 1.0) +" +
                 "  dr03HcalTowerSumEt)/pt < 0.1) ||"+
                 "( isEE && (dr03TkSumPt +" +
                 "  dr03EcalRecHitSumEt +" +
                 "  dr03HcalTowerSumEt)/pt < 0.1) ) ")

ELE_NOCONV =(" (abs(userFloat('convValueMapProd:dist')) > 0.02   ||" +
             "  abs(userFloat('convValueMapProd:dcot')) > 0.02 ) &&"
             "  gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0")


ELE_IP =(" abs(userFloat('dxyPV')) < 0.02 && " +
            " abs(userFloat('dzPV'))  < 1.0 " )

