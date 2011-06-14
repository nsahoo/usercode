import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *

MUON_ID_CUT=("(( (isGlobalMuon() && "
             "    globalTrack.normalizedChi2 <10 &&" +
             "    globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
             "    numberOfMatches > 1 ) || " + 
             "   (isTrackerMuon() && muonID('TMLastStationTight')) ) && " + 
             " innerTrack.found >10 &&" +
             " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
             " abs(track.ptError / pt) < 0.10 )")

MUON_ID_CUT_OLD=("(isGlobalMuon && isTrackerMuon &&" +
                 " innerTrack.found >10 &&" +
                 " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
                 " globalTrack.normalizedChi2 <10 &&" +
                 " globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
                 " numberOfMatches > 1 && " + 
                 " abs(track.ptError / pt) < 0.10 )")

#SMURF_ISO = ("( userFloat('smurfCharged') + userFloat('smurfPhoton') + userFloat('smurfNeutral') )/ pt")
SMURF_ISO = ("( userFloat('muSmurfPF') )/ pt")
MUON_MERGE_ISO  =   ("( (abs(eta) < 1.479 && pt >  20 && " + SMURF_ISO + " < 0.13) || ( abs(eta) >= 1.479 && pt >  20 && " + SMURF_ISO + " < 0.09 ) || " + 
                     "  (abs(eta) < 1.479 && pt <= 20 && " + SMURF_ISO + " < 0.06) || ( abs(eta) >= 1.479 && pt <= 20 && " + SMURF_ISO + " < 0.05 ) )  ")

MUON_MERGE_IP  = ("( ( (pt >= 20 && abs(userFloat('tip')) < 0.02) || (pt < 20 && abs(userFloat('tip')) < 0.01) ) && " +
                  "  abs(userFloat('dzPV'))  < 0.1 )" )
              
              


MUON_ISO_CUT = ("(isolationR03().emEt +" +
                " isolationR03().hadEt +" +
                " isolationR03().sumPt - userFloat('rhoMu')*3.14159265*0.3*0.3)/pt < 0.15 ");

MUON_ISO_CUT_TIGHT=("( ( pt > 20 && (isolationR03().emEt + isolationR03().hadEt + " +
                        " isolationR03().sumPt - userFloat('rhoMu')*3.14159265*0.3*0.3)/pt < 0.15 ) || " + 
                    "  ( pt <= 20 && (isolationR03().emEt + isolationR03().hadEt +" +
                        " isolationR03().sumPt - userFloat('rhoMu')*3.14159265*0.3*0.3)/pt < 0.12 ) )");

MUON_ISOPF_CUT = ("( (userFloat('pfCharged')+userFloat('pfPhoton')+userFloat('pfNeutral')-userFloat('rhoMuNoPU')*3.14159265*0.4*0.4) / pt < 0.20)")


MUON_IP_CUT=( "( abs(userFloat('tip2')) < 0.01 && " +
              "  abs(userFloat('dzPV'))  < 0.05    )" )


MUON_ID_CUT_4VETO=("(isTrackerMuon &&" +
                   " muonID('TMLastStationAngTight') &&" +
                   " innerTrack.found >10 && abs(userFloat('tip')) < 0.2 && abs(userFloat('dzPV')) < 0.1 &&" +
                   " ( (pt <= 20) || " +
                   "   (pt >  20  && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt" +
                   "                 )/pt > 0.10) ) )")

selectedRefPatMuons = cms.EDFilter("PATMuonViewRefSelector",
   src = cms.InputTag("input"),
   cut = cms.string("")
)

wwMuMatch = selectedRefPatMuons.clone()
wwMuMatch.src = "boostedMuons"
wwMuMatch.filter = cms.bool(False)
wwMuMatch.cut = ( "pt > 10 && abs(eta)<2.4")
#wwMuMatch.cut = ( "pt > 10 && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==13 && abs(genParticleRef.get().mother().mother().pdgId()) ==24")


wwMuonsID = selectedRefPatMuons.clone()
wwMuonsID.src = "wwMuMatch"
wwMuonsID.filter = cms.bool(False)
wwMuonsID.cut = ( MUON_ID_CUT )

wwMuonsISO = selectedRefPatMuons.clone()
wwMuonsISO.src = "wwMuonsID"
wwMuonsISO.filter = cms.bool(False)
wwMuonsISO.cut = ( MUON_ISO_CUT )

wwMuonsISOT = selectedRefPatMuons.clone()
wwMuonsISOT.src = "wwMuonsID"
wwMuonsISOT.filter = cms.bool(False)
wwMuonsISOT.cut = ( MUON_ISO_CUT_TIGHT )

wwMuonsMergeISO = selectedRefPatMuons.clone()
wwMuonsMergeISO.src = "wwMuonsID"
wwMuonsMergeISO.filter = cms.bool(False)
wwMuonsMergeISO.cut = ( MUON_MERGE_ISO )

wwMuonsISOPF = selectedRefPatMuons.clone()
wwMuonsISOPF.src = "wwMuonsID"
wwMuonsISOPF.filter = cms.bool(False)
wwMuonsISOPF.cut = ( MUON_ISOPF_CUT )

wwMuonsIP = selectedRefPatMuons.clone()
wwMuonsIP.src = "wwMuonsISO"
wwMuonsIP.filter = cms.bool(False)
wwMuonsIP.cut = ( MUON_IP_CUT )

wwMuonsIPT = selectedRefPatMuons.clone()
wwMuonsIPT.src = "wwMuonsISOT"
wwMuonsIPT.filter = cms.bool(False)
wwMuonsIPT.cut = ( MUON_IP_CUT )

wwMuonsMergeIP = selectedRefPatMuons.clone()
wwMuonsMergeIP.src = "wwMuonsMergeISO"
wwMuonsMergeIP.filter = cms.bool(False)
wwMuonsMergeIP.cut = ( MUON_MERGE_IP )

wwMuonsIPPF = selectedRefPatMuons.clone()
wwMuonsIPPF.src = "wwMuonsISOPF"
wwMuonsIPPF.filter = cms.bool(False)
wwMuonsIPPF.cut = ( MUON_IP_CUT )

wwMuons4Veto = selectedRefPatMuons.clone()
wwMuons4Veto.src = "boostedMuons"
wwMuons4Veto.filter = cms.bool(False)
wwMuons4Veto.cut = ( "pt > 3 && " +
                     MUON_ID_CUT_4VETO )


wwMuonSequence = cms.Sequence( 
    wwMuMatch * 
    wwMuonsID * 
    wwMuonsISO * 
    wwMuonsIP * 
    wwMuonsMergeISO * 
    wwMuonsMergeIP * 
#     wwMuonsISOT * 
#     wwMuonsIPT * 
#     wwMuonsISOPF * 
#     wwMuonsIPPF * 
    wwMuons4Veto
)

