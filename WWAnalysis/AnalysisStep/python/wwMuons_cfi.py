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

MUON_ISO_CUT=("(isolationR03().emEt +" +
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
                   " innerTrack.found >10 && abs(userFloat('dxyPV')) < 0.2 && " +
                   " ( (pt <= 20) || " +
                   "   (pt >  20  && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt" +
                   "                  - userFloat('rhoMu')*3.14159265*0.3*0.3)/pt > 0.10) ) )")

wwMuMatch = selectedPatMuons.clone()
wwMuMatch.src = "boostedMuons"
wwMuMatch.filter = cms.bool(False)
wwMuMatch.cut = ( "pt > 10 && abs(eta)<2.4")
#wwMuMatch.cut = ( "pt > 10 && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==13 && abs(genParticleRef.get().mother().mother().pdgId()) ==24")


wwMuonsID = selectedPatMuons.clone()
wwMuonsID.src = "wwMuMatch"
wwMuonsID.filter = cms.bool(False)
wwMuonsID.cut = ( MUON_ID_CUT )

wwMuonsISO = selectedPatMuons.clone()
wwMuonsISO.src = "wwMuonsID"
wwMuonsISO.filter = cms.bool(False)
wwMuonsISO.cut = ( MUON_ISO_CUT )

wwMuonsISOT = selectedPatMuons.clone()
wwMuonsISOT.src = "wwMuonsID"
wwMuonsISOT.filter = cms.bool(False)
wwMuonsISOT.cut = ( MUON_ISO_CUT_TIGHT )

wwMuonsISOPF = selectedPatMuons.clone()
wwMuonsISOPF.src = "wwMuonsID"
wwMuonsISOPF.filter = cms.bool(False)
wwMuonsISOPF.cut = ( MUON_ISOPF_CUT )

wwMuonsIP = selectedPatMuons.clone()
wwMuonsIP.src = "wwMuonsISO"
wwMuonsIP.filter = cms.bool(False)
wwMuonsIP.cut = ( MUON_IP_CUT )

wwMuonsIPT = selectedPatMuons.clone()
wwMuonsIPT.src = "wwMuonsISOT"
wwMuonsIPT.filter = cms.bool(False)
wwMuonsIPT.cut = ( MUON_IP_CUT )

wwMuonsIPPF = selectedPatMuons.clone()
wwMuonsIPPF.src = "wwMuonsISOPF"
wwMuonsIPPF.filter = cms.bool(False)
wwMuonsIPPF.cut = ( MUON_IP_CUT )

wwMuons4Veto = selectedPatMuons.clone()
wwMuons4Veto.src = "boostedMuons"
wwMuons4Veto.filter = cms.bool(False)
wwMuons4Veto.cut = ( "pt > 3 && " +
                     MUON_ID_CUT_4VETO )


wwMuonSequence = cms.Sequence( 
    wwMuMatch * 
    wwMuonsID * 
    wwMuonsISO * 
    wwMuonsIP * 
    wwMuonsISOT * 
    wwMuonsIPT * 
    wwMuonsISOPF * 
    wwMuonsIPPF * 
    wwMuons4Veto
)

