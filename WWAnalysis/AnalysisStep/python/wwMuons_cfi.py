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


MUON_IP_CUT=( "( abs(userFloat('dxyPV')) < 0.02 && " +
              "  abs(userFloat('dzPV'))  < 1.0    )" )


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


wwMuons = selectedPatMuons.clone()
wwMuons.src = "wwMuMatch"
wwMuons.filter = cms.bool(False)
wwMuons.cut = ( MUON_ID_CUT + "&&" + MUON_ISO_CUT + "&&" + MUON_IP_CUT )


wwMuons4Veto = selectedPatMuons.clone()
wwMuons4Veto.src = "boostedMuons"
wwMuons4Veto.filter = cms.bool(False)
wwMuons4Veto.cut = ( "pt > 3 && " +
                     MUON_ID_CUT_4VETO )


wwMuonSequence = cms.Sequence( 
    wwMuMatch*wwMuons*wwMuons4Veto
)

