import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *

boostedMuons = cms.EDProducer("PatMuonBooster",
  muonTag = cms.untracked.InputTag("patMuonsWithTrigger"), 
  vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
)

MUON_ISO_CUT=("(isolationR03().emEt +" +
            " isolationR03().hadEt +" +
            " isolationR03().sumPt)/pt < 0.15 ");

MUON_ID_CUT=("(isGlobalMuon && isTrackerMuon &&" +
             " innerTrack.found >10 &&" +
             " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
             " globalTrack.normalizedChi2 <10 &&" +
             " globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
             " numberOfMatches > 1 && " + 
             " abs(track.ptError / pt) < 0.10 && " +
             " abs(userFloat('dxyPV')) < 0.02 && " +
             " abs(userFloat('dzPV'))  < 1.0  )  " )


MUON_ID_CUT_4VETO=("(isTrackerMuon &&" +
                   " muonID('TMLastStationAngTight') &&" +
                   " innerTrack.found >10 && abs(userFloat('dxyPV')) < 0.2 && " +
                   " ( (pt <= 20) || " +
                   "   (pt >  20  && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt)/pt > 0.10) ) )")

wwMuons = selectedPatMuons.clone()
wwMuons.src = "boostedMuons"
wwMuons.filter = cms.bool(False)
wwMuons.cut = ( "pt > 10 && " +
                MUON_ID_CUT + " && " + 
                MUON_ISO_CUT)

wwMuons4Veto = selectedPatMuons.clone()
wwMuons4Veto.src = "boostedMuons"
wwMuons4Veto.filter = cms.bool(False)
wwMuons4Veto.cut = ( "pt > 3 && " +
                     MUON_ID_CUT_4VETO )


wwMuonSequence = cms.Sequence(
    boostedMuons*
    wwMuons+wwMuons4Veto
)

