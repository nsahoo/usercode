import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *


patMuIPSelector = cms.EDProducer("PatMuIPSelector",
#  src = cms.InputTag("patMuonsWithTrigger"), ## for new skim
  src = cms.InputTag("patMuons"), 
  vertexs = cms.InputTag("offlinePrimaryVertices"),
  d0Cut = cms.double(0.020),                              
)

ipSelMuon4Veto = cms.EDProducer("PatMuIPSelector",
#  src = cms.InputTag("patMuonsWithTrigger"), ## for new skim
  src = cms.InputTag("patMuons"),
  vertexs = cms.InputTag("offlinePrimaryVertices"),
  d0Cut = cms.double(0.20),                              
)


MUON_ISO_CUT=("(isolationR03().emEt +" +
            " isolationR03().hadEt +" +
            " isolationR03().sumPt)/pt < 0.15 ");

MUON_ID_CUT=("(isGlobalMuon && isTrackerMuon &&" +
             " innerTrack.found >10 &&" +
             " globalTrack.normalizedChi2 <10 &&" +
             " combinedMuon.hitPattern.numberOfValidMuonHits > 0)" 
             );

MUON_ID_CUT_4VETO=("(isTrackerMuon &&" +
                   " muonID('TMLastStationAngTight') &&" +
                   " innerTrack.found >10)");

wwMuons = selectedPatMuons.clone()
wwMuons.src = "patMuIPSelector"
wwMuons.filter = cms.bool(False)
wwMuons.cut = ( "pt > 10 && " +
                MUON_ID_CUT + " && " + 
                MUON_ISO_CUT)

wwMuons4Veto = selectedPatMuons.clone()
wwMuons4Veto.src = "ipSelMuon4Veto"
wwMuons4Veto.filter = cms.bool(False)
wwMuons4Veto.cut = ( "pt > 3 && " +
                     MUON_ID_CUT_4VETO )


wwMuonSequence = cms.Sequence(
    patMuIPSelector*ipSelMuon4Veto*
    wwMuons*wwMuons4Veto
)

