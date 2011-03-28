import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *

#boostedMuons = cms.EDProducer("PatMuonBooster",
#  muonTag = cms.untracked.InputTag("cleanPatMuonsTriggerMatch"), 
#  vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
#  deposits = cms.VPSet()
#)

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

# tempPSet =  cms.PSet(
#         mode = cms.string('sum'),
#         src = cms.InputTag("muIsoDepositTk"),
#         weight = cms.string('1'),
#         deltaR = cms.double(0.3),
#         vetos = cms.vstring(),
#         vetos = cms.vstring('0.01','Threshold(0.7)'),
#         skipDefaultVeto = cms.bool(False),
#         label = cms.string('isoTk')
# ) 

# boostedMuons.deposits.append( tempPSet.clone() )
# boostedMuons.deposits[-1].label = cms.string('isoTk')
# boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositTk')
# boostedMuons.deposits.append( tempPSet.clone() )
# boostedMuons.deposits[-1].label = cms.string('isoEcal')
# boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','ecal')
# boostedMuons.deposits.append( tempPSet.clone() )
# boostedMuons.deposits[-1].label = cms.string('isoHcal')
# boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','hcal')


# from RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi import muIsoDepositTk as myMuIsoDepositTk
# from RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi import muIsoDepositCalByAssociatorTowers as myMuIsoDepositCalByAssociatorTowers
# myMuIsoDepositTk.IOPSet.inputMuonCollection = 'patMuonsWithTrigger'
# myMuIsoDepositCalByAssociatorTowers.IOPSet.inputMuonCollection = 'patMuonsWithTrigger' 

wwMuonSequence = cms.Sequence( 
    #boostedMuons*
    wwMuons+wwMuons4Veto
)

