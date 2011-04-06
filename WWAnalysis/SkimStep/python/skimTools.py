import FWCore.ParameterSet.Config as cms

def addIsolationInformation(process):
    print "*********************************************************************************"
    print "WARNING:  Adding isolation information!  This means you MUST run on RECO"
    print "*********************************************************************************"

    process.load("RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi")
    process.muIsoDepositCalByAssociatorTowers.IOPSet.inputMuonCollection = 'boostedMuons' 
    
    process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequencePAT_cff")
    process.eleIsoDepositEcalFromHits.src = 'boostedElectrons'
    process.eleIsoDepositHcalFromTowers.src = 'boostedElectrons'
    process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE")
    process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB")
    process.eleIsoDepositEcalFromHits.ExtractorPSet.energyMin = 0.0

    process.isoSequence = cms.Sequence(
        process.eleIsoDepositEcalFromHits +
        process.eleIsoDepositHcalFromTowers +
        process.muIsoDepositCalByAssociatorTowers
    )
    process.postPatSequence += process.isoSequence

    process.out.outputCommands.append('keep recoTracks_generalTracks_*_*')
    process.out.outputCommands.append('keep *_eleIsoDepositEcalFromHits_*_'+process.name_())
    process.out.outputCommands.append('keep *_eleIsoDepositHcalFromTowers_*_'+process.name_())
    process.out.outputCommands.append('keep *_muIsoDepositCalByAssociatorTowers_ecal_'+process.name_())
    process.out.outputCommands.append('keep *_muIsoDepositCalByAssociatorTowers_hcal_'+process.name_())
    process.out.outputCommands.append('keep recoPFCandidates_particleFlow__RECO')

