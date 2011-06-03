import FWCore.ParameterSet.Config as cms

def addRhoVariables(process,seq,eleTag='gsfElectrons',muTag='muons',pfNoPUTag='pfNoPileUp'):
    process.load('RecoJets.JetProducers.kt4PFJets_cfi')
    
    process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
    process.kt6PFJets.Rho_EtaMax = cms.double(4.5)
    
    process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
    process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)
    process.kt6PFJetsForIso.Ghost_EtaMax = cms.double(2.5)
    
    process.kt6PFJetsNoPU = process.kt6PFJets.clone( src = pfNoPUTag )
    process.kt6PFJetsForIsoNoPU = process.kt6PFJetsForIso.clone( src = pfNoPUTag )
    
    # Re-cluster ak5PFJets w/ Area calculation on
    process.ak5PFJets.doAreaFastjet = True
    process.ak5PFJets.Rho_EtaMax = cms.double(4.5)
    
    # Re-cluster jets w/ pfNoPileUp
    process.ak5PFJetsNoPU = process.ak5PFJets.clone( src = pfNoPUTag )
    
    process.load("WWAnalysis.Tools.betaValueMapProducer_cfi")
    process.load("WWAnalysis.Tools.rhoValueMapProducer_cfi")
    
    process.betaMu = process.betaValueMapProducer.clone()
    process.betaEl = process.betaValueMapProducer.clone()
    process.betaEl.leptonTag = eleTag
    process.betaEl.dEtaVeto = 0.015
    process.betaEl.dRVeto = 0.0
    
    process.rhoMu = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIso","rho",process.name_()), leptonTag = muTag)
    process.rhoEl = process.rhoMu.clone(leptonTag = eleTag)
    process.rhoMuNoPU = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIsoNoPU","rho",process.name_()), leptonTag = muTag)
    process.rhoElNoPU = process.rhoMuNoPU.clone(leptonTag = eleTag)
    
    print "="*80
    print "Still need to add DeltaBeta!!!"
    print "="*80
    process.valueMaps = cms.Sequence(
        process.betaMu +
        process.betaEl +
        process.rhoMu +
        process.rhoEl +
        process.rhoMuNoPU +
        process.rhoElNoPU 
    )

    seq += (
        process.ak5PFJets + 
        process.kt6PFJets +
        process.kt6PFJetsForIso + 
        process.pfPileUp *
        process.pfNoPileUp * (
            process.ak5PFJetsNoPU +
            process.kt6PFJetsNoPU + 
            process.kt6PFJetsForIsoNoPU ) *
        process.valueMaps 
    )
        
