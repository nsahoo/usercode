import FWCore.ParameterSet.Config as cms
import math

def swapOutPSetValues(taker,giver):
    for x in giver.parameterNames_():
        setattr(taker,x,getattr(giver,x))


#can be used for wwelel or wwmumu
defaultWW = cms.PSet(
    ptMax       = cms.double(20),
    etaMu       = cms.double(2.4),
    etaEl       = cms.double(2.5),
    d0          = cms.double(0.020),
    dZ          = cms.double(1.0),
    isoEl       = cms.double(0.10),
    isoMu       = cms.double(0.15),
    met         = cms.double(20.0),
    mll         = cms.double(12.0),
    jetPt       = cms.double(25),
    jetEta      = cms.double(5.0),
    nCentralJet = cms.int32(0),
    nSoftMu     = cms.int32(0),
    bValue      = cms.double(2.1),
    nBtagJets   = cms.int32(0),

    #below needs to change for emu channel
    pMet        = cms.double(35.0),
    mZ          = cms.double(15.0),

    #below needs to change for higgs
    nExtraLep   = cms.int32(0),
    ptMin       = cms.double(10),
    useJEC      = cms.bool(True),

    #final selection cuts follow for higgs
    #these are currenltly set for the ww selection
    mllMaxFinal = cms.double(9999.0),
    ptMaxFinal  = cms.double(0.0),
    ptMinFinal  = cms.double(0.0),
    pMetFinal   = cms.double(0.0),
    metOverPtLL = cms.double(0.0),
    deltaPhiLL  = cms.double(2*math.pi),
)

#swap this out for wwelmu
oppositeFlavor = cms.PSet(
    pMet        = cms.double(20.0),
    pMetFinal   = cms.double(20.0),
    mZ          = cms.double(-1.0),
)


h120 = cms.PSet( ptMin=cms.double(15), useJEC=cms.bool(True), mllMaxFinal=cms.double( 40), ptMaxFinal=cms.double( 20), ptMinFinal=cms.double(15), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 60./180*math.pi))  
h130 = cms.PSet( ptMin=cms.double(15), useJEC=cms.bool(True), mllMaxFinal=cms.double( 45), ptMaxFinal=cms.double( 25), ptMinFinal=cms.double(15), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 60./180*math.pi))  
h140 = cms.PSet( ptMin=cms.double(15), useJEC=cms.bool(True), mllMaxFinal=cms.double( 45), ptMaxFinal=cms.double( 25), ptMinFinal=cms.double(15), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 60./180*math.pi))  
h150 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 50), ptMaxFinal=cms.double( 27), ptMinFinal=cms.double(20), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 50./180*math.pi))  
h160 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 50), ptMaxFinal=cms.double( 30), ptMinFinal=cms.double(28), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 45./180*math.pi))  
h170 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 50), ptMaxFinal=cms.double( 35), ptMinFinal=cms.double(28), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 50./180*math.pi))  
h180 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 60), ptMaxFinal=cms.double( 37), ptMinFinal=cms.double(28), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 70./180*math.pi))  
h190 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 80), ptMaxFinal=cms.double( 39), ptMinFinal=cms.double(28), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double( 90./180*math.pi))  
h200 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 90), ptMaxFinal=cms.double( 40), ptMinFinal=cms.double(28), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(100./180*math.pi))  
h250 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(120), ptMaxFinal=cms.double( 48), ptMinFinal=cms.double(28), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(120./180*math.pi))    
h300 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(140), ptMaxFinal=cms.double( 55), ptMinFinal=cms.double(30), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(140./180*math.pi))  
h350 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(170), ptMaxFinal=cms.double( 60), ptMinFinal=cms.double(30), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(140./180*math.pi))    
h400 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(200), ptMaxFinal=cms.double( 65), ptMinFinal=cms.double(30), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(145./180*math.pi))  
h450 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(250), ptMaxFinal=cms.double( 80), ptMinFinal=cms.double(30), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(160./180*math.pi))  
h500 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(300), ptMaxFinal=cms.double( 90), ptMinFinal=cms.double(30), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(170./180*math.pi))  
h550 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(350), ptMaxFinal=cms.double(110), ptMinFinal=cms.double(30), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(170./180*math.pi))  
h600 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(400), ptMaxFinal=cms.double(130), ptMinFinal=cms.double(30), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.6), deltaPhiLL=cms.double(170./180*math.pi))  

hReOpt120 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 40), ptMaxFinal=cms.double( 20), ptMinFinal=cms.double(20), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 60./180*math.pi))  
hReOpt130 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 45), ptMaxFinal=cms.double( 25), ptMinFinal=cms.double(20), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 60./180*math.pi))  
hReOpt140 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 45), ptMaxFinal=cms.double( 25), ptMinFinal=cms.double(20), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 60./180*math.pi))  
hReOpt150 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 50), ptMaxFinal=cms.double( 27), ptMinFinal=cms.double(20), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 60./180*math.pi))  
hReOpt160 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 50), ptMaxFinal=cms.double( 30), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 60./180*math.pi))  
hReOpt170 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 50), ptMaxFinal=cms.double( 34), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 60./180*math.pi))  
hReOpt180 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 60), ptMaxFinal=cms.double( 36), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 70./180*math.pi))  
hReOpt190 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 80), ptMaxFinal=cms.double( 38), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double( 90./180*math.pi))  
hReOpt200 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double( 90), ptMaxFinal=cms.double( 40), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(100./180*math.pi))  
hReOpt210 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(110), ptMaxFinal=cms.double( 44), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(110./180*math.pi))  
hReOpt220 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(120), ptMaxFinal=cms.double( 48), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(120./180*math.pi))  
hReOpt230 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(130), ptMaxFinal=cms.double( 52), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(130./180*math.pi))  
hReOpt250 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(150), ptMaxFinal=cms.double( 55), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(140./180*math.pi))    
hReOpt300 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(200), ptMaxFinal=cms.double( 70), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(175./180*math.pi))  
hReOpt350 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(250), ptMaxFinal=cms.double( 80), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(175./180*math.pi))    
hReOpt400 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(300), ptMaxFinal=cms.double( 90), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(175./180*math.pi))  
hReOpt450 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(350), ptMaxFinal=cms.double(110), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(175./180*math.pi))  
hReOpt500 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(400), ptMaxFinal=cms.double(120), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(175./180*math.pi))  
hReOpt550 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(450), ptMaxFinal=cms.double(130), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(175./180*math.pi))  
hReOpt600 = cms.PSet( ptMin=cms.double(20), useJEC=cms.bool(True), mllMaxFinal=cms.double(500), ptMaxFinal=cms.double(140), ptMinFinal=cms.double(25), pMetFinal=cms.double(35), metOverPtLL=cms.double(0.0), deltaPhiLL=cms.double(175./180*math.pi))  

