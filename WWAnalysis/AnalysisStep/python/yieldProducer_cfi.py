import FWCore.ParameterSet.Config as cms
import math

FWLiteParams = cms.PSet(
    ptMin       = cms.double(20),
    ptMax       = cms.double(20),
    etaMu       = cms.double(2.4),
    etaEl       = cms.double(2.5),
    d0          = cms.double(0.020),
    dZ          = cms.double(1.0),
    isoEl       = cms.double(0.10),
    isoMu       = cms.double(0.15),
    met         = cms.double(20.0),
    mll         = cms.double(12.0),
    mZ          = cms.double(15.0),
    pMet        = cms.double(20.0),
    jetPt       = cms.double(25),
    jetEta      = cms.double(5.0),
    nCentralJet = cms.int32(0),
    nSoftMu     = cms.int32(0),
    nExtraLep   = cms.int32(0),
    bValue      = cms.double(2.1),
    nBtagJets   = cms.int32(0),

    inputFolder = cms.string('/nfs/bluearc/group/trees/hwwme2nu/skimEvent_'), #depends on mode (ELEL,ELMU,MUMU)
    input       = cms.vstring('input1','input2'), #.root automatically appended
    fileOutName = cms.string('output.root')
)
