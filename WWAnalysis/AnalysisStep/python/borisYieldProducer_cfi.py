import FWCore.ParameterSet.Config as cms
import math

FWLiteParams = cms.PSet(
    maxEvents   =     cms.int32(-1),
    skipDuplicates   =     cms.bool(True),
    printEventNumber =     cms.bool(False),

    ### pre-selection cuts
    nMu            =     cms.int32(1),    #depends on mode (ELEL,ELMU,MUMU)
    nEl            =     cms.int32(1),    #depends on mode (ELEL,ELMU,MUMU)
    useNewExpectHits  =  cms.bool(False),
    ptMinHighPre   =     cms.double(20),
    ptMinLowPre    =     cms.double(10),
    minProjMetPre  =     cms.double(20),
    minMet         =     cms.double(20),
    minMll         =     cms.double(12),
    minDiffMz      =     cms.double(-1),  #depends on mode (ELEL,ELMU,MUMU)
    vetoSoftMuons  =     cms.bool(True),
    nCentralJet    =     cms.int32(0),
    jetVetoEt      =     cms.double(30),
    jetVetoEta     =     cms.double(5.0),


    ### final-selection cuts
    maxMll         =     cms.double(50),
    ptMinHigh      =     cms.double(30),
    ptMinLow       =     cms.double(25),
    minProjMet     =     cms.double(-1),
    minRatioMetPtLL =    cms.double(0.6),
    maxDPhiLLDegrees  =    cms.double(60.),

    #samples
    inputFolder =     cms.string(""),
    input       =     cms.vstring(),
    #fileOutName =     cms.string('output.root')

)
