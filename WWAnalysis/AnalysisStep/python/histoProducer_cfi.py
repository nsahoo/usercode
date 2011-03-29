import FWCore.ParameterSet.Config as cms
import math

FWLiteParams = cms.PSet(
    maxEvents   =     cms.int32(-1),

    ### pre-selection cuts
    nMu            =     cms.int32(1),    #depends on mode (ELEL,ELMU,MUMU)
    nEl            =     cms.int32(1),    #depends on mode (ELEL,ELMU,MUMU)
    #useNewExpectHits  =  cms.bool(True),
    useNewExpectHits  =  cms.bool(False),
    ptMinHighPre   =     cms.double(20),
    ptMinLowPre    =     cms.double(10),
    minProjMetPre  =     cms.double(20),
    minMll         =     cms.double(12),
    minDiffMz      =     cms.double(-1),  #depends on mode (ELEL,ELMU,MUMU)
    vetoSoftMuons  =     cms.bool(True),
    nCentralJet    =     cms.int32(0),
    jetVetoEt      =     cms.double(25),
    jetVetoEta     =     cms.double(5.0),


    ### final-selection cuts
    maxMll         =     cms.double(50),
    ptMinHigh      =     cms.double(30),
    ptMinLow       =     cms.double(28),
    minProjMet     =     cms.double(-1),
    minRatioMetPtLL =    cms.double(0.6),
    maxDPhiLL      =    cms.double(math.pi/4.0),

    #samples
    inputFolder =     cms.string('/nfs/bluearc/group/trees/hwwme2nu/skimEvent_'), #depends on mode (ELEL,ELMU,MUMU)
    input       =     cms.vstring('file:input1.root','file:input2.root'),
    fileOutName =     cms.string('output.root')

)
