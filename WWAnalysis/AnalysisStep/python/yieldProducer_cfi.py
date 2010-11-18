import FWCore.ParameterSet.Config as cms
import math

FWLiteParams = cms.PSet(
    selectionParams = cms.PSet(
        wwelel = cms.PSet(
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
            pMet        = cms.double(35.0),
            jetPt       = cms.double(25),
            jetEta      = cms.double(5.0),
            nCentralJet = cms.int32(0),
            nSoftMu     = cms.int32(0),
            nExtraLep   = cms.int32(0),
            bValue      = cms.double(2.1),
            nBtagJets   = cms.int32(0),
        ),
    ),


    inputParams = cms.PSet(
#         ww = cms.PSet( 
#             files = cms.vstring( '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd/id13.root' ),
#             scale = cms.double(0.001428),
#             printEvents = cms.bool(False),
#         ),
    ),

    l2File      = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2Relative.txt"),
    l3File      = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L3Absolute.txt"),
    resFile     = cms.string(""),  #set to nothing if running on MC
#     resFile     = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2L3Residual.txt"),  #set to this if running on Data

    fileOutName = cms.string('output.root')
)

