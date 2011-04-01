import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.cutPSets_cfi import defaultWW

#default yield params setup for wwelel
FWLiteParams = cms.PSet(
    selectionParams = cms.PSet(
        wwelelISO = cms.PSet(defaultWW)
    ),

    inputParams = cms.PSet(
    ),

#     l2File      = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2Relative.txt"),
#     l3File      = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L3Absolute.txt"),
#     resFile     = cms.string(""),  #set to nothing if running on MC
#     resFile     = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2L3Residual.txt"),  #set to this if running on Data

    fileOutName = cms.string('output.root'),

    histParams = cms.PSet(
#         mll = cms.PSet(
#             variable = cms.string("mll"),
#             nbins = cms.uint32(300),
#             low = cms.double(0),
#             high = cms.double(300),
#             xtitle = cms.string("m_{ll} [GeV]")
#         ),
#         ptMax = cms.PSet(
#             variable = cms.string("ptMax"),
#             nbins = cms.uint32(200),
#             low = cms.double(0),
#             high = cms.double(200),
#             xtitle = cms.string("p_{T} [GeV]")
#         ),
#         ptMin = cms.PSet(
#             variable = cms.string("ptMin"),
#             nbins = cms.uint32(200),
#             low = cms.double(0),
#             high = cms.double(200),
#             xtitle = cms.string("p_{T} [GeV]")
#         ),
#         dPhill = cms.PSet(
#             variable = cms.string("dPhillInDegrees"),
#             nbins = cms.uint32(180),
#             low = cms.double(0),
#             high = cms.double(180),
#             xtitle = cms.string("#Delta#phi_{ll} [deg]")
#         ),
#         mTHiggs = cms.PSet(
#             variable = cms.string("mTHiggs"),
#             nbins = cms.uint32(300),
#             low = cms.double(0),
#             high = cms.double(300),
#             xtitle = cms.string("m_{T} [GeV]")
#         ),
    )

)





