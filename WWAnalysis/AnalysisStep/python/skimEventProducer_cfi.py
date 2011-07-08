import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.wwElectrons_cfi import *
from WWAnalysis.AnalysisStep.wwMuons_cfi import *


skimEventProducer = cms.EDProducer('SkimEventProducer',
    muTag         = cms.InputTag("wwMuons"),
    elTag         = cms.InputTag("wwElectrons"),
    softMuTag     =   cms.InputTag("wwMuons4Veto"),
    # extraElTag    = cms.InputTag("wwElectrons"),
    jetTag        = cms.InputTag("slimPatJetsTriggerMatch"), 
    tagJetTag     = cms.InputTag("slimPatJetsTriggerMatch"), 
    pfMetTag      = cms.InputTag("pfMet"),
    tcMetTag      = cms.InputTag("tcMet"),
    chargedMetTag = cms.InputTag("trackMetProducer"),
    vtxTag        = cms.InputTag("goodPrimaryVertices"),
    chCandsTag    = cms.InputTag("reducedPFCands"),
    sptTag        = cms.InputTag("vertexMapProd","sumPt"),
    spt2Tag       = cms.InputTag("vertexMapProd","sumPt2"),
    # branchAlias   = cms.string("wwelmu"),
    hypoType      = cms.string("WWELMU"),
    
    triggerTag        = cms.InputTag("TriggerResults","","HLT"),
    singleMuDataPaths = cms.vstring(
        "1-163261:HLT_Mu15_v*",
        "163262-164237:HLT_Mu24_v*",
        "165085-166967:HLT_Mu30_v*",
        "163262-166967:HLT_IsoMu17_v*"
        "167039-999999:HLT_IsoMu20_eta2p1_v1*"
    ),
    doubleMuDataPaths = cms.vstring(
        "1-164237:HLT_DoubleMu7_v*",
        "165085-999999:HLT_Mu13_Mu8_v*"
    ),
    doubleElDataPaths = cms.vstring(
        "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*",
        #"HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*"
    ),
    muEGDataPaths     = cms.vstring(
        "HLT_Mu8_Ele17_CaloIdL_v*",
        "HLT_Mu17_Ele8_CaloIdL_v*"
    ),
    singleElDataPaths = cms.vstring("*"
#         "1-164237:HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
#         "165085-999999:HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*"
    ),
    singleMuMCPaths   = cms.vstring("*"),
    singleElMCPaths   = cms.vstring("*"),
    doubleMuMCPaths   = cms.vstring("*"),
    doubleElMCPaths   = cms.vstring("*"),
    muEGMCPaths       = cms.vstring("*"),
    # singleMuMCPaths   = cms.vstring("HLT_Mu21_v*"),
    # singleElMCPaths   = cms.vstring("FILLME"),
    # doubleMuMCPaths   = cms.vstring("HLT_DoubleMu5_v*"),
    # doubleElMCPaths   = cms.vstring("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*"),
    # muEGMCPaths       = cms.vstring("HLT_Mu5_Ele17_v*","HLT_Mu11_Ele8_v*"),
)

def addEventHypothesis(process,label,thisMuTag,thisEleTag,thisSoftMuTag='wwMuons4Veto',peakingType=None):
    hypos = ['mumu','muel','elmu','elel']
    process.peakingFilter = cms.EDFilter("GenFilterDiBosons")

    tempSkimEventFilter = cms.EDFilter("SkimEventSelector",
       src = cms.InputTag(""),
       filter = cms.bool(True),
       cut = cms.string("nLep >=2 "),                                   
    )

    for hypo in hypos:
        #create the four hypothesis:
        setattr(process,'ww'+hypo+label,process.skimEventProducer.clone(hypoType='WW'+hypo.upper(),muTag=thisMuTag,elTag=thisEleTag))
        #create SkimEventSelectors (asking for nLep >=2) 
        setattr(process,'skim'+hypo+label,tempSkimEventFilter.clone(src='ww'+hypo+label))
        # create sequence
        p = cms.Path()
        if peakingType == 'peaking':     p = cms.Path( process.peakingFilter)
        if peakingType == 'non-peaking': p = cms.Path(~process.peakingFilter)
        p += ( 
            getattr(process,thisMuTag)  +
            getattr(process,thisEleTag)  +
            getattr(process,thisSoftMuTag)  +
            getattr(process,'ww'+hypo+label)  +
            getattr(process,'skim'+hypo+label)
        )
        setattr(process,'sel'+hypo+label,p)
        # add to scheduler
        process.schedule.append( getattr(process,'sel'+hypo+label) )
        # add to pooloutput module
        process.out.outputCommands.append( 'keep *_{0}_*_*'.format( 'ww'+hypo+label ) )
        process.out.SelectEvents.SelectEvents.append( 'sel'+hypo+label )


