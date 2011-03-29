import FWCore.ParameterSet.Config as cms

process = cms.Process("WWEvents")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

process.GlobalTag.globaltag = 'START311_V2::All'
# process.GlobalTag.globaltag = 'SET_GLOBALTAG'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
   fileNames = cms.untracked.vstring(
        "file:/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.lowPU.root",
#        "file:/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.highPU.root",
#        "file:/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.veryHighPU.root",
  )
)
##############################################

#from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/nov12Skim/WWTo2L2NuZ26noPU/*.root') ]
#process.source.fileNames += [ 'file:%s'%x for x in glob('/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.root') ]

# process.source.fileNames = ['file:hwwSkim.root']


process.source.inputCommands = cms.untracked.vstring( 
    "keep *", 
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")


### HERE I SET THE SKIMEVENT PRODUCER CONFIGURATION ###
process.wwelmu0 = process.skimEventProducer.clone(hypoType = cms.string("WWELMU"))
process.wwelel0 = process.skimEventProducer.clone(hypoType = cms.string("WWELEL"))
process.wwmumu0 = process.skimEventProducer.clone(hypoType = cms.string("WWMUMU"))
process.wwelmu0.muTag = cms.InputTag("wwMuMatch")
process.wwelel0.muTag = cms.InputTag("wwMuMatch")
process.wwmumu0.muTag = cms.InputTag("wwMuMatch")
process.wwelmu0.elTag = cms.InputTag("wwEleMatch")
process.wwelel0.elTag = cms.InputTag("wwEleMatch")
process.wwmumu0.elTag = cms.InputTag("wwEleMatch")


process.wwelmuID = process.skimEventProducer.clone(hypoType = cms.string("WWELMU"))
process.wwelelID = process.skimEventProducer.clone(hypoType = cms.string("WWELEL"))
process.wwmumuID = process.skimEventProducer.clone(hypoType = cms.string("WWMUMU"))
process.wwelmuID.muTag = cms.InputTag("wwMuID")
process.wwelelID.muTag = cms.InputTag("wwMuID")
process.wwmumuID.muTag = cms.InputTag("wwMuID")
process.wwelmuID.elTag = cms.InputTag("wwEleID")
process.wwelelID.elTag = cms.InputTag("wwEleID")
process.wwmumuID.elTag = cms.InputTag("wwEleID")

process.wwelmuISO = process.skimEventProducer.clone(hypoType = cms.string("WWELMU"))
process.wwelelISO = process.skimEventProducer.clone(hypoType = cms.string("WWELEL"))
process.wwmumuISO = process.skimEventProducer.clone(hypoType = cms.string("WWMUMU"))
process.wwelmuISO.muTag = cms.InputTag("wwMuISO")
process.wwelelISO.muTag = cms.InputTag("wwMuISO")
process.wwmumuISO.muTag = cms.InputTag("wwMuISO")
process.wwelmuISO.elTag = cms.InputTag("wwEleISO")
process.wwelelISO.elTag = cms.InputTag("wwEleISO")
process.wwmumuISO.elTag = cms.InputTag("wwEleISO")

process.wwelmuNOCONV = process.skimEventProducer.clone(hypoType = cms.string("WWELMU"))
process.wwelelNOCONV = process.skimEventProducer.clone(hypoType = cms.string("WWELEL"))
process.wwmumuNOCONV = process.skimEventProducer.clone(hypoType = cms.string("WWMUMU"))
process.wwelmuNOCONV.muTag = cms.InputTag("wwMuISO")
process.wwelelNOCONV.muTag = cms.InputTag("wwMuISO")
process.wwmumuNOCONV.muTag = cms.InputTag("wwMuISO")
process.wwelmuNOCONV.elTag = cms.InputTag("wwEleNOCONV")
process.wwelelNOCONV.elTag = cms.InputTag("wwEleNOCONV")
process.wwmumuNOCONV.elTag = cms.InputTag("wwEleNOCONV")

process.wwelmuIP = process.skimEventProducer.clone(hypoType = cms.string("WWELMU"))
process.wwelelIP = process.skimEventProducer.clone(hypoType = cms.string("WWELEL"))
process.wwmumuIP = process.skimEventProducer.clone(hypoType = cms.string("WWMUMU"))
process.wwelmuIP.muTag = cms.InputTag("wwMuIP")
process.wwelelIP.muTag = cms.InputTag("wwMuIP")
process.wwmumuIP.muTag = cms.InputTag("wwMuIP")
process.wwelmuIP.elTag = cms.InputTag("wwEleIP")
process.wwelelIP.elTag = cms.InputTag("wwEleIP")
process.wwmumuIP.elTag = cms.InputTag("wwEleIP")


##############################################

#0
process.skimElMu0 = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmu0"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimElEl0 = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelel0"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimMuMu0 = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumu0"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)


#ID
process.skimElMuID = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmuID"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimElElID = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelelID"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimMuMuID = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumuID"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)

#ISO
process.skimElMuISO = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmuISO"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimElElISO = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelelISO"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimMuMuISO = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumuISO"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)

#NOCONV
process.skimElMuNOCONV = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmuNOCONV"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimElElNOCONV = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelelNOCONV"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimMuMuNOCONV = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumuNOCONV"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)

#IP
process.skimElMuIP = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmuIP"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimElElIP = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelelIP"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimMuMuIP = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumuIP"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)


process.out0 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.0.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_slimPatJetsTriggerMatch_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_wwelmu0_*_*',
        'keep *_wwelel0_*_*',
        'keep *_wwmumu0_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMu0','selElEl0','selMuMu0' ))
)


process.outID = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.ID.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_slimPatJetsTriggerMatch_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_wwelmuID_*_*',
        'keep *_wwelelID_*_*',
        'keep *_wwmumuID_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMuID','selElElID','selMuMuID' ))
)

process.outISO = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.ISO.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_slimPatJetsTriggerMatch_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_wwelmuISO_*_*',
        'keep *_wwelelISO_*_*',
        'keep *_wwmumuISO_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMuISO','selElElISO','selMuMuISO' ))
)

process.outNOCONV = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.NOCONV.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_slimPatJetsTriggerMatch_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_wwelmuNOCONV_*_*',
        'keep *_wwelelNOCONV_*_*',
        'keep *_wwmumuNOCONV_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMuNOCONV','selElElNOCONV','selMuMuNOCONV' ))
)

process.outIP = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.IP.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_slimPatJetsTriggerMatch_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_wwelmuIP_*_*',
        'keep *_wwelelIP_*_*',
        'keep *_wwmumuIP_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMuIP','selElElIP','selMuMuIP' ))
)

  
process.p = cms.Path( 
    process.wwElectronSequence *
    process.wwMuonSequence *
    process.wwelmu0*process.wwelel0*process.wwmumu0*
    process.wwelmuID*process.wwelelID*process.wwmumuID*
    process.wwelmuISO*process.wwelelISO*process.wwmumuISO*
    process.wwelmuNOCONV*process.wwelelNOCONV*process.wwmumuNOCONV*
    process.wwelmuIP*process.wwelelIP*process.wwmumuIP
)

process.selElMu0 = cms.Path(process.skimElMu0)
process.selElEl0 = cms.Path(process.skimElEl0)
process.selMuMu0 = cms.Path(process.skimMuMu0)

process.selElMuID = cms.Path(process.skimElMuID)
process.selElElID = cms.Path(process.skimElElID)
process.selMuMuID = cms.Path(process.skimMuMuID)

process.selElMuISO = cms.Path(process.skimElMuISO)
process.selElElISO = cms.Path(process.skimElElISO)
process.selMuMuISO = cms.Path(process.skimMuMuISO)

process.selElMuNOCONV = cms.Path(process.skimElMuNOCONV)
process.selElElNOCONV = cms.Path(process.skimElElNOCONV)
process.selMuMuNOCONV = cms.Path(process.skimMuMuNOCONV)

process.selElMuIP = cms.Path(process.skimElMuIP)
process.selElElIP = cms.Path(process.skimElElIP)
process.selMuMuIP = cms.Path(process.skimMuMuIP)



process.e = cms.EndPath(process.out0*process.outID*process.outISO*process.outNOCONV*process.outIP)

process.sched = cms.Schedule(process.p,
                             process.selElMu0,process.selElEl0,process.selMuMu0,
                             process.selElMuID,process.selElElID,process.selMuMuID,
                             process.selElMuISO,process.selElElISO,process.selMuMuISO,
                             process.selElMuNOCONV,process.selElElNOCONV,process.selMuMuNOCONV,
                             process.selElMuIP,process.selElElIP,process.selMuMuIP,
                             process.e
)



