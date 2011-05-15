import FWCore.ParameterSet.Config as cms

process = cms.Process("Step2Lite")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

process.GlobalTag.globaltag = 'RMMEGlobalTag'
#process.GlobalTag.globaltag = 'START311_V2::All'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:RMMEFN'),
#     fileNames = cms.untracked.vstring('file:/data/mangano/MC/Spring11/WW_414_SKIM_V04/ggToH150toWWto2L2Nu.root'),
    inputCommands = cms.untracked.vstring( "keep *" )
)
##############################################

#from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/nov12Skim/WWTo2L2NuZ26noPU/*.root') ]
#process.source.fileNames += [ 'file:%s'%x for x in glob('/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.root') ]

# process.source.fileNames = ['file:hwwSkim.root']

# import commands
# myDir = "/castor/cern.ch/user/m/mwlebour/hww/trees/WW_414_SKIM_V00/mwlebour/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WW_414_SKIM_V00_ID026_WJetsToLNuMad_Xfer/2bd6c1b7661d2d45dc2e9f49eb14f680/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")


##############################################
process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")

process.onlyHiggsGen = cms.EDProducer( "GenParticlePruner",
    src = cms.InputTag("prunedGen"),
    select = cms.vstring(
        "drop  *  ",
        "keep pdgId =   {h0}",
    )
)


# 0
process.wwelmu0 = process.skimEventProducer.clone(hypoType = "WWELMU", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmuel0 = process.skimEventProducer.clone(hypoType = "WWMUEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwelel0 = process.skimEventProducer.clone(hypoType = "WWELEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmumu0 = process.skimEventProducer.clone(hypoType = "WWMUMU", muTag="wwMuMatch", elTag = "wwEleMatch")

# LHT
process.wwelmuIDLHT   = process.wwelmu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwmuelIDLHT   = process.wwmuel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwelelIDLHT   = process.wwelel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwmumuIDLHT   = process.wwmumu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )

process.wwelmuISOLHT  = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwmuelISOLHT  = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwelelISOLHT  = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwmumuISOLHT  = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )

process.wwelmuCONVLHT = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwmuelCONVLHT = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwelelCONVLHT = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwmumuCONVLHT = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )

process.wwelmuIPLHT = process.wwelmu0.clone(muTag = "wwMuonsIP", elTag = "wwEleIPVLHT" )
process.wwmuelIPLHT = process.wwmuel0.clone(muTag = "wwMuonsIP", elTag = "wwEleIPVLHT" )
process.wwelelIPLHT = process.wwelel0.clone(muTag = "wwMuonsIP", elTag = "wwEleIPVLHT" )
process.wwmumuIPLHT = process.wwmumu0.clone(muTag = "wwMuonsIP", elTag = "wwEleIPVLHT" )


process.skimMuMuCONVMET = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumuCONVLHT"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 && pfMet() > 20"),                                   
)
process.skimMuElCONVMET = process.skimMuMuCONVMET.clone( src = "wwmuelCONVLHT" )
process.skimElMuCONVMET = process.skimMuMuCONVMET.clone( src = "wwelmuCONVLHT" )
process.skimElElCONVMET = process.skimMuMuCONVMET.clone( src = "wwelelCONVLHT" )


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RMMEFN'),
#     fileName = cms.untracked.string('step2Lite.root'),
     dropMetaData = cms.untracked.string('ALL'),
     outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_onlyHiggsGen_*_*',
        'keep GenEventInfoProduct_generator__HLT',
        'keep patJets_slimPatJetsTriggerMatch__Yield',
        'keep recoVertexs_offlinePrimaryVertices__Yield',
        'keep PileupSummaryInfos_addPileupInfo__REDIGI311X',
        'keep recoPFMETs_pfMet__REDIGI311X',
        'keep recoMETs_tcMet__REDIGI311X',
        'keep edmTriggerResults_*_*_*',
        'keep recoSkimEvents_ww*CONVLHT_*_*',
        #'keep recoSkimEvents_ww*ISOLHT_*_*',
        #'keep recoSkimEvents_ww*IPLHT_*_*',
        #'keep recoSkimEvents_ww*IPLHT_*_*',
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
)
process.out.SelectEvents.SelectEvents += ['selElMuCONVMET','selMuElCONVMET','selElElCONVMET','selMuMuCONVMET']
#process.out.SelectEvents.SelectEvents += ['selElMuIDMET','selMuElIDMET','selElElIDMET','selMuMuIDMET']
#process.out.SelectEvents.SelectEvents += ['selElMuISOMET','selMuElISOMET','selElElISOMET','selMuMuISOMET']
#process.out.SelectEvents.SelectEvents += ['selElMuIPMET','selMuElIPMET','selElElIPMET','selMuMuIPMET']


process.selMuMuCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 *  process.wwmumu0MET * process.skimMuMu0MET) 
process.selMuElCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 *  process.wwmuel0MET * process.skimMuEl0MET) 
process.selElMuCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 *  process.wwelmu0MET * process.skimElMu0MET) 
process.selElElCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 *  process.wwelel0MET * process.skimElEl0MET) 

process.selMuMuCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 *  process.wwmumuCONVLHT * process.skimMuMuCONVMET) 
process.selMuElCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 *  process.wwmuelCONVLHT * process.skimMuElCONVMET) 
process.selElMuCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 *  process.wwelmuCONVLHT * process.skimElMuCONVMET) 
process.selElElCONVMET = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 *  process.wwelelCONVLHT * process.skimElElCONVMET) 

process.selMuMuIPMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 *  process.wwmumuIPLHT * process.skimMuMuIPMET) 
process.selMuElIPMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 *  process.wwmuelIPLHT * process.skimMuElIPMET) 
process.selElMuIPMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 *  process.wwelmuIPLHT * process.skimElMuIPMET) 
process.selElElIPMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 *  process.wwelelIPLHT * process.skimElElIPMET) 

process.selMuMuIDMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 *  process.wwmumuIDLHT * process.skimMuMuIDMET) 
process.selMuElIDMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 *  process.wwmuelIDLHT * process.skimMuElIDMET) 
process.selElMuIDMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 *  process.wwelmuIDLHT * process.skimElMuIDMET) 
process.selElElIDMET   = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 *  process.wwelelIDLHT * process.skimElElIDMET) 

process.selMuMuISOMET  = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 *  process.wwmumuISOLHT * process.skimMuMuISOMET) 
process.selMuElISOMET  = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 *  process.wwmuelISOLHT * process.skimMuElISOMET) 
process.selElMuISOMET  = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 *  process.wwelmuISOLHT * process.skimElMuISOMET) 
process.selElElISOMET  = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 *  process.wwelelISOLHT * process.skimElElISOMET) 




process.e = cms.EndPath(process.out)

process.genPath = cms.Path(process.onlyHiggsGen)

process.schedule = cms.Schedule(
    process.genPath,
    process.selMuMuCONVMET,process.selMuElCONVMET,process.selElMuCONVMET,process.selElElCONVMET,
    # process.selMuMuIPMET,process.selMuElIPMET,process.selElMuIPMET,process.selElElIPMET,
    # process.selMuMuIDMET,process.selMuElIDMET,process.selElMuIDMET,process.selElElIDMET,
    # process.selMuMuISOMET,process.selMuElISOMET,process.selElMuISOMET,process.selElElISOMET,
    # end
    process.e
)




