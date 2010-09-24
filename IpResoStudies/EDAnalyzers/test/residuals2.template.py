import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
        '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEFC70B6-F53D-DF11-B57E-003048679150.root'
]);



secFiles.extend( [ ]);


process = cms.Process("IpResiduals")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'SET_GLOBALTAG'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = source

#======================================
# HLT
#======================================
process.HLTMinBias = cms.EDFilter("HLTHighLevel",
#                               TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI36"),
                               TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#                               HLTPaths = cms.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus'), 
                               HLTPaths = cms.vstring('HLT_L1Jet6U'), 
                               eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                               andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                               throw = cms.bool(True)    # throw exception on unknown path names
                               )


process.load('IpResoStudies.EDAnalyzers.residuals_cfi')

#### ---------- HERE IS THE SNIPET OF CODE FOR NEW ALIGNMENT/APE --------- ####

#################################################################
##refit tracks with new alignment geometry from Jula Draeger
#################################################################

### Apply corrections to local reco - needed only if you run on TK-DECO data (i.e., MinBias) !!!! 
### For cosmics in PEAK mode, comment them out
process.load("RecoLocalTracker.SiStripRecHitConverter.OutOfTime_cff")
process.OutOfTime.TOBlateBP=0.071
process.OutOfTime.TIBlateBP=0.036

from CondCore.DBCommon.CondDBSetup_cfi import *
from CalibTracker.Configuration.Common.PoolDBESSource_cfi import poolDBESSource
##include private db object
##
import CalibTracker.Configuration.Common.PoolDBESSource_cfi

process.stripLorentzAngle = cms.ESSource("PoolDBESSource",CondDBSetup,
                                        connect = cms.string('sqlite_file:SiStripLorentzAngle_Deco.db'),
                                        toGet = cms.VPSet(cms.PSet(record = cms.string('SiStripLorentzAngleRcd'),
                                                                   tag = cms.string('SiStripLorentzAngle_Deco') ))
                                        )
process.es_prefer_stripLorentzAngle = cms.ESPrefer("PoolDBESSource", "stripLorentzAngle")

### Load the new geometry
process.trackerAlignmentICHEP2010 = cms.ESSource("PoolDBESSource",
                                       CondDBSetup,
                                       toGet = cms.VPSet(cms.PSet(
                                               record = cms.string('TrackerAlignmentRcd'),
                                               tag = cms.string('Alignments')
                                               )),
                                       connect = cms.string('sqlite_file:TOBCenteredObjectICHEP2010.db')
                                       )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignmentICHEP2010")

### Load the new Tracker APE
process.trackerAPEICHEP2010 = cms.ESSource("PoolDBESSource",CondDBSetup,
                          connect = cms.string('sqlite_file:APEforICHEP30umFPIX.db'),
                          toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
                          tag = cms.string('AlignmentErrors') ))
#connect = cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X'),
#toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
#tag = cms.string('TrackerIdealGeometryErrors210_mc') ))
                          )
process.es_prefer_trackerAPE = cms.ESPrefer("PoolDBESSource", "trackerAPEICHEP2010")

# refitting
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.generalTracks = process.TrackRefitter.clone(      src ='generalTracks',#'ALCARECOTkAlCosmicsCTF0T',
   TrajectoryInEvent = True,
   NavigationSchool = "",
   TTRHBuilder = "WithAngleAndTemplate" #default
)



#### --------------------------------------------------------------------- ####

## Primary Vertex
process.offlinePrimaryVerticesWithBS.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minSiliconHits = 6
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minPixelHits = 1
process.offlinePrimaryVerticesWithBS.TkClusParameters.zSeparation = 1
process.offlinePrimaryVertices.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVertices.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVertices.TkFilterParameters.minSiliconHits = 6
process.offlinePrimaryVertices.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVertices.TkFilterParameters.minPixelHits = 1
process.offlinePrimaryVertices.TkClusParameters.zSeparation = 1


process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("SET_OUTPUT"),
      closeFileFast = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

process.p = cms.Path(process.HLTMinBias*
                     process.generalTracks*process.offlinePrimaryVertices*
                     #process.offlineBeamSpot*
                     process.residuals
                     )
