import FWCore.ParameterSet.Config as cms

skimEventProducer = cms.EDProducer('SkimEventProducer',
   muTag         = cms.InputTag("wwMuons"),
   elTag         = cms.InputTag("wwElectrons"),
   softMuTag     =   cms.InputTag("wwMuons4Veto"),
#   extraElTag    = cms.InputTag("wwElectrons"),
   jetTag        = cms.InputTag("slimPatJetsTriggerMatch"), 
   tagJetTag     = cms.InputTag("slimPatJetsTriggerMatch"), 
   pfMetTag      = cms.InputTag("pfMet"),
   tcMetTag      = cms.InputTag("tcMet"),
   chargedMetTag = cms.InputTag("trackMetProducer"),
   vtxTag        = cms.InputTag("goodPrimaryVertices"),
   sptTag        = cms.InputTag("vertexMapProd","sumPt"),
   spt2Tag       = cms.InputTag("vertexMapProd","sumPt2"),
#    branchAlias   = cms.string("wwelmu"),
   hypoType      = cms.string("WWELMU"),

   triggerTag        = cms.InputTag("TriggerResults","","HLT"),
   singleMuDataPaths = cms.vstring("HLT_Mu24_v*"),
   doubleMuDataPaths = cms.vstring("HLT_DoubleMu7_v*"),
   doubleElDataPaths = cms.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*"),
   muEGDataPaths     = cms.vstring("HLT_Mu8_Ele17_CaloIdL_v*","HLT_Mu17_Ele8_CaloIdL_v*"),
   singleMuMCPaths   = cms.vstring("HLT_Mu21_v*"),
   doubleMuMCPaths   = cms.vstring("HLT_DoubleMu5_v*"),
   doubleElMCPaths   = cms.vstring("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*"),
   muEGMCPaths       = cms.vstring("HLT_Mu5_Ele17_v*","HLT_Mu11_Ele8_v*"),
)

