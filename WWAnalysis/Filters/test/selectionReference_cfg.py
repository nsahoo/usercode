import FWCore.ParameterSet.Config as cms

process = cms.Process("WWSel2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START36_V9::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_10_1_u74.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_11_1_lXF.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_12_1_ozc.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_13_1_Pmh.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_1_1_u0x.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_2_1_KQI.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_3_1_fTy.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_4_1_XjO.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_5_1_j8A.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_6_1_xKg.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_7_1_IKi.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_8_1_33L.root",
        "file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_9_1_Rxw.root"
    ),
)


#======================================
# HLT
#======================================
process.selectionHLT = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI"),
        HLTPaths = cms.vstring('HLT_Mu9','HLT_Ele15_LW_L1R'), 
        eventSetupPathsKey = cms.string(''), 
        andOr = cms.bool(True),            
        throw = cms.bool(True),
)

process.load('WWAnalysis.Filters.electronSelections_cff')
process.load('WWAnalysis.Filters.muonSelections_cff')
process.load('WWAnalysis.Filters.eleMuSelections_cff')

process.load('WWAnalysis.Filters.daugtherListCleaner_cfi')
process.daugtherListCleanerDiMuon = process.daugtherListCleaner.clone(src = cms.InputTag("diMuonSel2"))
process.daugtherListCleanerDiEle  = process.daugtherListCleaner.clone(src = cms.InputTag("diEleSel2"))
process.daugtherListCleanerEleMu  = process.daugtherListCleaner.clone(src = cms.InputTag("diLeptSel2"))

process.load('WWAnalysis.Filters.metFilter_cfi')
process.metFilterDiMuon = process.metFilter.clone(srcLeptons = "daugtherListCleanerDiMuon")
process.metFilterDiEle  = process.metFilter.clone(srcLeptons = "daugtherListCleanerDiEle" )
process.metFilterEleMu  = process.metFilter.clone(srcLeptons = "daugtherListCleanerEleMu", skipMEtCut = True )
process.load('WWAnalysis.Filters.jetVetoFilter_cfi')
process.load('WWAnalysis.Filters.softMuonVeto_cfi')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

process.DiMuSequence = cms.Sequence(
    process.muPlusSelection1*process.muMinusSelection1*        
    process.muPlusIPSelection*process.muMinusIPSelection*    
    process.muPlusIsoSelection*process.muMinusIsoSelection*  
    process.muPlusID*process.muMinusID*                      
    process.diMuonSel1*process.diMuonFilter1*                
    process.diMuonSel2*                                      
    process.daugtherListCleanerDiMuon
)
process.DiEleSequence = cms.Sequence(
    process.elePlusSelection1*process.eleMinusSelection1*       
    process.elePlusIPSelection*process.eleMinusIPSelection*   
    process.elePlusIsoSelection*process.eleMinusIsoSelection* 
    process.elePlusID*process.eleMinusID*                     
    process.convRejectionPlus*process.convRejectionMinus* 
    process.diEleSel1*process.diEleFilter1*                   
    process.diEleSel2*                                        
    process.daugtherListCleanerDiEle                             
)
process.EleMuSequence = cms.Sequence(
    process.eleSelection1*process.muSelection1*
    process.selection1Pairs*process.selection1PairsCounter* 
    process.eleIPSelection*process.muIPSelection*
    process.ipSelPairs*process.ipSelPairsCounter*           
    process.eleIsoSelection*process.muIsoSelection*
    process.isoSelPairs*process.isoSelPairsCounter*         
    process.eleID*process.muID*
    process.idSelPairs*process.idSelPairsCounter*           
    process.eleConvRejection*
    process.convSelPairs*process.convSelPairsCounter*       
    process.diLeptSel1*process.diLeptFilter1*               
    process.diLeptSel2*                                     
    process.daugtherListCleanerEleMu
)
process.jetVetoSequence = cms.Sequence(
    process.eleForVetoSequence*
    process.muForVetoSequence*
    process.jetVetoFilter                                   
)


process.DiMuon_V0 = cms.Path(
    process.selectionHLT*
    process.DiMuSequence*
    process.metSel1*                                          
    process.metFilterDiMuon*                                        
    process.jetVetoSequence*
    process.softMuonVeto*                                     
    process.extraLeptonVetoForMu                             
)
process.DiEle_V0 = cms.Path(
    process.selectionHLT*
    process.DiEleSequence*
    process.metSel1*                                         
    process.metFilterDiEle* 
    process.jetVetoSequence*
    process.eleSoftMuonVeto*                                  
    process.extraLeptonVetoForEle                             
)

process.EleMu_V0 = cms.Path(
    process.selectionHLT*
    process.EleMuSequence*
    process.metFilterEleMu*                                      
    process.jetVetoSequence*
    process.leptSoftMuonVeto*                               
    process.extraLeptonVetoForLept                          
)

# ---- endPath ----
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    outputCommands = cms.untracked.vstring('drop *',
                                           'keep *_diLeptSel1_*_WWSel2',
                                           'keep *_diLeptSel2_*_WWSel2'
                                           ),
#    outputCommands = process.AODEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p2"))
)
#process.end = cms.EndPath(process.out)
process.outDebug = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("selected_EleMu.root"),
    outputCommands = cms.untracked.vstring('keep *', 'drop *_*_*_WWSel2', 'keep *_ak5PFJets_*_*','keep *_daugtherListCleaner_*_*', 'keep *_TriggerResults_*_*'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p2"))#,"p2nv"))
)
#process.endDebug = cms.EndPath(process.outDebug)

if False:
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.p2, "tcMet", "pfMet")
    process.p2.replace(process.selectionHLT, process.ak5PFJets + process.selectionHLT)
    from glob import glob
    process.source.fileNames = [ 'file:%s' % s for s in glob('/nfs/bluearc/group/skims/ww/TTbarJets-madgraph-v2/*.root') ]
    #process.source.fileNames = [ 'file:%s' % s for s in glob('/nfs/bluearc/group/skims/ww/WJets-madgraph-v2/*.root') ]
    process.outDebug.fileName = "selected_EleMu_ttbar.root"

if False:
    #from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    #massSearchReplaceAnyInputTag(process.p2, "tcMet", "pfMet")
    #process.p2.replace(process.selectionHLT, process.ak5PFJets + process.selectionHLT)
    #from glob import glob
    process.source.fileNames = [ 'file:pf2pat_debug_ttbar.root' ]
    #process.source.fileNames = [ 'file:%s' % s for s in glob('/nfs/bluearc/group/skims/ww/WJets-madgraph-v2/*.root') ]
    process.outDebug.fileName = "selected_EleMu_ttbar_onpf2pat.root"
    process.outDebug.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p2","p2nv"))
