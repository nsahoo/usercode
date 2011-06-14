import FWCore.ParameterSet.Config as cms

        
from PhysicsTools.PatAlgos.patSequences_cff import *

selectedRefPatElectrons = cms.EDFilter("PATElectronViewRefSelector",
   src = cms.InputTag("input"),
   cut = cms.string("")
)


wwEleMatch = selectedRefPatElectrons.clone()
wwEleMatch.src = "boostedElectrons"
wwEleMatch.filter = cms.bool(False)
wwEleMatch.cut = ( "pt > 10 && abs(eta)<2.5 ")
#wwEleMatch.cut = ( "pt > 10 && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==11 && abs(genParticleRef.get().mother().mother().pdgId()) ==24")


from WWAnalysis.AnalysisStep.electronIDs_cff import *


# LHL
wwEleIDLHL = selectedRefPatElectrons.clone()
wwEleIDLHL.src = "wwEleMatch"
wwEleIDLHL.filter = cms.bool(False)
wwEleIDLHL.cut = ( ELE_ID_LH_95_2011 )

wwEleISOLHL = selectedRefPatElectrons.clone()
wwEleISOLHL.src = "wwEleIDLHL"
wwEleISOLHL.filter = cms.bool(False)
wwEleISOLHL.cut = ( ELE_ISO_LH_95_2011 )

wwEleCONVLHL = selectedRefPatElectrons.clone()
wwEleCONVLHL.src = "wwEleISOLHL"
wwEleCONVLHL.filter = cms.bool(False)
wwEleCONVLHL.cut = ( ELE_NOCONV )

wwEleIPLHL = selectedRefPatElectrons.clone()
wwEleIPLHL.src = "wwEleCONVLHL"
wwEleIPLHL.filter = cms.bool(False)
wwEleIPLHL.cut = ( ELE_IP )

# LHT
wwEleIDLHT = selectedRefPatElectrons.clone()
wwEleIDLHT.src = "wwEleMatch"
wwEleIDLHT.filter = cms.bool(False)
wwEleIDLHT.cut = ( ELE_ID_LH_90_2011 )

wwEleISOLHT = selectedRefPatElectrons.clone()
wwEleISOLHT.src = "wwEleIDLHT"
wwEleISOLHT.filter = cms.bool(False)
wwEleISOLHT.cut = ( ELE_ISO_LH_90_2011 )

wwEleCONVLHT = selectedRefPatElectrons.clone()
wwEleCONVLHT.src = "wwEleISOLHT"
wwEleCONVLHT.filter = cms.bool(False)
wwEleCONVLHT.cut = ( ELE_NOCONV )

wwEleIPLHT = selectedRefPatElectrons.clone()
wwEleIPLHT.src = "wwEleCONVLHT"
wwEleIPLHT.filter = cms.bool(False)
wwEleIPLHT.cut = ( ELE_IP )

wwEleISOPFLHT = selectedRefPatElectrons.clone()
wwEleISOPFLHT.src = "wwEleIDLHT"
wwEleISOPFLHT.filter = cms.bool(False)
wwEleISOPFLHT.cut = ( ELE_ISOPF_LH_90_2011 )

wwEleCONVPFLHT = selectedRefPatElectrons.clone()
wwEleCONVPFLHT.src = "wwEleISOPFLHT"
wwEleCONVPFLHT.filter = cms.bool(False)
wwEleCONVPFLHT.cut = ( ELE_NOCONV )

wwEleIPPFLHT = selectedRefPatElectrons.clone()
wwEleIPPFLHT.src = "wwEleCONVPFLHT"
wwEleIPPFLHT.filter = cms.bool(False)
wwEleIPPFLHT.cut = ( ELE_IP )

# CBL
wwEleIDCBL = selectedRefPatElectrons.clone()
wwEleIDCBL.src = "wwEleMatch"
wwEleIDCBL.filter = cms.bool(False)
wwEleIDCBL.cut = ( ELE_ID_CB_95_2011 )

wwEleISOCBL = selectedRefPatElectrons.clone()
wwEleISOCBL.src = "wwEleIDCBL"
wwEleISOCBL.filter = cms.bool(False)
wwEleISOCBL.cut = ( ELE_ISO_CB_95_2011 )

wwEleCONVCBL = selectedRefPatElectrons.clone()
wwEleCONVCBL.src = "wwEleISOCBL"
wwEleCONVCBL.filter = cms.bool(False)
wwEleCONVCBL.cut = ( ELE_NOCONV )

wwEleIPCBL = selectedRefPatElectrons.clone()
wwEleIPCBL.src = "wwEleCONVCBL"
wwEleIPCBL.filter = cms.bool(False)
wwEleIPCBL.cut = ( ELE_IP )

# CBT
wwEleIDCBT = selectedRefPatElectrons.clone()
wwEleIDCBT.src = "wwEleMatch"
wwEleIDCBT.filter = cms.bool(False)
wwEleIDCBT.cut = ( ELE_ID_CB_90_2011 )

wwEleISOCBT = selectedRefPatElectrons.clone()
wwEleISOCBT.src = "wwEleIDCBT"
wwEleISOCBT.filter = cms.bool(False)
wwEleISOCBT.cut = ( ELE_ISO_CB_90_2011 )

wwEleCONVCBT = selectedRefPatElectrons.clone()
wwEleCONVCBT.src = "wwEleISOCBT"
wwEleCONVCBT.filter = cms.bool(False)
wwEleCONVCBT.cut = ( ELE_NOCONV )

wwEleIPCBT = selectedRefPatElectrons.clone()
wwEleIPCBT.src = "wwEleCONVCBT"
wwEleIPCBT.filter = cms.bool(False)
wwEleIPCBT.cut = ( ELE_IP )

# Merge
wwEleIDMerge = selectedRefPatElectrons.clone()
wwEleIDMerge.src = "wwEleMatch"
wwEleIDMerge.filter = cms.bool(False)
wwEleIDMerge.cut = ( ELE_MERGE_ID )

wwEleISOMerge = selectedRefPatElectrons.clone()
wwEleISOMerge.src = "wwEleIDMerge"
wwEleISOMerge.filter = cms.bool(False)
wwEleISOMerge.cut = ( ELE_MERGE_ISO )

wwEleCONVMerge = selectedRefPatElectrons.clone()
wwEleCONVMerge.src = "wwEleISOMerge"
wwEleCONVMerge.filter = cms.bool(False)
wwEleCONVMerge.cut = ( ELE_MERGE_CONV )

wwEleIPMerge = selectedRefPatElectrons.clone()
wwEleIPMerge.src = "wwEleCONVMerge"
wwEleIPMerge.filter = cms.bool(False)
wwEleIPMerge.cut = ( ELE_MERGE_IP )

# Merge2
wwEleIDMerge2 = selectedRefPatElectrons.clone()
wwEleIDMerge2.src = "wwEleMatch"
wwEleIDMerge2.filter = cms.bool(False)
wwEleIDMerge2.cut = ( ELE_MERGE_ID2 )

wwEleISOMerge2 = selectedRefPatElectrons.clone()
wwEleISOMerge2.src = "wwEleIDMerge2"
wwEleISOMerge2.filter = cms.bool(False)
wwEleISOMerge2.cut = ( ELE_MERGE_ISO )

wwEleCONVMerge2 = selectedRefPatElectrons.clone()
wwEleCONVMerge2.src = "wwEleISOMerge2"
wwEleCONVMerge2.filter = cms.bool(False)
wwEleCONVMerge2.cut = ( ELE_MERGE_CONV )

wwEleIPMerge2 = selectedRefPatElectrons.clone()
wwEleIPMerge2.src = "wwEleCONVMerge2"
wwEleIPMerge2.filter = cms.bool(False)
wwEleIPMerge2.cut = ( ELE_MERGE_IP )


wwElectronSequence = cms.Sequence(  
    wwEleMatch *
#     wwEleIDLHL *
#     wwEleISOLHL *
#     wwEleCONVLHL *
#     wwEleIPLHL *
#    wwEleIDLHT *
#    wwEleISOLHT *
#    wwEleCONVLHT *
#    wwEleIPLHT *
    wwEleIDMerge *
    wwEleISOMerge *
    wwEleCONVMerge *
    wwEleIPMerge *
    wwEleIDMerge2 *
    wwEleISOMerge2 *
    wwEleCONVMerge2 *
    wwEleIPMerge2 
#     wwEleISOPFLHT *
#     wwEleCONVPFLHT *
#     wwEleIPPFLHT 
#     wwEleIDCBL *
#     wwEleISOCBL *
#     wwEleCONVCBL *
#     wwEleIPCBL *
#     wwEleIDCBT *
#     wwEleISOCBT *
#     wwEleCONVCBT *
#     wwEleIPCBT 
)


