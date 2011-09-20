import FWCore.ParameterSet.Config as cms

        
from PhysicsTools.PatAlgos.patSequences_cff import *

selectedRefPatElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("boostedElectrons"),
   filter = cms.bool(False),
   cut = cms.string(""),
)


ELE_BASE = "( pt > 10 && abs(eta)<2.5 )"
wwEleMatch = selectedRefPatElectrons.clone( cut = ELE_BASE )

from WWAnalysis.AnalysisStep.electronIDs_cff import *

# Loose
wwEleLoose      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LOOSE_2011 )


# LHL
wwEleIDLHL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 )
wwEleISOLHL     = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 + " && " + ELE_ISO_LH_95_2011 ) 
wwEleCONVLHL    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 + " && " + ELE_ISO_LH_95_2011  + " && " + ELE_NOCONV ) 
wwEleIPLHL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 + " && " + ELE_ISO_LH_95_2011  + " && " + ELE_NOCONV + " && " + ELE_IP ) 

# LHT
wwEleIDLHT      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 ) 
wwEleISOLHT     = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISO_LH_90_2011 ) 
wwEleCONVLHT    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISO_LH_90_2011 + " && " + ELE_NOCONV ) 
wwEleIPLHT      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISO_LH_90_2011 + " && " + ELE_NOCONV + " && " + ELE_IP ) 

wwEleISOPFLHT   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISOPF_LH_90_2011 ) 
wwEleCONVPFLHT  = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISOPF_LH_90_2011 + " && " + ELE_NOCONV ) 
wwEleIPPFLHT    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISOPF_LH_90_2011 + " && " + ELE_NOCONV + " && " + ELE_IP )

# CBL
wwEleIDCBL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 ) 
wwEleISOCBL     = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 + " && " + ELE_ISO_CB_95_2011 ) 
wwEleCONVCBL    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 + " && " + ELE_ISO_CB_95_2011 + " && " + ELE_NOCONV ) 
wwEleIPCBL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 + " && " + ELE_ISO_CB_95_2011 + " && " + ELE_NOCONV + " && " + ELE_IP )

# CBT
wwEleIDCBT      = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 ) 
wwEleISOCBT     = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 + " && " + ELE_ISO_CB_90_2011 ) 
wwEleCONVCBT    = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 + " && " + ELE_ISO_CB_90_2011 + " && " + ELE_NOCONV ) 
wwEleIPCBT      = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 + " && " + ELE_ISO_CB_90_2011 + " && " + ELE_NOCONV + " && " + ELE_IP )

# Merge
wwEleIDMerge    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID ) 
wwEleISOMerge   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO ) 
wwEleCONVMerge  = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV ) 
wwEleIPMerge    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP )

# Merge2
wwEleIDMerge2   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 ) 
wwEleISOMerge2  = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 + " && " + ELE_MERGE_ISO ) 
wwEleCONVMerge2 = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV ) 
wwEleIPMerge2   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP ) 

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


