import FWCore.ParameterSet.Config as cms

        
from PhysicsTools.PatAlgos.patSequences_cff import *

wwEleMatch = selectedPatElectrons.clone()
wwEleMatch.src = "boostedElectrons"
wwEleMatch.filter = cms.bool(False)
wwEleMatch.cut = ( "pt > 10 && abs(eta)<2.5 ")
#wwEleMatch.cut = ( "pt > 10 && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==11 && abs(genParticleRef.get().mother().mother().pdgId()) ==24")


from WWAnalysis.AnalysisStep.electronIDs_cff import *


wwEleIDLHL = selectedPatElectrons.clone()
wwEleIDLHL.src = "wwEleMatch"
wwEleIDLHL.filter = cms.bool(False)
wwEleIDLHL.cut = ( ELE_NOCONV + "&&" + ELE_IP + "&&" + ELE_ID_LH_90_2011 )

wwEleISOLHL = selectedPatElectrons.clone()
wwEleISOLHL.src = "wwEleIDLHL"
wwEleISOLHL.filter = cms.bool(False)
wwEleISOLHL.cut = ( ELE_ISO_LH_90_2011 )

wwEleIDLHT = selectedPatElectrons.clone()
wwEleIDLHT.src = "wwEleMatch"
wwEleIDLHT.filter = cms.bool(False)
wwEleIDLHT.cut = ( ELE_NOCONV + "&&" + ELE_IP + "&&" + ELE_ID_LH_85_2011 )

wwEleISOLHT = selectedPatElectrons.clone()
wwEleISOLHT.src = "wwEleIDLHT"
wwEleISOLHT.filter = cms.bool(False)
wwEleISOLHT.cut = ( ELE_ISO_LH_85_2011 )

wwEleIDCBL = selectedPatElectrons.clone()
wwEleIDCBL.src = "wwEleMatch"
wwEleIDCBL.filter = cms.bool(False)
wwEleIDCBL.cut = ( ELE_NOCONV + "&&" + ELE_IP + "&&" + ELE_ID_CB_90_2011 )

wwEleISOCBL = selectedPatElectrons.clone()
wwEleISOCBL.src = "wwEleIDCBL"
wwEleISOCBL.filter = cms.bool(False)
wwEleISOCBL.cut = ( ELE_ISO_CB_90_2011 )

wwEleIDCBT = selectedPatElectrons.clone()
wwEleIDCBT.src = "wwEleMatch"
wwEleIDCBT.filter = cms.bool(False)
wwEleIDCBT.cut = ( ELE_NOCONV + "&&" + ELE_IP + "&&" + ELE_ID_CB_85_2011 )

wwEleISOCBT = selectedPatElectrons.clone()
wwEleISOCBT.src = "wwEleIDCBT"
wwEleISOCBT.filter = cms.bool(False)
wwEleISOCBT.cut = ( ELE_ISO_CB_85_2011 )

wwElectronSequence = cms.Sequence(  
    wwEleMatch *
    wwEleIDLHL *
    wwEleISOLHL *
    wwEleIDLHT *
    wwEleISOLHT *
    wwEleIDCBL *
    wwEleISOCBL *
    wwEleIDCBT *
    wwEleISOCBT 
)


