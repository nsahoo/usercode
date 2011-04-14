import FWCore.ParameterSet.Config as cms

        
from PhysicsTools.PatAlgos.patSequences_cff import *

wwEleMatch = selectedPatElectrons.clone()
wwEleMatch.src = "boostedElectrons"
wwEleMatch.filter = cms.bool(False)
wwEleMatch.cut = ( "pt > 10 && abs(eta)<2.5 ")
#wwEleMatch.cut = ( "pt > 10 && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==11 && abs(genParticleRef.get().mother().mother().pdgId()) ==24")


from WWAnalysis.AnalysisStep.electronIDs_cff import *


wwEleID = selectedPatElectrons.clone()
wwEleID.src = "wwEleMatch"
wwEleID.filter = cms.bool(False)
wwEleID.cut = ( ELE_ID_CB_85_2011 )

wwEleISO = selectedPatElectrons.clone()
wwEleISO.src = "wwEleID"
wwEleISO.filter = cms.bool(False)
wwEleISO.cut = ( ELE_ISO_CB_85_2011 )

wwEleNOCONV = selectedPatElectrons.clone()
wwEleNOCONV.src = "wwEleISO"
wwEleNOCONV.filter = cms.bool(False)
wwEleNOCONV.cut = ( ELE_NOCONV )

wwEleIP = selectedPatElectrons.clone()
wwEleIP.src = "wwEleNOCONV"
wwEleIP.filter = cms.bool(False)
wwEleIP.cut = ( ELE_IP )



wwElectronSequence = cms.Sequence(  
    wwEleMatch*wwEleID*wwEleISO*wwEleNOCONV*wwEleIP
)


