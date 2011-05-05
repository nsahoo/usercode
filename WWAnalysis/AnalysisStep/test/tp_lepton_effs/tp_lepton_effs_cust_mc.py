# global tag
process.GlobalTag.globaltag = cms.string('START311_V2::All')

# triggers
# 1. tags
process.tagElectrons.cut = process.tagElectrons.cut.value().replace("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT","HLT_Ele22_SW_TighterCaloIdIsol_L1R")
process.tagMuons.cut = process.tagMuons.cut.value().replace("HLT_Mu24","HLT_Mu25")
# 2. probes
from WWAnalysis.AnalysisStep.tp_lepton_effs_cff import EleTriggersMC, MuTriggersMC
for pn, pv in EleTriggersMC.parameters_().items(): 
    setattr(process.tpTreeElEl.flags,    pn, pv)
    setattr(process.tpTreeElEl.tagFlags, pn, pv)
for pn, pv in  MuTriggersMC.parameters_().items(): 
    setattr(process.tpTreeMuMu.flags, pn, pv)

# leave some lines blank at the end please


