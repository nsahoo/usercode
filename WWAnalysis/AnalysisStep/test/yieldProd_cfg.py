import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.histoProducer_cfi")

process.HistoProducerParams = process.FWLiteParams.clone()
process.HistoProducerParams.eventHypo = cms.string("wwmumu")
# process.HistoProducerParams.eventHypo = cms.string("wwelel")
# process.HistoProducerParams.eventHypo = cms.string("wwelmu")
process.HistoProducerParams.nMu = cms.int32(2)
process.HistoProducerParams.nEl = cms.int32(0)
process.HistoProducerParams.mZ = cms.double(15.0)
process.HistoProducerParams.pMet = cms.double(35.0)
# process.HistoProducerParams.mZ = cms.double(-1.0)        #for elmu
# process.HistoProducerParams.pMet = cms.double(20.0)      #for elmu
process.HistoProducerParams.inputFolder =     cms.string("data/")
process.HistoProducerParams.input       =     cms.vstring('fWWAllEvents') #note: no .root at the end, don't ask me, ask boris

