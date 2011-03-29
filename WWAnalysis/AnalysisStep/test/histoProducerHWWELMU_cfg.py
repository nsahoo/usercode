import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.histoProducer_cfi")

process.HistoProducerParams = process.FWLiteParams.clone()
#process.HistoProducerParams.eventHypo = cms.string("wwelmu0")
#process.HistoProducerParams.eventHypo = cms.string("wwelmuID")
#process.HistoProducerParams.eventHypo = cms.string("wwelmuISO")
process.HistoProducerParams.eventHypo = cms.string("wwelmuNOCONV")
#process.HistoProducerParams.eventHypo = cms.string("wwelmuIP")
process.HistoProducerParams.nMu = cms.int32(1)
process.HistoProducerParams.nEl = cms.int32(1)
process.HistoProducerParams.minDiffMz = cms.double(-1)
process.HistoProducerParams.minProjMet = cms.double(20.0)

