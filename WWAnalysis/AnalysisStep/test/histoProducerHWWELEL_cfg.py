import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.histoProducer_cfi")

process.HistoProducerParams = process.FWLiteParams.clone()
#process.HistoProducerParams.eventHypo = cms.string("wwelel0")
#process.HistoProducerParams.eventHypo = cms.string("wwelelID")
#process.HistoProducerParams.eventHypo = cms.string("wwelelISO")
#process.HistoProducerParams.eventHypo = cms.string("wwelelNOCONV")
process.HistoProducerParams.eventHypo = cms.string("wwelelIP")
process.HistoProducerParams.nMu = cms.int32(0)
process.HistoProducerParams.nEl = cms.int32(2)
process.HistoProducerParams.minDiffMz = cms.double(15.0)
process.HistoProducerParams.minProjMet = cms.double(35.0)

