import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.histoProducer_cfi")

process.HistoProducerParams = process.FWLiteParams.clone()
process.HistoProducerParams.eventHypo = cms.string("wwelmu")
process.HistoProducerParams.nMu = cms.int32(1)
process.HistoProducerParams.nEl = cms.int32(1)
process.HistoProducerParams.minDiffMz = cms.double(-1)
process.HistoProducerParams.minProjMet = cms.double(20.0)
process.HistoProducerParams.inputFolder =     cms.string("/nfs/bluearc/group/trees/ww/SkimPassREPLACE_PASS/")
process.HistoProducerParams.input       =     cms.vstring(REPLACEME_INPUT)
process.HistoProducerParams.fileOutName =     cms.string(REPLACEME_OUTPUT)

