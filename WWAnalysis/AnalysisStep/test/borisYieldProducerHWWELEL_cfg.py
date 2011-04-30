import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.borisYieldProducer_cfi")


process.HistoProducerParams = process.FWLiteParams.clone()
process.HistoProducerParams.eventHypo = cms.string("label")
process.HistoProducerParams.nMu = cms.int32(0)
process.HistoProducerParams.nEl = cms.int32(2)
process.HistoProducerParams.minDiffMz = cms.double(15.0)
process.HistoProducerParams.minProjMet = cms.double(35.0)
process.HistoProducerParams.fileOutName =     cms.string("")
process.HistoProducerParams.eventHypo0      = cms.string("label0")
process.HistoProducerParams.eventHypoID     = cms.string("labelID")
process.HistoProducerParams.eventHypoISO    = cms.string("labelISO")
process.HistoProducerParams.eventHypoNOCONV = cms.string("labelNOCONV")
process.HistoProducerParams.eventHypoIP     = cms.string("labelIP")

isOnCastor = replace_isOnCastor
if isOnCastor:
    import commands
    myDir="replace_inputFolder/"
    process.HistoProducerParams.input = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ]
else:
    from glob import glob
    process.HistoProducerParams.input += [ '%s'%x for x in glob('replace_inputFolder/*.root') ]

