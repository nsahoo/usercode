import FWCore.ParameterSet.Config as cms

from WWAnalysis.SkimStep.preYieldFilter_cff import addPreYieldFilter
from WWAnalysis.SkimStep.preFakeFilter_cff import addPreFakeFilter

def addPreFilters(process,addGenFilter):
    addPreYieldFilter(process,addGenFilter)
    addPreFakeFilter(process)


