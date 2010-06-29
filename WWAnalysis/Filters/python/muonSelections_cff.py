import FWCore.ParameterSet.Config as cms

# --- selections of muon plus
from WWAnalysis.Filters.muPlusSelection1_cfi import *

# --- selections of muon minus
import WWAnalysis.Filters.muPlusSelection1_cfi
muMinusSelection1 = WWAnalysis.Filters.muPlusSelection1_cfi.muPlusSelection1.clone(
    charge = cms.int32(-1)
)



