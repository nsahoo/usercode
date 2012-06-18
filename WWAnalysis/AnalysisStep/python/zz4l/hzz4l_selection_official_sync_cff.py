import FWCore.ParameterSet.Config as cms

## Include MELA cut

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import SEL_ZZ4L_STEP_6
SEL_ZZ4L_STEP_6 = "userFloat('melaSMH')/(userFloat('melaSMH')+userFloat('melaQQZZ')) > 0.1"
