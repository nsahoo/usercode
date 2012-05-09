import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *

# make sure it's  5/3
EL_PT_MIN=7; EL_PT_MIN_LOOSE=7
MU_PT_MIN=5; MU_PT_MIN_LOOSE=5

MUID_LOOSE = "pt > %f && %s" % (MU_PT_MIN_LOOSE, MUID_LOOSE_NO_PT_CUT)
ELID_LOOSE = "pt > %f && %s" % (EL_PT_MIN_LOOSE, ELID_LOOSE_NO_PT_CUT)

MUID_GOOD = " && ".join(["pt > %f" % MU_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])
ELID_GOOD = " && ".join(["pt > %f" % EL_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])

MUID_LOOSE_CR = " && ".join(["pt > %f" % MU_PT_MIN, MU_PRESELECTION, SINGLE_SIP_CUT])
ELID_LOOSE_CR = " && ".join(["pt > %f" % EL_PT_MIN, EL_PRESELECTION, SINGLE_SIP_CUT])

boostedElectronsEAPFIso.separatePhotonAndHadronEAs = False
boostedElectronsEAPFIso.truncateAtZero = "yes"

ARBITRATE_EARLY = True
FOUR_LEPTON_FILTER_PRE_Z = True # plug a 4-lepton count filter before making Z's

SEL_ZZ4L_STEP_1 = "4 < mz(1) < 120"
SEL_ZZ4L_STEP_2 = "lByPt(0).pt > 20 && lByPt(1).pt > 10"
SEL_ZZ4L_STEP_3 = "nGoodPairs(\"mass > 4\", 1) >= 6"
SEL_ZZ4L_STEP_4 = "mass > 70"
SEL_ZZ4L_STEP_5 = "mass > 100"

TRIGGER_FILTER = 'triggerFilter7TeV'

## We modify also the skim sequence
from WWAnalysis.AnalysisStep.zz4l.reSkim_cff import *
muons4skim.cut = "(isTrackerMuon||isGlobalMuon) && abs(eta) < 2.4"
electrons4skim.cut = "abs(eta) < 2.5"

SKIM_SEQUENCE = 'reskimNoOS'   # our default (20/10  and a 40 SF pair with mll > 40 with no OS requirement)
