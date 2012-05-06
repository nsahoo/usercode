import FWCore.ParameterSet.Config as cms

MU_ID_PF="userInt('pfMuId')>0"
MU_ID_PRL="isGlobalMuon && track.numberOfValidHits > 10"

EL_BDT="bdtnontrig"
EL_ID_BDT=("(pt <= 10 && (       abs(eta) <  0.8    && userFloat('bdtID') > 0.47  ||"+
                        " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') > 0.004 || "+
                        "        abs(eta) >= 1.479  && userFloat('bdtID') > 0.295) || "+
            "pt > 10 && (        abs(eta) <  0.8    && userFloat('bdtID') > 0.5  ||"+
                        " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') > 0.12 || "+
                        "        abs(eta) >=  1.479 && userFloat('bdtID') > 0.6))")
EL_ID_NEW=(EL_ID_BDT+" && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")
EL_ID_PRL=("test_bit(electronID('cicTight'),0) == 1 && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")

## ==== Alternative set of cuts optimized for MVA iso
MU_MVA_ISO=("(isGlobalMuon && isTrackerMuon && (" +
                 "pt <  10 && abs(eta) <  1.5 && userFloat('bdtIso') > -0.593 ||"+
                 "pt >= 10 && abs(eta) <  1.5 && userFloat('bdtIso') >  0.337 ||"+
                 "pt <  10 && abs(eta) >= 1.5 && userFloat('bdtIso') > -0.767 ||"+
                 "pt >= 10 && abs(eta) >= 1.5 && userFloat('bdtIso') >  0.410) ||"+
             " isGlobalMuon && !isTrackerMuon && userFloat('bdtIso') > -0.989 ||"+
             "!isGlobalMuon &&  isTrackerMuon && userFloat('bdtIso') > -0.995)")
EL_MVA_ID=("(pt <  10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.369 ||"+
                        " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') > -0.025 || "+
                        "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.531 ) || "+
            "pt >= 10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.735 ||"+
                        " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') >  0.467 || "+
                        "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.795 ))")
EL_MVA_ISO=("(pt < 10 && (       abs(eta) <  0.8    && userFloat('bdtIso') >  0.385 ||"+
                         "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') > -0.083 || "+
                         "       abs(eta) >= 1.479  && userFloat('bdtIso') > -0.573) || "+
             "pt >= 10 && (      abs(eta) <  0.8    && userFloat('bdtIso') >  0.413||"+
                         "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') >  0.271|| "+
                         "       abs(eta) >= 1.479  && userFloat('bdtIso') >  0.135))")


boostedMuonsUpdatedBDTIso = cms.EDProducer("PatMuonBoosterBDTIso", 
      useExistingIsoValues = cms.bool(True),
      src = cms.InputTag("boostedMuons"),
      rho = cms.string("rhoMu"),
      outputName = cms.string("bdtisoNew"),
      effectiveAreaTarget = cms.string("Data2011"),
)

boostedElectronsEAPFIso = cms.EDProducer("PatElectronEffAreaIso",
    src = cms.InputTag("boostedElectrons"),
    rho = cms.string("rhoEl"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    effectiveAreaTarget = cms.string("Data2011"),
    #separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    #truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
    separatePhotonAndHadronEAs = cms.bool(True), # use separte EA for photons and hadrons
    truncateAtZero = cms.string("both"), # (yes|no) for total EA, (both|sum|no) for separate EA
)
boostedMuonsEAPFIso = cms.EDProducer("PatMuonEffAreaIso",
    src = cms.InputTag("boostedMuonsUpdatedBDTIso"),
    rho = cms.string("rhoMu"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    effectiveAreaTarget = cms.string("Data2011"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
)
boostedElectrons = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectronsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('electronPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('electronPFIsoChHad04')"),
        bdtID  = cms.string("userFloat('%s')" % EL_BDT),
        bdtIso = cms.string("userFloat('bdtisonontrig')"),
        sip   = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    flags = cms.PSet(
        bdtID = cms.string(EL_ID_BDT), 
        newID = cms.string(EL_ID_NEW), 
        prlID = cms.string(EL_ID_PRL),
        mvaID  = cms.string(EL_MVA_ID), 
        mvaIso = cms.string(EL_MVA_ISO), 
    )
)
boostedMuons = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('muonPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('muonPFIsoChHad04')"),
        bdtIso = cms.string("userFloat('bdtisoNew')"),
        sip    = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    flags = cms.PSet(
        prlID = cms.string(MU_ID_PRL),
        newID = cms.string(MU_ID_PF),
       #mvaID  = cms.string(MU_MVA_ID), # not for muons, at least not yet
        mvaIso = cms.string(MU_MVA_ISO), 
    )
)

reboosting = cms.Sequence(
    boostedElectronsEAPFIso   * boostedElectrons +
    boostedMuonsUpdatedBDTIso * boostedMuonsEAPFIso * boostedMuons 
)

#### CUT FLOW BUILDING BLOCKS
MU_PT_MIN=5
EL_PT_MIN=7

MU_PRESELECTION = ("abs(eta) < 2.4 && (isGlobalMuon || isTrackerMuon)") ## the PT cut
EL_PRESELECTION = ("abs(eta) < 2.5")                                    ## is below
SINGLE_SIP_CUT       = "userFloat('sip') < 4"
SINGLE_SIP_CUT_LOOSE = "userFloat('sip') < 100"
PAIR_SIP_CUT       = "luserFloat(0,'sip') < 4   && luserFloat(1,'sip') < 4  "
PAIR_SIP_CUT_LOOSE = "luserFloat(0,'sip') < 100 && luserFloat(1,'sip') < 100"

SINGLE_ID_NEW = "userInt('newID')"
SINGLE_ID_MVA = "userInt('mvaID')" # to go jointly with MVA iso 
PAIR_ID_NEW   = "luserInt(0, 'newID') && luserInt(1, 'newID')"

SINGLE_PFISO_1D_LOOSE="userFloat('pfChHadRelIso04') < 0.7"
SINGLE_PFISO_1D="userFloat('pfCombRelIso04EACorr') < 0.25"
SINGLE_MVA_ISO="userInt('mvaIso')"
PAIR_PFISO_1D="luserFloat(0,'pfCombRelIso04EACorr') < 0.25 && luserFloat(1,'pfCombRelIso04EACorr') < 0.25"
PAIR_PFISO_2D="luserFloat(0,'pfCombRelIso04EACorr') + luserFloat(1,'pfCombRelIso04EACorr') < 0.35"

#### CUT FLOW
MUID_LOOSE_NO_PT_CUT = " && ".join([MU_PRESELECTION, SINGLE_SIP_CUT_LOOSE])
ELID_LOOSE_NO_PT_CUT = " && ".join([EL_PRESELECTION, SINGLE_SIP_CUT_LOOSE])
MUID_LOOSE = "pt > %f && %s" % (MU_PT_MIN, MUID_LOOSE_NO_PT_CUT)
ELID_LOOSE = "pt > %f && %s" % (EL_PT_MIN, ELID_LOOSE_NO_PT_CUT)

MUID_GOOD = " && ".join([SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])
ELID_GOOD = " && ".join([SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])
#MUID_GOOD = " && ".join([SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_MVA_ISO]) ## MVA WP
#ELID_GOOD = " && ".join([SINGLE_SIP_CUT, SINGLE_ID_MVA, SINGLE_MVA_ISO]) ## MVA WP

SEL_ANY_Z = "daughter(0).pdgId = - daughter(1).pdgId"

ARBITRATE_EARLY = False # True = PRL-logic; False = keep all candidates until the end

SEL_BEST_Z1 = "" # Not used if ARBITRATE_EARLY = False

SEL_ZZ4L_STEP_1 = "lByPt(0).pt > 20 && lByPt(1).pt > 10"
SEL_ZZ4L_STEP_2 = "40 < mz(0) < 120"
SEL_ZZ4L_STEP_3 = "4 < mz(1) < 120"
SEL_ZZ4L_STEP_4 = "nGoodPairs(\"mass > 4\", 1) >= 6"
SEL_ZZ4L_STEP_5 = "mass > 100"

SEL_ZZ4L_ARBITRATION_1 = "-abs(mz(0)-91.188)"
SEL_ZZ4L_ARBITRATION_2 = "daughter(1).daughter(0).pt + daughter(1).daughter(1).pt"

#### CUTS RELATIVE TO CONTROL REGION ONLY
SEL_BEST_Z1 = ""

MUID_LOOSE_CR = " && ".join([MU_PRESELECTION, SINGLE_SIP_CUT])
ELID_LOOSE_CR = " && ".join([EL_PRESELECTION, SINGLE_SIP_CUT])


SEL_ANY_CR_Z = "mass > 4"
