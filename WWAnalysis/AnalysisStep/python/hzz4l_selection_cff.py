import FWCore.ParameterSet.Config as cms

MU_PT_MIN=5
EL_PT_MIN=7

MU_ID_PF="userInt('pfMuId')>0"
MU_ID_PRL="isGlobalMuon && track.numberOfValidHits > 10"

EL_BDT="bdtnontrig"
EL_ID_BDT=("(pt <= 10 && ( abs(eta) < 0.8 && userFloat('bdtID') > 0.47  ||"+
                        " 0.8 <= abs(eta) < 1.479 && userFloat('bdtID') > 0.004 || "+
                        " abs(eta) >= 1.479       && userFloat('bdtID') > 0.295) || "+
            "pt > 10 && ( abs(eta) < 0.8 && userFloat('bdtID') > 0.5  ||"+
                        " 0.8 <= abs(eta) < 1.479 && userFloat('bdtID') > 0.12 || "+
                        " abs(eta) >= 1.479       && userFloat('bdtID') > 0.6))")
EL_ID_NEW=(EL_ID_BDT+" && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")
EL_ID_PRL=("test_bit(electronID('cicTight'),0) == 1 && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")

boostedElectronsEAPFIso = cms.EDProducer("PatElectronEffAreaIso",
    src = cms.InputTag("boostedElectrons"),
    rho = cms.string("rhoEl"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    effectiveAreaTarget = cms.string("Data2011"),
)
boostedMuonsEAPFIso = cms.EDProducer("PatMuonEffAreaIso",
    src = cms.InputTag("boostedMuons"),
    rho = cms.string("rhoMu"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    effectiveAreaTarget = cms.string("Data2011"),
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
    )
)
boostedMuons = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('muonPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('muonPFIsoChHad04')"),
        bdtIso = cms.string("userFloat('bdtisonontrigPFNOPU')"),
        sip    = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    flags = cms.PSet(
        prlID = cms.string(MU_ID_PRL),
        newID = cms.string(MU_ID_PF),
    )
)

reboosting = cms.Sequence(
    boostedElectronsEAPFIso * boostedElectrons +
    boostedMuonsEAPFIso * boostedMuons 
)

#### CUT FLOW BUILDING BLOCKS
MU_PRESELECTION = ("pt > %f && abs(eta) < 2.4 && (isGlobalMuon || isTrackerMuon)" % MU_PT_MIN)
EL_PRESELECTION = ("pt > %f && abs(eta) < 2.5 " % EL_PT_MIN)
SINGLE_SIP_CUT       = "userFloat('sip') < 4"
SINGLE_SIP_CUT_LOOSE = "userFloat('sip') < 100"
PAIR_SIP_CUT       = "luserFloat(0,'sip') < 4   && luserFloat(1,'sip') < 4  "
PAIR_SIP_CUT_LOOSE = "luserFloat(0,'sip') < 100 && luserFloat(1,'sip') < 100"

SINGLE_ID_NEW = "userInt('newID')"
PAIR_ID_NEW   = "luserInt(0, 'newID') && luserInt(1, 'newID')"

SINGLE_PFISO_1D_LOOSE="userFloat('pfChHadRelIso04') < 0.7"
SINGLE_PFISO_1D="userFloat('pfCombRelIso04EACorr') < 0.2"
PAIR_PFISO_1D="luserFloat(0,'pfCombRelIso04EACorr') < 0.25 && luserFloat(1,'pfCombRelIso04EACorr') < 0.2"
PAIR_PFISO_2D="luserFloat(0,'pfCombRelIso04EACorr') + luserFloat(1,'pfCombRelIso04EACorr') < 0.35"

#### CUT FLOW
MUID_LOOSE = " && ".join([MU_PRESELECTION, SINGLE_SIP_CUT_LOOSE])
ELID_LOOSE = " && ".join([EL_PRESELECTION, SINGLE_SIP_CUT_LOOSE])

MUID_GOOD = " && ".join([SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])
ELID_GOOD = " && ".join([SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])

SEL_ANY_Z = "mass > 4 && daughter(0).pdgId = - daughter(1).pdgId"

SEL_ZZ4L_PRE_BEST_1 = "nGoodPairs(\"mass > 4\", 1) >= 6"
SEL_ZZ4L_PRE_BEST_2 = "lByPt(0).pt > 20 && lByPt(1).pt > 10"
SEL_ZZ4L_PRE_BEST_3 = "40 < mz(0) < 120"
SEL_ZZ4L_PRE_BEST_4 = "mass > 100"

SEL_ZZ4L_ARBITRATION_1 = "-abs(mz(0)-91.188)"
SEL_ZZ4L_ARBITRATION_2 = "daughter(1).daughter(0).pt + daughter(1).daughter(1).pt"

SEL_ZZ4L_POST_BEST_1 = "mz(0) > 50 && mz(1) > 12"
SEL_ZZ4L_POST_BEST_2 = ""
SEL_ZZ4L_POST_BEST_3 = ""
SEL_ZZ4L_POST_BEST_4 = ""
SEL_ZZ4L_POST_FILTER = True

#### CUTS RELATIVE TO CONTROL REGION ONLY
SEL_BEST_Z1 = ""

MUID_LOOSE_CR = " && ".join([MU_PRESELECTION, SINGLE_SIP_CUT])
ELID_LOOSE_CR = " && ".join([EL_PRESELECTION, SINGLE_SIP_CUT])


SEL_ANY_CR_Z = "mass > 4"
