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

## ==== Alternative set of cuts optimized for MVA iso (and electron ID) from Si and company
MU_MVA_ISO=("(isGlobalMuon && isTrackerMuon && (" +
                 "pt <  10 && abs(eta) <  1.5 && userFloat('bdtIso') > -0.593 ||"+
                 "pt >= 10 && abs(eta) <  1.5 && userFloat('bdtIso') >  0.337 ||"+
                 "pt <  10 && abs(eta) >= 1.5 && userFloat('bdtIso') > -0.767 ||"+
                 "pt >= 10 && abs(eta) >= 1.5 && userFloat('bdtIso') >  0.410) ||"+
              " !isGlobalMuon && isTrackerMuon && userFloat('bdtIso') > -0.989 ||"+
              "isGlobalMuon &&  !isTrackerMuon && userFloat('bdtIso') > -0.995)")


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

EL_MVA_ID_TIGHT=("(pt <  10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.093 ||"+
                              " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') >  0.451 || "+
                              "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.595 ) || "+
                  "pt >= 10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.881 ||"+
                              " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') >  0.731 || "+
                              "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.819 ))")
EL_MVA_ISO_TIGHT=("(pt < 10 && (       abs(eta) <  0.8    && userFloat('bdtIso') >  0.553 ||"+
                               "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') > -0.237 || "+
                               "       abs(eta) >= 1.479  && userFloat('bdtIso') > -0.573) || "+
                   "pt >= 10 && (      abs(eta) <  0.8    && userFloat('bdtIso') >  0.521||"+
                               "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') >  0.531|| "+
                               "       abs(eta) >= 1.479  && userFloat('bdtIso') >  0.493))")



## Re-run Muon BDT Iso (tag has changed w.r.t. skim production)
boostedMuonsUpdatedBDTIso = cms.EDProducer("PatMuonBoosterBDTIso", 
      useExistingIsoValues = cms.bool(True),
      src = cms.InputTag("boostedMuons"),
      rho = cms.string("rhoMu"),
      outputName = cms.string("bdtisoNew"),
      effectiveAreaTarget = cms.string("Data2011"),
)

## Re-run POG PF isolation for neutrals with threshold at pt > 0.5 GeV
updatedMuonPFIsoNHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuonsUpdatedBDTIso"),
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 130 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
updatedMuonPFIsoPhoton04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuonsUpdatedBDTIso"),
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 22 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
boostedMuonsUpdatedPFIso = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsUpdatedBDTIso"),
    valueMaps = cms.PSet(
        muonPFIsoPhoton04pt05_step1 = cms.InputTag("updatedMuonPFIsoPhoton04"),
        muonPFIsoNHad04pt05_step1   = cms.InputTag("updatedMuonPFIsoNHad04"),
    ),
)

fixupMuonSequence = cms.Sequence(
    boostedMuonsUpdatedBDTIso *
    updatedMuonPFIsoNHad04 * updatedMuonPFIsoPhoton04 * boostedMuonsUpdatedPFIso
)

#===========================================================================

## Recompute iso without vetos that are not in POG code
updatedElectronPFIsoChHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedElectrons"),
    endcapDefinition = cms.string("abs(superCluster.eta) > 1.479"), # for electrons we use supercluster eta, not the default momentum eta
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge != 0 && abs(pdgId) == 211"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.015), 
    directional = cms.bool(False),
)
updatedElectronPFIsoNHad04 = updatedElectronPFIsoChHad04.clone(
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 130"), # neutral hadrons
    vetoConeEndcaps = cms.double(0.0), # No veto for NH, EGamma has abandoned it
)
updatedElectronPFIsoPhoton04 = updatedElectronPFIsoChHad04.clone(
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 22"), # neutral hadrons
    vetoConeEndcaps = cms.double(0.08), # Preserve veto for Photons
)
updatedElectronPFIsoChHad01 = updatedElectronPFIsoChHad04.clone(deltaR = 0.1)
updatedElectronPFIsoChHad02 = updatedElectronPFIsoChHad04.clone(deltaR = 0.2)
updatedElectronPFIsoChHad03 = updatedElectronPFIsoChHad04.clone(deltaR = 0.3)
updatedElectronPFIsoChHad05 = updatedElectronPFIsoChHad04.clone(deltaR = 0.5)
updatedElectronPFIsoNHad01 = updatedElectronPFIsoNHad04.clone(deltaR = 0.1)
updatedElectronPFIsoNHad02 = updatedElectronPFIsoNHad04.clone(deltaR = 0.2)
updatedElectronPFIsoNHad03 = updatedElectronPFIsoNHad04.clone(deltaR = 0.3)
updatedElectronPFIsoNHad05 = updatedElectronPFIsoNHad04.clone(deltaR = 0.5)
updatedElectronPFIsoPhoton01 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.1)
updatedElectronPFIsoPhoton02 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.2)
updatedElectronPFIsoPhoton03 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.3)
updatedElectronPFIsoPhoton05 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.5)

boostedElectronsUpdatedPFIso = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectrons"),
    valueMaps = cms.PSet(
        electronPFIsoChHad01_v2 = cms.InputTag("updatedElectronPFIsoChHad01"),
        electronPFIsoChHad02_v2 = cms.InputTag("updatedElectronPFIsoChHad02"),
        electronPFIsoChHad03_v2 = cms.InputTag("updatedElectronPFIsoChHad03"),
        electronPFIsoChHad04_v2 = cms.InputTag("updatedElectronPFIsoChHad04"),
        electronPFIsoChHad05_v2 = cms.InputTag("updatedElectronPFIsoChHad05"),
        electronPFIsoNHad01_v2 = cms.InputTag("updatedElectronPFIsoNHad01"),
        electronPFIsoNHad02_v2 = cms.InputTag("updatedElectronPFIsoNHad02"),
        electronPFIsoNHad03_v2 = cms.InputTag("updatedElectronPFIsoNHad03"),
        electronPFIsoNHad04_v2 = cms.InputTag("updatedElectronPFIsoNHad04"),
        electronPFIsoNHad05_v2 = cms.InputTag("updatedElectronPFIsoNHad05"),
        electronPFIsoPhoton01_v2 = cms.InputTag("updatedElectronPFIsoPhoton01"),
        electronPFIsoPhoton02_v2 = cms.InputTag("updatedElectronPFIsoPhoton02"),
        electronPFIsoPhoton03_v2 = cms.InputTag("updatedElectronPFIsoPhoton03"),
        electronPFIsoPhoton04_v2 = cms.InputTag("updatedElectronPFIsoPhoton04"),
        electronPFIsoPhoton05_v2 = cms.InputTag("updatedElectronPFIsoPhoton05"),
    )
)

## Fix BDT ID bug (userInt vs userFloat)
boostedElectronsFixBDTID = cms.EDProducer("PatElectronBoosterBDTID", 
    src = cms.InputTag("boostedElectronsUpdatedPFIso"),
    postfix = cms.string("_Fix"),
)
boostedElectronsFixBDTIso = cms.EDProducer("PatElectronBoosterBDTIso", 
    src = cms.InputTag("boostedElectronsFixBDTID"),
    effectiveAreaTarget = cms.string("Data2011"),
    chargedOption = cms.string("_v2"), # postfix to the userFloat values
    neutralsOption = cms.string("_v2"), # postfix to the userFloat value
    outputName = cms.string("bdtisoFix"),
)

fixupElectronSequence = cms.Sequence(
    ( updatedElectronPFIsoPhoton01 + updatedElectronPFIsoPhoton02 + updatedElectronPFIsoPhoton03 + updatedElectronPFIsoPhoton04 + updatedElectronPFIsoPhoton05 +
      updatedElectronPFIsoNHad01 + updatedElectronPFIsoNHad02 + updatedElectronPFIsoNHad03 + updatedElectronPFIsoNHad04 + updatedElectronPFIsoNHad05 +
      updatedElectronPFIsoChHad01 + updatedElectronPFIsoChHad02 + updatedElectronPFIsoChHad03 + updatedElectronPFIsoChHad04 + updatedElectronPFIsoChHad05) *
    boostedElectronsUpdatedPFIso * 
    boostedElectronsFixBDTID * 
    boostedElectronsFixBDTIso  
)
 
##================================================================================================0

## Compute EA-corrected isolations
boostedElectronsEAPFIso = cms.EDProducer("PatElectronEffAreaIso",
    src = cms.InputTag("boostedElectronsFixBDTIso"),
    rho = cms.string("rhoEl"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    chargedOption = cms.string("_v2"), # postfix to the userFloat values
    neutralsOption = cms.string("_v2"), # postfix to the userFloat value
    effectiveAreaTarget = cms.string("Data2011"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
)
boostedMuonsEAPFIso = cms.EDProducer("PatMuonEffAreaIso",
    src = cms.InputTag("boostedMuonsUpdatedPFIso"),
    rho = cms.string("rhoMu"),
    deltaR = cms.string("04"),
    neutralsOption = cms.string("pt05_step1"), # postfix to the userFloat value
    label = cms.string("pfCombIso04EACorr"),
    effectiveAreaTarget = cms.string("Data2011"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
)

##================================================================================================0

## Embed final values and discriminators. Beyond here, fixups should be invisible
boostedElectrons = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectronsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('electronPFIsoChHad04_v2')/pt"),
        pfChHadIso04 = cms.string("userFloat('electronPFIsoChHad04_v2')"),
        pfNHadIso04 = cms.string("userFloat('electronPFIsoNHad04_v2')"),
        pfPhotonIso04 = cms.string("userFloat('electronPFIsoPhoton04_v2')"),
        bdtID     = cms.string("userFloat('%s_Fix')" % EL_BDT),
        bdtIDTrig = cms.string("userFloat('bdttrig_Fix')"),
        #bdtIso = cms.string("userFloat('bdtisonontrig')"),
        bdtIso = cms.string("userFloat('bdtisoFix')"),
        sip   = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    flags = cms.PSet(
        bdtID = cms.string(EL_ID_BDT), 
        newID = cms.string(EL_ID_NEW), 
        prlID = cms.string(EL_ID_PRL),
        mvaID  = cms.string(EL_MVA_ID), 
        mvaIso = cms.string(EL_MVA_ISO), 
        mvaIDTight  = cms.string(EL_MVA_ID_TIGHT), 
        mvaIsoTight = cms.string(EL_MVA_ISO_TIGHT), 
    )
)
boostedMuons = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('muonPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('muonPFIsoChHad04')"),
        pfNHadIso04 = cms.string("userFloat('muonPFIsoNHad04pt05_step1')"),
        pfPhotonIso04 = cms.string("userFloat('muonPFIsoPhoton04pt05_step1')"),
        bdtIso = cms.string("userFloat('bdtisoNew')"),
        #bdtIso = cms.string("userFloat('bdtisonontrigPFNOPU')"),
        sip    = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    flags = cms.PSet(
        prlID = cms.string(MU_ID_PRL),
        newID = cms.string(MU_ID_PF),
        mvaIso = cms.string(MU_MVA_ISO), 
    )
)

reboosting = cms.Sequence(
    fixupElectronSequence     *
    boostedElectronsEAPFIso   *  boostedElectrons +
    fixupMuonSequence   *
    boostedMuonsEAPFIso * boostedMuons 
)

#### CUT FLOW BUILDING BLOCKS

SKIM_SEQUENCE = 'reskimNoOS'   # their default (20/10  and a 40 SF pair with mll > 40 but no OS request)
TRIGGER_FILTER = None 

EL_PT_MIN=7;  EL_PT_MIN_LOOSE=5
MU_PT_MIN=5;  MU_PT_MIN_LOOSE=3

MU_PRESELECTION = ("abs(eta) < 2.4 && (isGlobalMuon || isTrackerMuon)") ## the PT cut
EL_PRESELECTION = ("abs(eta) < 2.5")                                    ## is below
SINGLE_SIP_CUT       = "userFloat('sip') < 4"
SINGLE_SIP_CUT_LOOSE = "userFloat('sip') < 100"
PAIR_SIP_CUT       = "luserFloat(0,'sip') < 4   && luserFloat(1,'sip') < 4  "
PAIR_SIP_CUT_LOOSE = "luserFloat(0,'sip') < 100 && luserFloat(1,'sip') < 100"

SINGLE_ID_NEW = "userInt('newID')"
SINGLE_ID_MVA = "userInt('mvaID')" # to go jointly with MVA iso 
SINGLE_ID_MVA_TIGHT = "userInt('mvaIDTight')" # to go jointly with MVA iso 
PAIR_ID_NEW   = "luserInt(0, 'newID') && luserInt(1, 'newID')"

SINGLE_PFISO_1D_LOOSE="userFloat('pfChHadRelIso04') < 0.7"
SINGLE_PFISO_1D="userFloat('pfCombRelIso04EACorr') < 0.25"
SINGLE_MVA_ISO="userInt('mvaIso')"
SINGLE_MVA_ISO_TIGHT="userInt('mvaIsoTight')"
PAIR_PFISO_1D="luserFloat(0,'pfCombRelIso04EACorr') < 0.25 && luserFloat(1,'pfCombRelIso04EACorr') < 0.25"
PAIR_PFISO_2D="luserFloat(0,'pfCombRelIso04EACorr') + luserFloat(1,'pfCombRelIso04EACorr') < 0.35"
PAIR_DETISO_2D="combinedPairRelativeIso() < 0.35"

#### CUT FLOW
MUID_LOOSE_NO_PT_CUT = " && ".join([MU_PRESELECTION, SINGLE_SIP_CUT_LOOSE])
ELID_LOOSE_NO_PT_CUT = " && ".join([EL_PRESELECTION, SINGLE_SIP_CUT_LOOSE])
MUID_LOOSE = "pt > %f && %s" % (MU_PT_MIN_LOOSE, MUID_LOOSE_NO_PT_CUT)
ELID_LOOSE = "pt > %f && %s" % (EL_PT_MIN_LOOSE, ELID_LOOSE_NO_PT_CUT)

# for later
#MUID_MID = "(" + (" || ".join([SINGLE_PFISO_1D_LOOSE])) + ")"
#ELID_MID = "(" + (" || ".join([SINGLE_PFISO_1D_LOOSE + " && "+ SINGLE_ID_MVA, SINGLE_MVA_ISO])) + ")"

MUID_GOOD = " && ".join(["pt > %f" % MU_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])
ELID_GOOD = " && ".join(["pt > %f" % EL_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])


SEL_ANY_Z = "daughter(0).pdgId = - daughter(1).pdgId"

ARBITRATE_EARLY = True # True = PRL-logic; False = keep all candidates until the end
                       # Baseline logic is True; Current synch twiki is False

FOUR_LEPTON_FILTER_PRE_Z  = False # plug a 4-lepton count filter before making Z's
FOUR_LEPTON_FILTER_POST_Z = True # plug a 4-lepton count filter after  making Z's

SEL_BEST_Z1 = "40 < mass < 120" # Not used if ARBITRATE_EARLY = False

SEL_ZZ4L_STEP_1 = "4 < mz(1) < 120"
SEL_ZZ4L_STEP_2 = "lByPt(0).pt > 20 && lByPt(1).pt > 10"
SEL_ZZ4L_STEP_3 = "nGoodPairs(\"mass > 4\", 1) >= 6"
SEL_ZZ4L_STEP_4 = "mass > 70"
SEL_ZZ4L_STEP_5 = "mass > 100"

SEL_ZZ4L_ARBITRATION_1 = "-abs(mz(0)-91.188)"
SEL_ZZ4L_ARBITRATION_2 = "daughter(1).daughter(0).pt + daughter(1).daughter(1).pt"

#### CUTS RELATIVE TO CONTROL REGION ONLY
MUID_LOOSE_CR = " && ".join(["pt > %f" % MU_PT_MIN, MU_PRESELECTION, SINGLE_SIP_CUT])
ELID_LOOSE_CR = " && ".join(["pt > %f" % EL_PT_MIN, EL_PRESELECTION, SINGLE_SIP_CUT])

SEL_ANY_CR_Z = "mass > 4"
