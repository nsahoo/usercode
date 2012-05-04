import FWCore.ParameterSet.Config as cms

MU_PT_MIN=5
EL_PT_MIN=7

MU_PRESELECTION = ("pt > %f && abs(eta) < 2.4 && (isGlobalMuon || isTrackerMuon)" % MU_PT_MIN)
EL_PRESELECTION = ("pt > %f && abs(eta) < 2.5 " % EL_PT_MIN)
SIP_CUT = "userFloat('ip')/userFloat('ipErr') < 4"

MU_ID_PF="userInt('pfMuId')>0"
MU_ID_PRL="isGlobalMuon && track.numberOfValidHits > 10"

EL_BDT="bdtnontrig"
EL_ID_BDT=("(pt <= 10 && ( abs(eta) < 0.8 && userFloat('{0}') > 0.47  ||"+
                        " 0.8 <= abs(eta) < 1.479 && userFloat('{0}') > 0.004 || "+
                        " abs(eta) >= 1.479       && userFloat('{0}') > 0.295) || "+
            "pt > 10 && ( abs(eta) < 0.8 && userFloat('{0}') > 0.5  ||"+
                        " 0.8 <= abs(eta) < 1.479 && userFloat('{0}') > 0.12 || "+
                        " abs(eta) >= 1.479       && userFloat('{0}') > 0.6))".format(EL_BDT))
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
        pfChHadIso04 = cms.string("userFloat('electronPFIsoChHad04')"),
        bdtID  = cms.string("userFloat('%s')" % EL_BDT),
        bdtIso = cms.string("userFloat('bdtisonontrig')"),
        sip   = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    cuts = cms.PSet(
        bdtID = cms.string(EL_ID_BDT), 
        newID = cms.string(EL_ID_BDT), 
        prlID = cms.string(EL_ID_PRL),
    )
)
boostedMuons = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsEAPFIso"),
    variables = cms.PSet(
        pfChHadIso04 = cms.string("userFloat('muonPFIsoChHad04')"),
        bdtIso = cms.string("userFloat('bdtisonontrigPFNOPU')"),
        sip    = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    cuts = cms.PSet(
        prlID = cms.string(MU_ID_PRL),
        newID = cms.string(MU_ID_PF),
    )
)

reboosting = cms.Sequence(
    boostedElectronsEAPFIso * boostedElectrons +
    boostedMuonsEAPFIso * boostedMuons 
)

