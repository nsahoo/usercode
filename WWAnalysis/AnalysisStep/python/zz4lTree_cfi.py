import FWCore.ParameterSet.Config as cms

zz4lTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("skimEvent4L"),
    cut = cms.string("mass > 0"),
    addRunLumiInfo = cms.bool(True),
    variables = cms.PSet(
        hypo = cms.string("hypo()"),
        channel = cms.string("channel()"),
        mass = cms.string("mass"),
        pt   = cms.string("pt"),
        rap  = cms.string("rapidity"),
        z1mass = cms.string("z(0).mass"),
        z1pt   = cms.string("z(0).pt"),
        z2mass = cms.string("z(1).mass"),
        z2pt   = cms.string("z(1).pt"),
        elPtMin = cms.string("elePtMin"),
        muPtMin  = cms.string("muPtMin"),
        pfmet    = cms.string("pfMet.pt()"),
        rho    = cms.string("getRho()"),
        worsePairCombRelIsoBaseline = cms.string("worsePairCombRelIsoBaseline()"),
        intimeSimVertices = cms.string("intimeSimVertices()"),
        recoVertices = cms.string("getNumRecoVertices()")
    ),
    flags = cms.PSet(
        # ===========================
        muId11Hits   = cms.string("nGoodLeptons('track.numberOfValidHits > 10', '') == 4"),
        muIdGlbOr2st = cms.string("nGoodLeptons('isGlobalMuon || numberOfMatchedStations > 1', '') == 4"),
        muIdBaseline  = cms.string("nGoodLeptons('isGlobalMuon && track.numberOfValidHits > 10', '') == 4"),
        muIdPF       = cms.string("nGoodLeptons('userInt(\"pfMuId\") > 0', '') == 4"),
        # ===========================
        elIdLoose   = cms.string("nGoodLeptons('', 'test_bit(electronID(\"cicLoose\"),0)  == 1') == 4"),
        elIdMedium  = cms.string("nGoodLeptons('', 'test_bit(electronID(\"cicMedium\"),0) == 1') == 4"),
        elIdBaseline = cms.string("nGoodLeptons('', 'test_bit(electronID(\"cicTight\"),0)  == 1') == 4"),
        elIdBDT     = cms.string("nGoodLeptons('', 'abs(eta) < 1.479 && electronID(\"dbBDT\") > -0.11 || "+
                                                   "abs(eta) > 1.479 && electronID(\"dbBDT\") > -0.05') == 4"),
        # ===========================
        #  For the 4e and 4mu final states, at least three of the four combinations of opposite-sign pairs must satisfy m(ll) > 12 GeV
        threeMassCutBaseline = cms.string("hypo == 4 || hypo == 5 || nGoodPairs(\"mass > 12\", 0) >= 3"),
        # ===========================
        mc_l1 = cms.string("genl(0,0).isNonnull"),
        mc_l2 = cms.string("genl(0,1).isNonnull"),
        mc_l3 = cms.string("genl(1,0).isNonnull"),
        mc_l4 = cms.string("genl(1,1).isNonnull"),
        mc_z1 = cms.string("genz(0).isNonnull"),
        mc_z2 = cms.string("genz(0).isNonnull"),
        mc_4lany = cms.string("nMatchedLeptons == 4"),
        mc_4l    = cms.string("hasCorrectAssignment"),
    ),
)
def addLepVar(zz4l, name, expr) :
    setattr(zz4l.variables, "l1"+name, cms.string(expr % (0, 0)));
    setattr(zz4l.variables, "l2"+name, cms.string(expr % (0, 1)));
    setattr(zz4l.variables, "l3"+name, cms.string(expr % (1, 0)));
    setattr(zz4l.variables, "l4"+name, cms.string(expr % (1, 1)));
def addLepCut(zz4l, name, expr) :
    setattr(zz4l.flags, "l1"+name, cms.string(expr % (0, 0)));
    setattr(zz4l.flags, "l2"+name, cms.string(expr % (0, 1)));
    setattr(zz4l.flags, "l3"+name, cms.string(expr % (1, 0)));
    setattr(zz4l.flags, "l4"+name, cms.string(expr % (1, 1)));
addLepVar(zz4lTree, "pt", "lpt(%d,%d)")
addLepVar(zz4lTree, "eta", "leta(%d,%d)")
addLepVar(zz4lTree, "phi", "lphi(%d,%d)")
addLepVar(zz4lTree, "q", "lq(%d,%d)")
addLepVar(zz4lTree, "dz", "ldz(%d,%d)")
addLepVar(zz4lTree, "ip2d", "lip2d(%d,%d)")
addLepVar(zz4lTree, "sip3d", "lsip3d(%d,%d)")
addLepVar(zz4lTree, "pdgId", "lpdgId(%d,%d)")
addLepVar(zz4lTree, "trig",  "lfiresTrigger(%d,%d)")
addLepVar(zz4lTree, "chIso", "lisoCh(%d,%d)")
addLepVar(zz4lTree, "neuIso", "lisoNeu(%d,%d)")
addLepVar(zz4lTree, "phoIso", "lisoPho(%d,%d)")
addLepVar(zz4lTree, "trkIsoBaseline", "lisoTrkBaseline(%d,%d)")
addLepVar(zz4lTree, "ecalIsoBaseline", "lisoEcalBaseline(%d,%d)")
addLepVar(zz4lTree, "hcalIsoBaseline", "lisoHcalBaseline(%d,%d)")
addLepVar(zz4lTree, "ecalIsoBaselineRaw", "luserFloat(%d,%d, 'ecalZZ4L')")
addLepVar(zz4lTree, "hcalIsoBaselineRaw", "luserFloat(%d,%d, 'hcalZZ4L')")
addLepVar(zz4lTree, "combRelIsoBaseline", "lisoCombRelBaseline(%d,%d)")
