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
        recoVertices = cms.string("getNumRecoVertices()"),
        massErr  = cms.string("userFloat('massErr')"),
        melaLD   = cms.string("userFloat('melaSMH')/(userFloat('melaSMH')+userFloat('melaQQZZ'))"),
    ),
    flags = cms.PSet(
        # ===========================
        threeMassCut12SF = cms.string("hypo == 4 || hypo == 5 || nGoodPairs(\"mass > 12\", 0) >= 3"),
        fourMassCut4SF   = cms.string("hypo == 4 || hypo == 5 || nGoodPairs(\"mass > 4\", 0) >= 4"),
        fourMassCut4Any  = cms.string("nGoodPairs(\"mass > 4\", 0) >= 4"),
        fourMassCut4AnyS = cms.string("nGoodPairs(\"mass > 4\", 1) >= 6"),
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
addLepVar(zz4lTree, "chIso", "lisoTrkCustom(%d,%d)")
addLepVar(zz4lTree, "neuIso", "lisoNeuCustom(%d,%d)")
addLepVar(zz4lTree, "phoIso", "lisoPhoCustom(%d,%d)")
addLepVar(zz4lTree, "trkIsoBaseline", "lisoTrkBaseline(%d,%d)")
addLepVar(zz4lTree, "ecalIsoBaseline", "lisoEcalBaseline(%d,%d)")
addLepVar(zz4lTree, "hcalIsoBaseline", "lisoHcalBaseline(%d,%d)")
addLepVar(zz4lTree, "ecalIsoBaselineRaw", "luserFloat(%d,%d, 'ecalZZ4L')")
addLepVar(zz4lTree, "hcalIsoBaselineRaw", "luserFloat(%d,%d, 'hcalZZ4L')")
addLepVar(zz4lTree, "combRelIsoBaseline", "lisoCombRelBaseline(%d,%d)")
addLepVar(zz4lTree, "combRelIsoBaseline", "lisoCombRelBaseline(%d,%d)")
addLepVar(zz4lTree, "idNew", "luserInt(%d,%d,'newID')")
addLepVar(zz4lTree, "idPRL", "luserInt(%d,%d,'prlID')")
addLepVar(zz4lTree, "pfIsoChHad04",      "luserFloat(%d,%d, 'pfChHadIso04')")
addLepVar(zz4lTree, "pfIsoComb04EACorr", "luserFloat(%d,%d, 'pfCombIso04EACorr')")
addLepVar(zz4lTree, "bdtIso", "luserFloat(%d,%d, 'bdtIso')")
