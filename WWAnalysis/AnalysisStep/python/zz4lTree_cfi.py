import FWCore.ParameterSet.Config as cms

zz4lTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("skimEvent4L"),
    cut = cms.string("mass > 0"),
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
        pfmet  = cms.string("pfMet.pt()"),
    ),
    flags = cms.PSet(
        # ===========================
        muId11Hits   = cms.string("nGoodLeptons('track.numberOfValidHits > 10', '') == 4"),
        muIdGlbOr2st = cms.string("nGoodLeptons('isGlobalMuon || numberOfMatchedStations > 1', '') == 4"),
        muIdNominal  = cms.string("nGoodLeptons('isGlobalMuon && track.numberOfValidHits > 10', '') == 4"),
        # ===========================
        elIdLoose   = cms.string("nGoodLeptons('', 'test_bit(electronID(\"cicLoose\"),0)  == 1') == 4"),
        elIdMedium  = cms.string("nGoodLeptons('', 'test_bit(electronID(\"cicMedium\"),0) == 1') == 4"),
        elIdNominal = cms.string("nGoodLeptons('', 'test_bit(electronID(\"cicTight\"),0)  == 1') == 4"),
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
addLepVar(zz4lTree, "ip3d", "lip3d(%d,%d)")
addLepVar(zz4lTree, "sip2d", "lsip2d(%d,%d)")
addLepVar(zz4lTree, "sip3d", "lsip3d(%d,%d)")
addLepVar(zz4lTree, "isoCh",   "userFloat('iso4lCharged_%d_%d')")
addLepVar(zz4lTree, "isoChNC", "userFloat('iso4lChargedNC_%d_%d')")
addLepVar(zz4lTree, "isoNeu",  "userFloat('iso4lNeutral_%d_%d')")
addLepVar(zz4lTree, "isoPho",  "userFloat('iso4lPhotons_%d_%d')")
addLepVar(zz4lTree, "isoPhoC", "userFloat('iso4lPhotonsC_%d_%d')")
