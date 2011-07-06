import FWCore.ParameterSet.Config as cms

nverticesModule = cms.EDProducer("VertexMultiplicityCounter",
    probes = cms.InputTag("REPLACE_ME"),
    objects = cms.InputTag("goodPrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
)

step3Tree = cms.EDFilter("ProbeTreeProducer",
    cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && "+
                     "leptEtaCut(2.4,2.5) && ptMax > 20 && ptMin > 10"
#                      " && triggerMatchingCut('DATASET')"
#                      "nExtraLep(10) == 0 "
#                     +" && passesIP"
#                    +(" && triggerMatchingCut('DATASET')")
    ),
    variables = cms.PSet(
        hypo = cms.string("hypo()"),
        mll  = cms.string("mll()"),
        ptll = cms.string("pTll()"),
        yll  = cms.string("yll()"),  #fixed! returns (p4a+p4b).Rapidity()
        pt1  = cms.string("ptMax"),
        pt2  = cms.string("ptMin"),
        met  = cms.string("pfMet"),
        peaking  = cms.string("peaking"),
        trigmatch = cms.string("triggerMatchingCut('DATASET')"),
        trigguil  = cms.string("guillelmoTrigger('DATASET')"),
        trigbits  = cms.string("triggerBitsCut('DATASET')"),
        nextra  = cms.string("nExtraLep(10)"),
        pmet = cms.string("projPfMet"),
        pchmet = cms.string("projChargedMetSmurf"), 
        pmmet = cms.string("min(projPfMet,projChargedMetSmurf)"), ##note: min of proj and proj of min are not the same
        dphill = cms.string("dPhill()"),
        drll   = cms.string("dRll()"),
        dphilljet  = cms.string("dPhillLeadingJet(5.0)"),
        dphillcjet = cms.string("dPhillLeadingJet(2.5)"),
        dphillmet  = cms.string("dPhillMet('PFMET')"),
        dphilmet = cms.string("dPhilMet('PFMET')"),
        mtw1 = cms.string("mTByPt(0,'PFMET')"),
        mtw2 = cms.string("mTByPt(1,'PFMET')"),
        mth  = cms.string("mTHiggs('PFMET')"),
        gammaMRStar = cms.string("gammaMRStar"),
        njet  = cms.string("nCentralJets(30,5.0)"),
        ncjet = cms.string("nCentralJets(30,2.5)"),
        nbjet = cms.string("bTaggedJetsOver(30,2.1)"),
        leadjetpt = cms.string("leadingJetPt(0,5.0)"),
        leadjeteta = cms.string("leadingJetEta(0,5.0)"),
        iso1 = cms.string("allIsoByPt(0)/ptByPt(0)"),
        iso2 = cms.string("allIsoByPt(1)/ptByPt(1)"),
        eta1 = cms.string("etaByPt(0)"),
        eta2 = cms.string("etaByPt(1)"),
        softbdisc = cms.string("highestSoftBDisc(30.0)"),
        hardbdisc = cms.string("highestHardBDisc(30.0)"),
        tightmu = cms.string("passesSmurfMuonID"),
        worstJetLepPt = cms.string("max(matchedJetPt(0, 0.5)/pt(0), matchedJetPt(1, 0.5)/pt(1))"),
        dataset = cms.string("REPLACE_ME"),
        puWeight   = cms.InputTag("puWeight"),
        ptWeight   = cms.InputTag("ptWeight"),
        lumiWeight = cms.string("REPLACE_ME"),
        fourWeight = cms.string("REPLACE_ME"),
    ),
    flags = cms.PSet(
        sameflav = cms.string("hypo == 3 || hypo == 6"),
        zveto    = cms.string("abs(mll-91.1876)>15 || hypo == 4 || hypo == 5"),
        bveto    = cms.string("bTaggedJetsUnder(30,2.1) == 0 && nSoftMu(3) == 0"),
        bveto_ip = cms.string("bTaggedJetsUnder(30,2.1) == 0"),
        bveto_mu = cms.string("nSoftMu(3) == 0"),
        dphiveto = cms.string("passesDPhillJet"),
    ),
    addRunLumiInfo = cms.bool(True)
)

from WWAnalysis.AnalysisStep.pileupReweighting_cfi import reWeightVector
puWeight     = cms.EDProducer("CombinedWeightProducer",
    baseWeight = cms.double(1.0),
    puWeight   = cms.vdouble(*reWeightVector[:]),
    puLabel    = cms.InputTag("addPileupInfo"),
    src        = cms.InputTag("REPLACE_ME"),
)
higgsPt = cms.EDProducer("HWWKFactorProducer",
    genParticlesTag = cms.InputTag("onlyHiggsGen"),
    inputFilename = cms.untracked.string("REPLACE_ME"),
    ProcessID = cms.untracked.int32(10010),
    Debug =cms.untracked.bool(False)
)
ptWeight     = cms.EDProducer("CombinedWeightProducer",
    baseWeight = cms.double(1.0),
    ptWeight   = cms.InputTag("higgsPt"),
    src        = cms.InputTag("REPLACE_ME"),
)
