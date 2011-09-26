import FWCore.ParameterSet.Config as cms

nverticesModule = cms.EDProducer("VertexMultiplicityCounter",
    probes = cms.InputTag("REPLACE_ME"),
    objects = cms.InputTag("goodPrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
)

# "( (abs(pdgId(0))==11) * ( 0.998 * (pt(0) < 20 && abs(eta(0)) < 1.479) + "   + 
# "                          0.948 * (pt(0) > 20 && abs(eta(0)) < 1.479) + "   + 
# "                          0.948 * (pt(0) < 20 && abs(eta(0)) > 1.479) + "   + 
# "                          0.948 * (pt(0) > 20 && abs(eta(0)) > 1.479) ) + " + 
# "  (abs(pdgId(0))==13) * ( 0.998 * (pt(1) < 20 && abs(eta(1)) < 1.479) + "   + 
# "                          0.948 * (pt(1) > 20 && abs(eta(1)) < 1.479) + "   + 
# "                          0.948 * (pt(1) < 20 && abs(eta(1)) > 1.479) + "   + 
# "                          0.948 * (pt(1) > 20 && abs(eta(1)) > 1.479) ) ) * " + 
# "( (abs(pdgId(1))==11) * ( 0.998 * (pt(0) < 20 && abs(eta(0)) < 1.479) + "   + 
# "                          0.948 * (pt(0) > 20 && abs(eta(0)) < 1.479) + "   + 
# "                          0.948 * (pt(0) < 20 && abs(eta(0)) > 1.479) + "   + 
# "                          0.948 * (pt(0) > 20 && abs(eta(0)) > 1.479) ) + " + 
# "  (abs(pdgId(1))==13) * ( 0.998 * (pt(1) < 20 && abs(eta(1)) < 1.479) + "   + 
# "                          0.948 * (pt(1) > 20 && abs(eta(1)) < 1.479) + "   + 
# "                          0.948 * (pt(1) < 20 && abs(eta(1)) > 1.479) + "   + 
# "                          0.948 * (pt(1) > 20 && abs(eta(1)) > 1.479) ) )" + 

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
        channel = cms.string("channel()"),
        mll  = cms.string("mll()"),
        ptll = cms.string("pTll()"),
        yll  = cms.string("yll()"),  #fixed! returns (p4a+p4b).Rapidity()
        pt1  = cms.string("ptMax"),
        pt2  = cms.string("ptMin"),
        peaking  = cms.string("peaking"),
        trigger  = cms.string("guillelmoTrigger('DATASET')"),
        nextra  = cms.string("nExtraLep(10)"),
        tcmet  = cms.string("tcMet"),
        tcmetphi  = cms.string("tcMetPhi"),
        ptcmet = cms.string("projTcMet"),
        pfmet  = cms.string("pfMet"),
        pfmetphi  = cms.string("pfMetPhi"),
        ppfmet = cms.string("projPfMet"),
        chmet = cms.string("chargedMetSmurf"), 
        chmetphi = cms.string("chargedMetSmurfPhi"), 
        pchmet = cms.string("projChargedMetSmurf"), 
        redmet = cms.string("-9999"), 
        predmet = cms.string("-9999"), 
        mpmet = cms.string("min(projPfMet,projChargedMetSmurf)"), ##note: min of proj and proj of min are not the same
        dphill = cms.string("dPhill()"),
        drll   = cms.string("dRll()"),
        dphilljet  = cms.string("dPhillLeadingJet(5.0)"),
        dphillmet  = cms.string("dPhillMet('PFMET')"),
        dphilmet = cms.string("dPhilMet('PFMET')"),
        dphilmet1 = cms.string("dPhilMet(0,'PFMET')"),
        dphilmet2 = cms.string("dPhilMet(1,'PFMET')"),
        mtw1 = cms.string("mTByPt(0,'PFMET')"),
        mtw2 = cms.string("mTByPt(1,'PFMET')"),
        mth  = cms.string("mTHiggs('PFMET')"),
        gammaMRStar = cms.string("gammaMRStar"),
        njet  = cms.string("nCentralJets(30,5.0)"),
        njetid  = cms.string("nCentralJets(30,5.0,1,0)"),
        nbjet = cms.string("bTaggedJetsOver(30,2.1)"),
        jetpt1 = cms.string("leadingJetPt(0,0,5.0)"),
        jetpt2 = cms.string("leadingJetPt(1,0,5.0)"),
        jeteta1 = cms.string("leadingJetEta(0,0,5.0)"),
        jeteta2 = cms.string("leadingJetEta(1,0,5.0)"),
        jetphi1 = cms.string("leadingJetPhi(0,0,5.0)"),
        jetphi2 = cms.string("leadingJetPhi(1,0,5.0)"),
        jettche1 = cms.string("leadingJetBtag(0,'trackCountingHighEffBJetTags',0,5.0)"),
        jettche2 = cms.string("leadingJetBtag(1,'trackCountingHighEffBJetTags',0,5.0)"),
        jettchp1 = cms.string("leadingJetBtag(0,'trackCountingHighPurBJetTags',0,5.0)"),
        jettchp2 = cms.string("leadingJetBtag(1,'trackCountingHighPurBJetTags',0,5.0)"),
        #jetcsv1 = cms.string("leadingJetBtag(0,'combinedSecondaryVertexBJetTags',0,5.0)"),
        #jetcsv2 = cms.string("leadingJetBtag(1,'combinedSecondaryVertexBJetTags',0,5.0)"),
        #jetcsvm1 = cms.string("leadingJetBtag(0,'combinedSecondaryVertexMVABJetTags',0,5.0)"),
        #jetcsvm2 = cms.string("leadingJetBtag(1,'combinedSecondaryVertexMVABJetTags',0,5.0)"),
        #jetjbpb1 = cms.string("leadingJetBtag(0,'jetBProbabilityBJetTags',0,5.0)"),
        #jetjbpb2 = cms.string("leadingJetBtag(1,'jetBProbabilityBJetTags',0,5.0)"),
        #jetjpb1 = cms.string("leadingJetBtag(0,'jetProbabilityBJetTags',0,5.0)"),
        #jetjpb2 = cms.string("leadingJetBtag(1,'jetProbabilityBJetTags',0,5.0)"),
        iso1 = cms.string("allIsoByPt(0)/ptByPt(0)"),
        iso2 = cms.string("allIsoByPt(1)/ptByPt(1)"),
        eta1 = cms.string("etaByPt(0)"),
        eta2 = cms.string("etaByPt(1)"),
        phi1 = cms.string("phiByPt(0)"),
        phi2 = cms.string("phiByPt(1)"),
        ch1 = cms.string("qByPt(0)"),
        ch2 = cms.string("qByPt(1)"),
        softbdisc = cms.string("highestSoftBDisc(30.0)"),
        hardbdisc = cms.string("highestHardBDisc(30.0)"),
        #softtche = cms.string("highestSoftBDisc(30.0,'trackCountingHighEffBJetTags')"),
        #hardtche = cms.string("highestHardBDisc(30.0,'trackCountingHighEffBJetTags')"),
        #softtchp = cms.string("highestSoftBDisc(30.0,'trackCountingHighPurBJetTags')"),
        #hardtchp = cms.string("highestHardBDisc(30.0,'trackCountingHighPurBJetTags')"),
        #softcsv  = cms.string("highestSoftBDisc(30.0,'combinedSecondaryVertexBJetTags')"),
        #hardcsv  = cms.string("highestHardBDisc(30.0,'combinedSecondaryVertexBJetTags')"),
        #softcsvm = cms.string("highestSoftBDisc(30.0,'combinedSecondaryVertexMVABJetTags')"),
        #hardcsvm = cms.string("highestHardBDisc(30.0,'combinedSecondaryVertexMVABJetTags')"),
        #softjbpb = cms.string("highestSoftBDisc(30.0,'jetBProbabilityBJetTags')"),
        #hardjbpb = cms.string("highestHardBDisc(30.0,'jetBProbabilityBJetTags')"),
        #softjpb  = cms.string("highestSoftBDisc(30.0,'jetProbabilityBJetTags')"),
        #hardjpb  = cms.string("highestHardBDisc(30.0,'jetProbabilityBJetTags')"),
        tightmu = cms.string("passesSmurfMuonID"),
        worstJetLepPt = cms.string("max(matchedJetPt(0, 0.5)/pt(0), matchedJetPt(1, 0.5)/pt(1))"),
        dataset = cms.string("REPLACE_ME"),
        puW   = cms.InputTag("puWeight"),
        kfW   = cms.InputTag("ptWeight"),
        baseW = cms.string("REPLACE_ME"),
        fourW = cms.string("REPLACE_ME"),
        effW = cms.string("1"),
        triggW = cms.string("1"),
        #vbf stuff:
        njetvbf = cms.string("nJetVBF(30,5.)"),
        mjj = cms.string("mjj(30,5.)"),
        detajj = cms.string("dEtajj(30,5.)"),
        #zep
    ),
    flags = cms.PSet(
        sameflav   = cms.string("hypo == 3 || hypo == 6"),
        zveto      = cms.string("abs(mll-91.1876)>15 || hypo == 4 || hypo == 5"),
        bveto      = cms.string("bTaggedJetsUnder(30,2.1) == 0 && nSoftMu(3) == 0"),
        bveto_ip   = cms.string("bTaggedJetsUnder(30,2.1) == 0"),
        bveto_mu   = cms.string("nSoftMu(3) == 0"),
        bveto_nj   = cms.string("bTaggedJetsUnder(30,2.1) == 0 && nSoftMu(3,1) == 0"),
        bveto_munj = cms.string("nSoftMu(3,1) == 0"),
        dphiveto   = cms.string("passesDPhillJet"),
        passTight1 = cms.string('passTightByPt(0)'),
        passTight2 = cms.string('passTightByPt(1)'),
        #passLoose1 = cms.string('passLooseByPt(0)'),
        #passLoose2 = cms.string('passLooseByPt(1)'),        
    ),
    addRunLumiInfo = cms.bool(True)
)

# from WWAnalysis.AnalysisStep.pileupReweighting_cfi import reWeightVector
from WWAnalysis.AnalysisStep.pileupReweighting_cfi import puS4,dataWeights
puWeight     = cms.EDProducer("CombinedWeightProducer",
    baseWeight = cms.double(1.0),
#     puWeight   = cms.vdouble(*reWeightVector[:]),
#     puLabel    = cms.InputTag("addPileupInfo"),
    s4Dist = cms.vdouble(puS4[:]),
    dataDist = cms.vdouble(dataWeights[:]),
    src        = cms.InputTag("REPLACE_ME"),
)
higgsPt = cms.EDProducer("HWWKFactorProducer",
    genParticlesTag = cms.InputTag("prunedGen"),
    inputFilename = cms.untracked.string("REPLACE_ME"),
    ProcessID = cms.untracked.int32(10010),
    Debug =cms.untracked.bool(False)
)
ptWeight     = cms.EDProducer("CombinedWeightProducer",
    baseWeight = cms.double(1.0),
    ptWeight   = cms.InputTag("higgsPt"),
    src        = cms.InputTag("REPLACE_ME"),
)
