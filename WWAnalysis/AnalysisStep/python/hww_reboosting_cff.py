import FWCore.ParameterSet.Config as cms

boostedMuonsUpdatedBDTIso = cms.EDProducer("PatMuonBoosterBDTIso", 
      useExistingIsoValues = cms.bool(True),
      src = cms.InputTag("boostedMuons"),
      rho = cms.string("rhoMu"),
      outputName = cms.string("bdtisoNew"),
      effectiveAreaTarget = cms.string("Data2011"),
)
updatedMuonPFIsoNHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuonsUpdatedBDTIso"),
#     pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfLabel     = cms.InputTag("reducedPFCands"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 130 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
updatedMuonPFIsoPhoton04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuonsUpdatedBDTIso"),
#     pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfLabel     = cms.InputTag("reducedPFCands"), 
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

# pfLabel: pfIso done on pfNoPU
updatedElectronPFIsoChHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedElectrons"),
    pfLabel     = cms.InputTag("reducedPFCandsPfNoPU"), 
#    pfLabel     = cms.InputTag("reducedPFCands"), 
    pfSelection = cms.string("charge != 0 && abs(pdgId) == 211"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.015), 
    directional = cms.bool(False),
)
updatedElectronPFIsoChHad01 = updatedElectronPFIsoChHad04.clone(deltaR = 0.1)
updatedElectronPFIsoChHad02 = updatedElectronPFIsoChHad04.clone(deltaR = 0.2)
updatedElectronPFIsoChHad03 = updatedElectronPFIsoChHad04.clone(deltaR = 0.3)
updatedElectronPFIsoChHad05 = updatedElectronPFIsoChHad04.clone(deltaR = 0.5)
boostedElectronsUpdatedPFIso = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectrons"),
    valueMaps = cms.PSet( 
        electronPFIsoChHad01_v2 = cms.InputTag("updatedElectronPFIsoChHad01"),
        electronPFIsoChHad02_v2 = cms.InputTag("updatedElectronPFIsoChHad02"),
        electronPFIsoChHad03_v2 = cms.InputTag("updatedElectronPFIsoChHad03"),
        electronPFIsoChHad04_v2 = cms.InputTag("updatedElectronPFIsoChHad04"),
        electronPFIsoChHad05_v2 = cms.InputTag("updatedElectronPFIsoChHad05"),
        #rhoElActiveArea = cms.InputTag("updatedElectronRhoAA"),
    )
)

boostedElectronsEAPFIso = cms.EDProducer("PatElectronEffAreaIso",
    src = cms.InputTag("boostedElectronsUpdatedPFIso"),
    rho = cms.string("rhoElFullEta"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    chargedOption = cms.string("_v2"), # postfix to the userFloat values
    effectiveAreaTarget = cms.string("Data2012"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
    #separatePhotonAndHadronEAs = cms.bool(True), # use separte EA for photons and hadrons
    #truncateAtZero = cms.string("both"), # (yes|no) for total EA, (both|sum|no) for separate EA
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
boostedElectronsFixBDTID = cms.EDProducer("PatElectronBoosterBDTID", 
    src = cms.InputTag("boostedElectronsEAPFIso"),
    postfix = cms.string("_Fix"),
)
boostedElectrons = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectronsFixBDTID"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('electronPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('electronPFIsoChHad04_v2')"),
        pfChHadIso04_old = cms.string("userFloat('electronPFIsoChHad04')"),
#         bdtID  = cms.string("userFloat('%s_Fix')" % EL_BDT),
        bdtIso = cms.string("userFloat('bdtisonontrig')"),
        sip   = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
#     flags = cms.PSet(
#         bdtID = cms.string(EL_ID_BDT), 
#         newID = cms.string(EL_ID_NEW), 
#         prlID = cms.string(EL_ID_PRL),
#         mvaID  = cms.string(EL_MVA_ID), 
#         mvaIso = cms.string(EL_MVA_ISO), 
#         mvaIDTight  = cms.string(EL_MVA_ID_TIGHT), 
#         mvaIsoTight = cms.string(EL_MVA_ISO_TIGHT), 
#     )
)

boostedMuons = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('muonPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('muonPFIsoChHad04')"),
        bdtIso = cms.string("userFloat('bdtisoNew')"),
        #bdtIso = cms.string("userFloat('bdtisonontrigDZ')"),
        #bdtIso = cms.string("userFloat('bdtisonontrigPFNOPU')"),
        sip    = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
#     flags = cms.PSet(
#         prlID = cms.string(MU_ID_PRL),
#         newID = cms.string(MU_ID_PF),
       #mvaID  = cms.string(MU_MVA_ID), # not for muons, at least not yet
#         mvaIso = cms.string(MU_MVA_ISO), 
#     )
)

reboosting = cms.Sequence(
    #updatedElectronRhoAA * 
    ( updatedElectronPFIsoChHad01 + updatedElectronPFIsoChHad02 + updatedElectronPFIsoChHad03 + updatedElectronPFIsoChHad04 + updatedElectronPFIsoChHad05) *
    boostedElectronsUpdatedPFIso * 
    boostedElectronsEAPFIso   * boostedElectronsFixBDTID * boostedElectrons +
    boostedMuonsUpdatedBDTIso * 
    (updatedMuonPFIsoNHad04 + updatedMuonPFIsoPhoton04) *
    boostedMuonsUpdatedPFIso *
    boostedMuonsEAPFIso * boostedMuons 
)


