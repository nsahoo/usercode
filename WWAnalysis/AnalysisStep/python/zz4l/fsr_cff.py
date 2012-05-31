import FWCore.ParameterSet.Config as cms

fsrPhotonsNoEleSC = cms.EDProducer("PATPFParticleCleaner",
    src = cms.InputTag("boostedFsrPhotons"),
    preselection = cms.string("pt > 2 && abs(eta) < 2.4"),
    checkOverlaps = cms.PSet(
        eleVeto = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string("abs(pdgId) == 11"),
            deltaR = cms.double(2.05), # wider than sqrt(2^2 + 0.05*2)
            pairCut = cms.string("abs(deltaPhi(ele.phi,pf.phi)) < 2.0 && abs(pf.eta - ele.eta) < 0.05 || deltaR(ele.eta,ele.phi,pf.eta,pf.phi) < 0.15"),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False), ## Will Veto these photons later
        ),
    ),
    finalCut = cms.string(""),
)

## Patrick's Isolation tune
_netural_iso_pj = "pf.userFloat('fsrPhotonPFIsoPhoton03') + pf.userFloat('fsrPhotonPFIsoNHad03')"
_neutral_iso_dbeta_pj = "max(0, %s - 0.5*pf.userFloat('fsrPhotonPFIsoChHadPU03pt02'))" % _netural_iso_pj
_comb_rel_iso_dbeta_pj = "(pf.userFloat('fsrPhotonPFIsoChHad03pt02') + %s)/pf.pt" % _neutral_iso_dbeta_pj
## Muon Isolation tune
_netural_iso_04 = "pf.userFloat('fsrPhotonPFIsoPhoton04') + pf.userFloat('fsrPhotonPFIsoNHad04')"
_neutral_iso_dbeta_04 = "max(0, %s - 0.5*pf.userFloat('fsrPhotonPFIsoChHadPU04'))" % _netural_iso_04
_comb_rel_iso_dbeta_04 = "(pf.userFloat('fsrPhotonPFIsoChHad04') + %s)/pf.pt" % _neutral_iso_dbeta_04
## DeltaR 03 using deltaBeta from 0.4 with AA = 9.6/16 (determined looking at some events in MC)
## Should reproduce farily accurately the deltaBeta correction for (R < 0.3, pt > 0.2)
_netural_iso_03aa = "pf.userFloat('fsrPhotonPFIsoPhoton03') + pf.userFloat('fsrPhotonPFIsoNHad03')"
_neutral_iso_dbeta_03aa = "max(0, %s - %f*0.5*pf.userFloat('fsrPhotonPFIsoChHadPU03pt02'))" % (_netural_iso_03aa, 9.7/16)
_comb_rel_iso_dbeta_03aa = "(pf.userFloat('fsrPhotonPFIsoChHad03pt02') + %s)/pf.pt" % _neutral_iso_dbeta_03aa

fsrPhotons = cms.EDProducer("PATPFParticleCleaner",
    src = cms.InputTag("fsrPhotonsNoEleSC"),
    preselection = cms.string(""),
    checkOverlaps = cms.PSet(
        goodLepNoIso = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            pairCut = cms.string("deltaR < 0.07 || pf.pt > 4"),
            deltaR = cms.double(0.5),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False),
        ),
        goodLepPJ = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            pairCut = cms.string("deltaR < 0.07 || (pf.pt > 4 && %s < 1.0)" % _comb_rel_iso_dbeta_pj),
            deltaR = cms.double(0.5),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False),
        ),
        goodLep04 = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            pairCut = cms.string("deltaR < 0.07 || (pf.pt > 4 && %s < 1.0)" % _comb_rel_iso_dbeta_04),
            deltaR = cms.double(0.5),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False),
        ),
        goodLep03aa = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            pairCut = cms.string("deltaR < 0.07 || (pf.pt > 4 && %s < 1.0)" % _comb_rel_iso_dbeta_03aa),
            deltaR = cms.double(0.5),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False),
        ),
    ),
    finalCut = cms.string(""),
)

fsrPhotonSeq = cms.Sequence(fsrPhotonsNoEleSC * fsrPhotons)
