import FWCore.ParameterSet.Config as cms

z1EE = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("bestZ"),
    cut = cms.string("abs(daughter(0).pdgId) == 11"),
    filter = cms.bool(True),
)
z1MM = z1EE.clone(
    cut = "abs(daughter(0).pdgId) == 13"
)

# for 3-path logic
oneZEE = cms.EDFilter("SkimEvent2LSelector",
    src = cms.InputTag("zll"),
    cut = cms.string("abs(daughter(0).pdgId) == 11"),
    filter = cms.bool(True),
)
oneZMM = oneZEE.clone(cut = "abs(daughter(0).pdgId) == 13")
bestZMM = cms.EDFilter("SkimEvent2LSorter",
    src = cms.InputTag("oneZMM"),
    sortBy = cms.string("abs(mass - 91.188)"),
    sortAscending = cms.bool(True),
    maxNumber = cms.uint32(1),
)
bestZEE = bestZMM.clone(src = "oneZEE")
selectedZEE = oneZEE.clone(src = "bestZEE", cut = "40 < mass < 120")
selectedZMM = oneZMM.clone(src = "bestZMM", cut = "40 < mass < 120")


countGoodLep4M = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodMu"),
    minNumber = cms.uint32(4),
)
countGoodLep4E = countGoodLep4M.clone(
    src = cms.InputTag("goodEl"),
)
countGoodLep2M = countGoodLep4M.clone(minNumber = 2)
countGoodLep3M = countGoodLep4M.clone(minNumber = 3)
countGoodLep2E = countGoodLep4E.clone(minNumber = 2)
countGoodLep3E = countGoodLep4E.clone(minNumber = 2)
countGoodLep2E2M = cms.Sequence(countGoodLep2E + countGoodLep2M)
countLooseLep2M = countGoodLep2M.clone(src = 'looseMu')
countLooseLep2E = countGoodLep2E.clone(src = 'looseEl')
countLooseLep3M = countGoodLep3M.clone(src = 'looseMu')
countLooseLep3E = countGoodLep3E.clone(src = 'looseEl')
countLooseLep4M = countGoodLep4M.clone(src = 'looseMu')
countLooseLep4E = countGoodLep4E.clone(src = 'looseEl')
countLooseLep2E2M = cms.Sequence(countLooseLep2E + countLooseLep2M)

zz4E = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("best4L"),
    cut = cms.string("abs(daughter(0).daughter(0).pdgId) == 11 && abs(daughter(1).daughter(0).pdgId) == 11"),
    filter = cms.bool(True),
)
zz4M = zz4E.clone(
    cut = cms.string("abs(daughter(0).daughter(0).pdgId) == 13 && abs(daughter(1).daughter(0).pdgId) == 13"),
)
zz2E2M = zz4E.clone(
    cut = cms.string("abs(daughter(0).daughter(0).pdgId) != abs(daughter(1).daughter(0).pdgId)"),
)

def makeSplittedPaths4L(process, name, triggerFilterName=None):
    from PhysicsTools.PatAlgos.tools.helpers import listModules
    path = getattr(process, name)
    modules = listModules(path)
    seq4E   = cms.Sequence(sum(modules[1:], modules[0]))
    seq4M   = cms.Sequence(sum(modules[1:], modules[0]))
    seq2E2M = cms.Sequence(sum(modules[1:], modules[0]))
    if process.bestZ in modules:
        seq4E.replace(process.bestZ, process.bestZ + process.z1EE)
        seq4M.replace(process.bestZ, process.bestZ + process.z1MM)
    if process.countGoodLep3 in modules:
        seq4E.replace(process.countGoodLep3, process.countGoodLep3E)
        seq4M.replace(process.countGoodLep3, process.countGoodLep3M)
    if process.countGoodLep in modules:
        seq4E.replace(process.countGoodLep, process.countGoodLep4E)
        seq4M.replace(process.countGoodLep, process.countGoodLep4M)
        seq2E2M.replace(process.countGoodLep, process.countGoodLep2E2M)
    if process.countLooseLep3 in modules:
        seq4E.replace(process.countLooseLep3, process.countLooseLep3E)
        seq4M.replace(process.countLooseLep3, process.countLooseLep3M)
    if process.countLooseLep in modules:
        seq4E.replace(process.countLooseLep, process.countLooseLep4E)
        seq4M.replace(process.countLooseLep, process.countLooseLep4M)
        seq2E2M.replace(process.countLooseLep, process.countLooseLep2E2M)
    if process.best4L in modules:
        seq4E.replace(process.best4L, process.best4L + process.zz4E)
        seq4M.replace(process.best4L, process.best4L + process.zz4M)
        seq2E2M.replace(process.best4L, process.best4L + process.zz2E2M)
    if triggerFilterName:
        seq4E.replace(getattr(process,triggerFilterName), getattr(process,triggerFilterName+"_El"))
        seq4M.replace(getattr(process,triggerFilterName), getattr(process,triggerFilterName+"_Mu"))
    setattr(process, name+"_4E", cms.Path(seq4E))
    setattr(process, name+"_4M", cms.Path(seq4M))
    setattr(process, name+"_2E2M", cms.Path(seq2E2M))
    ## the three-path logic should come here, but not immediately
    modules4E = listModules(getattr(process, name+"_4E"))
    modules4M = listModules(getattr(process, name+"_4M"))
    seq4E_3P = cms.Sequence(sum(modules4E[1:], modules4E[0]))
    seq4M_3P = cms.Sequence(sum(modules4M[1:], modules4M[0]))
    if process.oneZ in modules:
        seq4E_3P.replace(process.oneZ, process.oneZEE + process.bestZEE + process.selectedZEE)
        seq4M_3P.replace(process.oneZ, process.oneZMM + process.bestZMM + process.selectedZMM)
    setattr(process, name+"_4E_3Path", cms.Path(seq4E_3P))
    setattr(process, name+"_4M_3Path", cms.Path(seq4M_3P))
