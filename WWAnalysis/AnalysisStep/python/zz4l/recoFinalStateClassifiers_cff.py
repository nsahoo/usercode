import FWCore.ParameterSet.Config as cms

z1EE = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("bestZ"),
    cut = cms.string("abs(daughter(0).pdgId) == 11"),
    filter = cms.bool(True),
)
z1MM = z1EE.clone(
    cut = "abs(daughter(0).pdgId) == 13"
)

countGoodLep4M = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodMu"),
    minNumber = cms.uint32(4),
)
countGoodLep4E = countGoodLep4M.clone(
    src = cms.InputTag("goodEl"),
)
countGoodLep2M = countGoodLep4M.clone(minNumber = 2)
countGoodLep2E = countGoodLep4E.clone(minNumber = 2)
countGoodLep2E2M = cms.Sequence(countGoodLep2E + countGoodLep2M)

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

def makeSplittedPaths4L(process, name):
    from PhysicsTools.PatAlgos.tools.helpers import listModules
    path = getattr(process, name)
    modules = listModules(path)
    seq4E   = cms.Sequence(sum(modules[1:], modules[0]))
    seq4M   = cms.Sequence(sum(modules[1:], modules[0]))
    seq2E2M = cms.Sequence(sum(modules[1:], modules[0]))
    if process.bestZ in modules:
        seq4E.replace(process.bestZ, process.bestZ + process.z1EE)
        seq4M.replace(process.bestZ, process.bestZ + process.z1MM)
    if process.countGoodLep in modules:
        seq4E.replace(process.countGoodLep, process.countGoodLep4E)
        seq4M.replace(process.countGoodLep, process.countGoodLep4M)
        seq2E2M.replace(process.countGoodLep, process.countGoodLep2E2M)
    if process.best4L in modules:
        seq4E.replace(process.best4L, process.best4L + process.zz4E)
        seq4M.replace(process.best4L, process.best4L + process.zz4M)
        seq2E2M.replace(process.best4L, process.best4L + process.zz2E2M)
    setattr(process, name+"_4E", cms.Path(seq4E))
    setattr(process, name+"_4M", cms.Path(seq4M))
    setattr(process, name+"_2E2M", cms.Path(seq2E2M))
