import FWCore.ParameterSet.Config as cms
import math

from WWAnalysis.AnalysisStep.cutMassDependPSets_cfi import *


defaultWW = cms.VPSet(
    cms.PSet(
        label = cms.string("10/10"),
        cut = cms.string("(1)"),
    ),
    cms.PSet(
        label = cms.string("20/10"),
        cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMin() > 10 && ptMax() > 20")
    ),
    cms.PSet(
        label = cms.string("Lepton Veto"),
        cut = cms.string("nExtraLep(10) == 0")
    ),
    cms.PSet(
        label = cms.string("Trigger"),
        cut = cms.string("triggerMatchingCut('MC')")
    ),
    cms.PSet(
        label = cms.string("#slash{E}_{T}"),
        cut = cms.string("pfMet() > 20")
    ),
    cms.PSet(
        label = cms.string("m_{ll} > 12"),
        cut = cms.string("mll() > 12")
    ),
    cms.PSet(
        label = cms.string("Z-veto"),
        cut = cms.string("abs(mll()-91.1876)>15.0") 
    ),
    cms.PSet(
        label = cms.string("minProj#slash{E}_{T}"),
        cut = cms.string("min(projPfMet(),projChargedMet()) > 35")
    ),
    cms.PSet(
        label = cms.string("Jet Veto"),
        cut = cms.string("nCentralJets( 30.0, 5.0, 1 ) == 0"),
    ),
    cms.PSet(
        label = cms.string("Soft #mu Veto"),
        cut = cms.string("nSoftMu( 3.0 ) == 0"),
    ),
    cms.PSet(
        label = cms.string("b-tag Veto"),
        cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) == 0")
    ),
)

ttBar = cms.VPSet(
    cms.PSet(
        label = cms.string("10/10"),
        cut = cms.string("(1)"),
    ),
    cms.PSet(
        label = cms.string("20/10"),
        cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMin() > 10 && ptMax() > 20")
    ),
    cms.PSet(
        label = cms.string("Lepton Veto"),
        cut = cms.string("nExtraLep(10) == 0")
    ),
    cms.PSet(
        label = cms.string("Trigger"),
        cut = cms.string("triggerMatchingCut('MC')")
    ),
    cms.PSet(
        label = cms.string("#slash{E}_{T}"),
        cut = cms.string("pfMet() > 20")
    ),
    cms.PSet(
        label = cms.string("m_{ll} > 12"),
        cut = cms.string("mll() > 12")
    ),
    cms.PSet(
        label = cms.string("Z-veto"),
        cut = cms.string("abs(mll()-91.1876)>15.0") 
    ),
    cms.PSet(
        label = cms.string("pMin#slash{E}_{T}"),
        cut = cms.string("min(projPfMet(),projChargedMet()) > 35")
    ),
)

highPtB = cms.PSet(
    label = cms.string("b-Jet"),
    cut = cms.string("bTaggedJetsOver( 30.0, 2.1 ) == 1 && nCentralJets( 30.0, 5.0, 1) == 1"),
)

softMu = cms.PSet(
    label = cms.string("Soft #mu Veto"),
    cut = cms.string("nSoftMu( 3.0 ) > 0"),
)

bTag = cms.PSet(
    label = cms.string("b-tag Veto"),
    cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) > 0")
)

NoBtag = cms.PSet(
    label = cms.string("b-tag Veto"),
    cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) == 0")
)

topTag = cms.PSet(
    label = cms.string("Top Veto"),
    cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) > 0 || nSoftMu( 3.0 ) > 0")
)


def cloneVPSet(oldvpset):
    newvpset = cms.VPSet()
    for pset in oldvpset:
        newvpset.append(pset.clone())
    return newvpset

def switchToOppoFlavor(vp):
    vp[6].cut = "(1)"
    vp[7].cut = "min(projPfMet(),projChargedMet()) > 20.0"

def switchToTrailingElectron(vp):
    vp[1].cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMin() > 15 && ptMax() > 20")

#swaps for PDs
MC =             "MC"
SingleMuon =     "SingleMuon"
DoubleMuon =     "DoubleMuon"
MuEG =           "MuEG"
DoubleElectron = "DoubleElectron"

def switchTrigger(vp,value):
    vp[3].cut = "triggerMatchingCut('{0}')".format(value)

def addMassDependentCuts(vpa,vpb):
    for pset in vpb:
        vpa.append(pset)




