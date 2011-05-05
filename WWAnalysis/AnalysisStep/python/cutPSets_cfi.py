import FWCore.ParameterSet.Config as cms
import math

# hReOpt120 = cms.PSet( mllMaxFinal=cms.double(40 ), ptMaxFinal=cms.double(20 ), ptMinFinal=cms.double(10), deltaPhiLL=cms.double(115./180*math.pi))  
# hReOpt130 = cms.PSet( mllMaxFinal=cms.double(45 ), ptMaxFinal=cms.double(25 ), ptMinFinal=cms.double(10), deltaPhiLL=cms.double(90. /180*math.pi))  
# hReOpt140 = cms.PSet( mllMaxFinal=cms.double(45 ), ptMaxFinal=cms.double(25 ), ptMinFinal=cms.double(15), deltaPhiLL=cms.double(90. /180*math.pi))  
# hReOpt150 = cms.PSet( mllMaxFinal=cms.double(50 ), ptMaxFinal=cms.double(27 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(90. /180*math.pi))  
# hReOpt160 = cms.PSet( mllMaxFinal=cms.double(50 ), ptMaxFinal=cms.double(30 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(60. /180*math.pi))  
# hReOpt170 = cms.PSet( mllMaxFinal=cms.double(50 ), ptMaxFinal=cms.double(34 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(60. /180*math.pi))  
# hReOpt180 = cms.PSet( mllMaxFinal=cms.double(60 ), ptMaxFinal=cms.double(36 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(70. /180*math.pi))  
# hReOpt190 = cms.PSet( mllMaxFinal=cms.double(80 ), ptMaxFinal=cms.double(38 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(90. /180*math.pi))  
# hReOpt200 = cms.PSet( mllMaxFinal=cms.double(90 ), ptMaxFinal=cms.double(40 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(100./180*math.pi))  
# hReOpt250 = cms.PSet( mllMaxFinal=cms.double(150), ptMaxFinal=cms.double(55 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(140./180*math.pi))    
# hReOpt300 = cms.PSet( mllMaxFinal=cms.double(200), ptMaxFinal=cms.double(70 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(175./180*math.pi))  
# hReOpt350 = cms.PSet( mllMaxFinal=cms.double(250), ptMaxFinal=cms.double(80 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(175./180*math.pi))    
# hReOpt400 = cms.PSet( mllMaxFinal=cms.double(300), ptMaxFinal=cms.double(90 ), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(175./180*math.pi))  
# hReOpt450 = cms.PSet( mllMaxFinal=cms.double(350), ptMaxFinal=cms.double(110), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(175./180*math.pi))  
# hReOpt500 = cms.PSet( mllMaxFinal=cms.double(400), ptMaxFinal=cms.double(120), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(175./180*math.pi))  
# hReOpt550 = cms.PSet( mllMaxFinal=cms.double(450), ptMaxFinal=cms.double(130), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(175./180*math.pi))  
# hReOpt600 = cms.PSet( mllMaxFinal=cms.double(500), ptMaxFinal=cms.double(140), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(175./180*math.pi))  
# 
# hReOpt210 = cms.PSet( mllMaxFinal=cms.double(110), ptMaxFinal=cms.double( 44), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(110./180*math.pi))  
# hReOpt220 = cms.PSet( mllMaxFinal=cms.double(120), ptMaxFinal=cms.double( 48), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(120./180*math.pi))  
# hReOpt230 = cms.PSet( mllMaxFinal=cms.double(130), ptMaxFinal=cms.double( 52), ptMinFinal=cms.double(25), deltaPhiLL=cms.double(130./180*math.pi))  

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
        label = cms.string("p#slash{E}_{T}"),
        cut = cms.string("projPfMet() > 35")
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
        label = cms.string("p#slash{E}_{T}"),
        cut = cms.string("projPfMet() > 35")
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

topTag = cms.PSet(
    label = cms.string("Top Veto"),
    cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) > 0 || nSoftMu( 3.0 ) > 0")
)


hReOpt160 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 30"), cut = cms.string("ptMax()  > 30.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 60"),  cut = cms.string("dPhill() < 60.0 /180.0*3.14159265")),
)

def cloneVPSet(oldvpset):
    newvpset = cms.VPSet()
    for pset in oldvpset:
        newvpset.append(pset.clone())
    return newvpset

def switchToOppoFlavor(vp):
    vp[6].cut = "(1)"
    vp[7].cut = "(1)"

#swaps for PDs
MC =             "MC"
SingleMuon =     "SingleMuo"
DoubleMuon =     "DoubleMuo"
MuEG =           "MuE"
DoubleElectron = "DoubleElectro"

def switchTrigger(vp,value):
    vp[3].cut = "triggerMatchingCut('{0}')".format(value)

def addMassDependentCuts(vpa,vpb):
    for pset in vpb:
        vpa.append(pset)




