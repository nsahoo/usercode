import FWCore.ParameterSet.Config as cms

hReOpt120 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 45"),      cut = cms.string("mll()    < 45.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 20"), cut = cms.string("ptMax()  > 20.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 10"), cut = cms.string("ptMin()  > 10.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 115"),  cut = cms.string("dPhillInDegrees() < 115.0")),
)

hReOpt130 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 45"),      cut = cms.string("mll()    < 45.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 25"), cut = cms.string("ptMax()  > 25.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 10"), cut = cms.string("ptMin()  > 10.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
)

hReOpt140 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 45"),      cut = cms.string("mll()    < 45.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 25"), cut = cms.string("ptMax()  > 25.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 15"), cut = cms.string("ptMin()  > 15.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
)

hReOpt150 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 27"), cut = cms.string("ptMax()  > 27.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
)


hReOpt160 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 30"), cut = cms.string("ptMax()  > 30.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 60"),  cut = cms.string("dPhillInDegrees() < 60.0")),
)

hReOpt170 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 34"), cut = cms.string("ptMax()  > 34.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 60"),  cut = cms.string("dPhillInDegrees() < 60.0")),
)

hReOpt180 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 60"),      cut = cms.string("mll()    < 60.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 36"), cut = cms.string("ptMax()  > 36.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 70"),  cut = cms.string("dPhillInDegrees() < 70.0")),
)

hReOpt190 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 80"),      cut = cms.string("mll()    < 80.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 38"), cut = cms.string("ptMax()  > 38.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
)

hReOpt200 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 90"),      cut = cms.string("mll()    < 90.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 40"), cut = cms.string("ptMax()  > 40.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 100"),  cut = cms.string("dPhillInDegrees() < 100.0")),
)

hReOpt250 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 150"),      cut = cms.string("mll()    < 150.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 55"), cut = cms.string("ptMax()  > 55.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 140"),  cut = cms.string("dPhillInDegrees() < 140.0")),
)

hReOpt300 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 200"),      cut = cms.string("mll()    < 200.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 70"), cut = cms.string("ptMax()  > 70.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
)

hReOpt350 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 250"),      cut = cms.string("mll()    < 250.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 80"), cut = cms.string("ptMax()  > 80.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
)

hReOpt400 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 300"),      cut = cms.string("mll()    < 300.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 90"), cut = cms.string("ptMax()  > 90.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
)

hReOpt450 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 350"),      cut = cms.string("mll()    < 350.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 110"), cut = cms.string("ptMax()  > 110.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
)

hReOpt500 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 400"),      cut = cms.string("mll()    < 400.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 120"), cut = cms.string("ptMax()  > 120.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
)

hReOpt550 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 450"),      cut = cms.string("mll()    < 450.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 130"), cut = cms.string("ptMax()  > 130.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
)

hReOpt600 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 500"),      cut = cms.string("mll()    < 500.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 140"), cut = cms.string("ptMax()  > 140.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
)
