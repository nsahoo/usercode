import FWCore.ParameterSet.Config as cms

hThursday120 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 40"),      cut = cms.string("mll()    < 40.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 20"), cut = cms.string("ptMax()  > 20.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 10"), cut = cms.string("ptMin()  > 10.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 115"), cut = cms.string("dPhillInDegrees() < 115.0")),
    cms.PSet( label = cms.string("70 < m_{T} < 120"), cut = cms.string("mTHiggs('PFMET') > 70 && mTHiggs('PFMET') < 120")),
)

hThursday130 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 45"),      cut = cms.string("mll()    < 45.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 25"), cut = cms.string("ptMax()  > 25.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 10"), cut = cms.string("ptMin()  > 10.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("75 < m_{T} < 125"), cut = cms.string("mTHiggs('PFMET') > 75 && mTHiggs('PFMET') < 125")),
)

hThursday140 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 45"),      cut = cms.string("mll()    < 45.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 25"), cut = cms.string("ptMax()  > 25.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 15"), cut = cms.string("ptMin()  > 15.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("80 < m_{T} < 130"), cut = cms.string("mTHiggs('PFMET') > 80 && mTHiggs('PFMET') < 130")),
)

hThursday150 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 27"), cut = cms.string("ptMax()  > 27.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("80 < m_{T} < 150"), cut = cms.string("mTHiggs('PFMET') > 80 && mTHiggs('PFMET') < 150")),
)

hThursday160 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 30"), cut = cms.string("ptMax()  > 30.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 60"),  cut = cms.string("dPhillInDegrees() < 60.0")),
    cms.PSet( label = cms.string("90 < m_{T} < 160"), cut = cms.string("mTHiggs('PFMET') > 90 && mTHiggs('PFMET') < 160")),
)

hThursday170 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 34"), cut = cms.string("ptMax()  > 34.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 60"),  cut = cms.string("dPhillInDegrees() < 60.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 170"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 170")),
)

hThursday180 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 60"),      cut = cms.string("mll()    < 60.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 36"), cut = cms.string("ptMax()  > 36.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 70"),  cut = cms.string("dPhillInDegrees() < 70.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 180"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 180")),
)

hThursday190 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 80"),      cut = cms.string("mll()    < 80.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 38"), cut = cms.string("ptMax()  > 38.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 190"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 190")),
)

hThursday200 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 90"),      cut = cms.string("mll()    < 90.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 40"), cut = cms.string("ptMax()  > 40.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 100"),  cut = cms.string("dPhillInDegrees() < 100.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 200"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 200")),
)

hThursday250 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 150"),      cut = cms.string("mll()    < 150.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 55"), cut = cms.string("ptMax()  > 55.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 140"),  cut = cms.string("dPhillInDegrees() < 140.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 250"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 250")),
)

hThursday300 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 200"),      cut = cms.string("mll()    < 200.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 70"), cut = cms.string("ptMax()  > 70.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 300"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 300")),
)

hThursday350 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 250"),      cut = cms.string("mll()    < 250.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 80"), cut = cms.string("ptMax()  > 80.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 350"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 350")),
)

hThursday400 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 300"),      cut = cms.string("mll()    < 300.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 90"), cut = cms.string("ptMax()  > 90.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 400"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 400")),
)

hThursday450 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 350"),      cut = cms.string("mll()    < 350.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 110"), cut = cms.string("ptMax()  > 110.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 450"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 450")),
)

hThursday500 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 400"),      cut = cms.string("mll()    < 400.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 120"), cut = cms.string("ptMax()  > 120.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 500"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 500")),
)

hThursday550 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 450"),      cut = cms.string("mll()    < 450.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 130"), cut = cms.string("ptMax()  > 130.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 550"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 550")),
)

hThursday600 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 500"),      cut = cms.string("mll()    < 500.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 140"), cut = cms.string("ptMax()  > 140.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("20 < m_{T} < 600"), cut = cms.string("mTHiggs('PFMET') > 20 && mTHiggs('PFMET') < 600")),
)

# begin gammamrstar sets

hGammaMR120 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 40"),      cut = cms.string("mll()    < 40.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 20"), cut = cms.string("ptMax()  > 20.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 10"), cut = cms.string("ptMin()  > 10.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 115"), cut = cms.string("dPhillInDegrees() < 115.0")),
    cms.PSet( label = cms.string("0.5 < #gamma*m_{R}^{*}/m_{H} < 1.1"), cut = cms.string("gammaMRStar()/120 > 0.5 && gammaMRStar()/120 < 1.1")),
)

hGammaMR130 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 45"),      cut = cms.string("mll()    < 45.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 25"), cut = cms.string("ptMax()  > 25.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 10"), cut = cms.string("ptMin()  > 10.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("0.5 < #gamma*m_{R}^{*}/m_{H} < 1.1"), cut = cms.string("gammaMRStar()/130 > 0.5 && gammaMRStar()/130 < 1.1")),
)

hGammaMR140 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 45"),      cut = cms.string("mll()    < 45.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 25"), cut = cms.string("ptMax()  > 25.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 15"), cut = cms.string("ptMin()  > 15.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("0.5 < #gamma*m_{R}^{*}/m_{H} < 1.1"), cut = cms.string("gammaMRStar()/140 > 0.5 && gammaMRStar()/140 < 1.1")),
)

hGammaMR150 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 27"), cut = cms.string("ptMax()  > 27.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("0.75 < #gamma*m_{R}^{*}/m_{H} < 1.15"), cut = cms.string("gammaMRStar()/150 > 0.75 && gammaMRStar()/150 < 1.15")),
)

hGammaMR160 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 30"), cut = cms.string("ptMax()  > 30.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 60"),  cut = cms.string("dPhillInDegrees() < 60.0")),
    cms.PSet( label = cms.string("0.75 < #gamma*m_{R}^{*}/m_{H} < 1.2"), cut = cms.string("gammaMRStar()/160 > 0.75 && gammaMRStar()/160 < 1.2")),
)

hGammaMR170 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 50"),      cut = cms.string("mll()    < 50.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 34"), cut = cms.string("ptMax()  > 34.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 60"),  cut = cms.string("dPhillInDegrees() < 60.0")),
    cms.PSet( label = cms.string("0.75 < #gamma*m_{R}^{*}/m_{H} < 1.3"), cut = cms.string("gammaMRStar()/170 > 0.75 && gammaMRStar()/170 < 1.3")),
)

hGammaMR180 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 60"),      cut = cms.string("mll()    < 60.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 36"), cut = cms.string("ptMax()  > 36.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 70"),  cut = cms.string("dPhillInDegrees() < 70.0")),
    cms.PSet( label = cms.string("0.8 < #gamma*m_{R}^{*}/m_{H} < 1.4"), cut = cms.string("gammaMRStar()/180 > 0.8 && gammaMRStar()/180 < 1.4")),
)

hGammaMR190 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 80"),      cut = cms.string("mll()    < 80.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 38"), cut = cms.string("ptMax()  > 38.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 90"),  cut = cms.string("dPhillInDegrees() < 90.0")),
    cms.PSet( label = cms.string("0.8 < #gamma*m_{R}^{*}/m_{H} < 1.4"), cut = cms.string("gammaMRStar()/190 > 0.8 && gammaMRStar()/190 < 1.4")),
)

hGammaMR200 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 90"),      cut = cms.string("mll()    < 90.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 40"), cut = cms.string("ptMax()  > 40.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 100"),  cut = cms.string("dPhillInDegrees() < 100.0")),
    cms.PSet( label = cms.string("0.8 < #gamma*m_{R}^{*}/m_{H} < 1.4"), cut = cms.string("gammaMRStar()/200 > 0.8 && gammaMRStar()/200 < 1.4")),
)

hGammaMR250 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 150"),      cut = cms.string("mll()    < 150.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 55"), cut = cms.string("ptMax()  > 55.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 140"),  cut = cms.string("dPhillInDegrees() < 140.0")),
    cms.PSet( label = cms.string("0.9 < #gamma*m_{R}^{*}/m_{H} < 1.4"), cut = cms.string("gammaMRStar()/250 > 0.9 && gammaMRStar()/250 < 1.4")),
)

hGammaMR300 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 200"),      cut = cms.string("mll()    < 200.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 70"), cut = cms.string("ptMax()  > 70.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("0.9 < #gamma*m_{R}^{*}/m_{H} < 1.5"), cut = cms.string("gammaMRStar()/300 > 0.9 && gammaMRStar()/300 < 1.5")),
)

hGammaMR350 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 250"),      cut = cms.string("mll()    < 250.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 80"), cut = cms.string("ptMax()  > 80.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("0.9 < #gamma*m_{R}^{*}/m_{H} < 1.5"), cut = cms.string("gammaMRStar()/350 > 0.9 && gammaMRStar()/350 < 1.5")),
)

hGammaMR400 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 300"),      cut = cms.string("mll()    < 300.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 90"), cut = cms.string("ptMax()  > 90.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("0.9 < #gamma*m_{R}^{*}/m_{H} < 1.5"), cut = cms.string("gammaMRStar()/400 > 0.9 && gammaMRStar()/400 < 1.5")),
)

hGammaMR450 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 350"),      cut = cms.string("mll()    < 350.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 110"), cut = cms.string("ptMax()  > 110.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("0.9 < #gamma*m_{R}^{*}/m_{H} < 1.5"), cut = cms.string("gammaMRStar()/450 > 0.9 && gammaMRStar()/450 < 1.5")),
)

hGammaMR500 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 400"),      cut = cms.string("mll()    < 400.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 120"), cut = cms.string("ptMax()  > 120.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("0.8 < #gamma*m_{R}^{*}/m_{H} < 1.5"), cut = cms.string("gammaMRStar()/500 > 0.8 && gammaMRStar()/500 < 1.5")),
)

hGammaMR550 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 450"),      cut = cms.string("mll()    < 450.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 130"), cut = cms.string("ptMax()  > 130.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("0.8 < #gamma*m_{R}^{*}/m_{H} < 1.5"), cut = cms.string("gammaMRStar()/550 > 0.8 && gammaMRStar()/550 < 1.5")),
)

hGammaMR600 = cms.VPSet(
    cms.PSet( label = cms.string("m_{ll} < 500"),      cut = cms.string("mll()    < 500.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MAX} > 140"), cut = cms.string("ptMax()  > 140.0 ")),
    cms.PSet( label = cms.string("p_{T}^{MIN} > 25"), cut = cms.string("ptMin()  > 25.0 ")),
    cms.PSet( label = cms.string("#Delta#phi < 175"),  cut = cms.string("dPhillInDegrees() < 175.0")),
    cms.PSet( label = cms.string("0.8 < #gamma*m_{R}^{*}/m_{H} < 1.5"), cut = cms.string("gammaMRStar()/600 > 0.8 && gammaMRStar()/600 < 1.5")),
)
