import FWCore.ParameterSet.Config as cms


#Weighting for the Winter10 Flat samples
winter10 = [ 0.0441997   , 0.051395    , 0.0487518   , 0.0478708   , 
             0.0508076   , 0.0499266   , 0.0456681   , 0.0508076   , 
             0.0553598   , 0.0487518   , 0.0499266   , 0.0505139   , 
             0.0465492   , 0.0461087   , 0.0478708   , 0.0483113   , 
             0.0424376   , 0.0437592   , 0.0453744   , 0.0368576   , 
             0.0267254   , 0.0132159   , 0.00602056  , 0.00249633  , 
             0.000293686 ]

newHildrethUnNorm = [ 3.14153e+06, 3.13983e+06, 3.13893e+06, 3.13798e+06,
                      3.12936e+06, 3.11072e+06, 3.07349e+06, 2.99678e+06,
                      2.88199e+06, 2.7231e+06 , 2.51695e+06, 2.27149e+06,
                      1.99438e+06, 1.71016e+06, 1.42121e+06, 1.1538e+06 ,
                      913435     , 707974     , 534684     , 394500     ,
                      285410     , 202605     , 141245     , 96448      ,
                      64853      , 42721      , 27673      , 17805      ,
                      11306      , 7027       , 4311       , 2658       ,
                      1458       , 934        , 572        ]

sum = 0
for x in newHildrethUnNorm:
    sum += x

newHildreth = [ x/sum for x in newHildrethUnNorm ][:25]


# Fall 11 values directly from the correct mixing module
# from SimGeneral.MixingModule.mix_E7TeV_FlatDist10_2011EarlyData_50ns_cfi import mix as mix11
# spring11 = [ x for x in mix11.input.nbPileupEvents.probValue ]

spring11 = [   0.069286816, 0.069286816, 0.069286816, 0.069286816,
             0.069286816, 0.069286816, 0.069286816, 0.069286816,
             0.069286816, 0.069286816, 0.069286816, 0.06518604 ,
             0.053861878, 0.040782032, 0.030135062, 0.019550796,
             0.012264707, 0.007449117, 0.004502075, 0.002194605,
             0.001166276, 0.000476543, 0.000188109, 7.52436E-05,
             1.25406E-05 
]


# Summer 11 PU_S4 recommended:
summer11 = [ 0.104109, 0.0703573, 0.0698445, 0.0698254,
             0.0697054, 0.0697907, 0.0696751, 0.0694486,
             0.0680332, 0.0651044, 0.0598036, 0.0527395,
             0.0439513, 0.0352202, 0.0266714, 0.019411,
             0.0133974, 0.00898536, 0.0057516, 0.00351493,
             0.00212087, 0.00122891, 0.00070592, 0.000384744,
             0.000219377 
]

puS3 = spring11[:]
puS4 = summer11[:]
mcNominal = summer11[:]



#poisson around 5
weight05 = [      0.0336897   , 0.0842243   , 0.140374    , 0.175467    , 
                  0.175467    , 0.146223    , 0.104445    , 0.065278    , 
                  0.0362656   , 0.0181328   , 0.00824218  , 0.00343424  , 
                  0.00132086  , 0.000471736 , 0.000157245 , 4.91392e-05 , 
                  1.44527e-05 , 4.01464e-06 , 1.05648e-06 , 2.64121e-07 , 
                  6.2886e-08  , 1.42923e-08 , 3.10701e-09 , 6.47295e-10 ,
                  1.585e-10 ]
                  

#poisson around 10
weight10 = [     0.000453999  , 0.00227      , 0.00756665   , 0.0189166    , 
                 0.0378333    , 0.0630555    , 0.0900792    , 0.112599     , 
                 0.12511      , 0.12511      , 0.113736     , 0.0947803    , 
                 0.0729079    , 0.0520771    , 0.0347181    , 0.0216988    , 
                 0.012764     , 0.00709111   , 0.00373216   , 0.00186608   , 
                 0.00088861   , 0.000403914  , 0.000175615  , 2.9269e-05   ,
                 4.6975e-05]
                 

#poisson around 15
weight15 = [     4.58853e-06  , 3.4414e-05   , 0.00017207   , 0.000645263  , 
                 0.00193579   , 0.00483947   , 0.0103703    , 0.0194443    , 
                 0.0324072    , 0.0486108    , 0.0662874    , 0.0828592    , 
                 0.0956068    , 0.102436     , 0.102436     , 0.0960336    , 
                 0.0847356    , 0.070613     , 0.0557471    , 0.0418103    , 
                 0.0298645    , 0.0203622    , 0.0132797    , 0.011171     ,
                 4.979876e-3]

# kinda similar to the distro in 2010A ...
weight2010 = [   0.145168     , 0.251419     , 0.251596     , 0.17943      , 
                 0.10         , 0.05         , 0.02         , 0.01         , 
                 0.005        , 0.002        , 0.001        , 0            ,
                 0            , 0            , 0            , 0            ,
                 0            , 0            , 0            , 0            ,
                 0            , 0            , 0            , 0            ,
                 0 ]


from WWAnalysis.Misc.certifiedPileUp_cfi import pu2011A,pu2011B,pu2011AB,puS4fromMC
# mcTemplate = puS4fromMC[:]

# mcWeights   = puS4fromMC
# dataWeights = puVectorCertified[:]

# currentNew = puVectorCertified[:]
#FIXME
# currentOld = puVectorCertified[:]
# lpNew = puVectorCertified[:]
# lpOld = puVectorCertified[:]

# to be used on our MC
# reWeightVector =  [ dataWeights[i]/mcWeights[i] for i in range(len(mcWeights)) ] 

mcWeights = puS4fromMC[:25]

puS4fromMC = puS4fromMC[:35]
pu2011A    = pu2011A[:35]
w2011A     =  [ pu2011A[i]/puS4fromMC[i] for i in range(35) ] 
pu2011B    = pu2011B[:35]
w2011B     =  [ pu2011B[i]/puS4fromMC[i] for i in range(35) ] 
pu2011AB   = pu2011AB[:35]
w2011AB    =  [ pu2011AB[i]/puS4fromMC[i] for i in range(35) ] 

#reweighting for 05
vtx05 = weight05[:]
for i in range(0,len(weight05)):
    vtx05[i] = vtx05[i]/mcWeights[i]

#reweighting for 10
vtx10 = weight10[:]
for i in range(0,len(weight10)):
    vtx10[i] = vtx10[i]/mcWeights[i]

#reweighting for 15
vtx15 = weight15[:]
for i in range(0,len(weight15)):
    vtx15[i] = vtx15[i]/mcWeights[i]

#reweighting for 2010
vtx2010 = weight2010[:]
for i in range(0,len(weight2010)):
    vtx2010[i] = vtx2010[i]/mcWeights[i]

