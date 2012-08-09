import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("BranchAdded")

process.General = cms.PSet(   
    inputFiles    =     cms.vstring('file:hzzTree_id2600.root'),
    branchAdders    =     cms.vstring('BaseWeightBranchAdder',
                                      'PileupBranchAdder',
                                      'LeptSfBranchAdder',
                                      'LeptPRBranchAdder',
                                      ),


)


process.BaseWeightBranchAdder = cms.PSet(
    #treesToProcess = cms.vstring('zz4lTree','muonTree'),
    treesToProcess = cms.vstring(), # if empty process all the trees in the file
    newBranchName = cms.string("baseW"),
    datasetConf = cms.string("8TeV"),
)


process.PileupBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zz4lTree','zxTree','anyZxTree',
                                 'zmmtree','zeetree',
                                 'zllmtree','zlletree'),
    newBranchName = cms.string("puW"),
    pileupConf = cms.string("2012"),
)

process.LeptSfBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zz4lTree','zxTree','anyZxTree'),
    newBranchName = cms.string("lepW"),
    datasetConf = cms.string("2012"),
)

process.LeptPRBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zxTree'),
    newBranchName = cms.string("lepPRW"),
    datasetConf = cms.string("2012"),
)



### HERE YOU HAVE TO ADD THE CONFIGURATION FOR ANY NEW BRANCHADDER
#process.myBranchAdder = cms.PSet(
#    treesToProcess = cms.vstring('zz4lTree'),
#    newBranchName = cms.string("puW"),
#)

#process.General.inputFiles = cms.vstring(
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id100.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id101.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id102.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id103.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id104.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id105.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id106.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id107.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id111.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1115.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1117.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1119.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1120.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1121.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1123.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1124.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1125.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1126.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1127.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id114.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id113.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1145.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id115.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1150.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id118.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id117.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id116.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id119.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1180.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id120.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1200.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id121.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1220.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1250.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1300.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1325.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1350.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1400.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1450.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1500.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1550.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id1600.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2115.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2117.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2119.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2120.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2121.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2123.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2124.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2125.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2127.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2126.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2130.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2150.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2180.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2200.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2220.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2250.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2300.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2400.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2450.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2500.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2550.root",
#"file:/data/hzz4l/step3/2012/MC/HZZ4L_S3_V00/hzzTree_id2600.root",
#)
