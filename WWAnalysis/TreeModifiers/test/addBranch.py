import FWCore.ParameterSet.Config as cms
import math

process = cms.Process("BranchAdded")

process.General = cms.PSet(   
    inputFiles    =     cms.vstring('file:hzzTree_id1115.root','file:hzzTree_id1117.root'),
    branchAdders    =     cms.vstring('BaseWeightBranchAdder',
                                      'PileupBranchAdder',
                                      'LeptSfBranchAdder'),


)


process.BaseWeightBranchAdder = cms.PSet(
    #treesToProcess = cms.vstring('zz4lTree','muonTree'),
    treesToProcess = cms.vstring(), # if empty process all the trees in the file
    newBranchName = cms.string("baseW"),
    datasetConf = cms.string("8TeV"),
)


process.PileupBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zz4lTree'),
    newBranchName = cms.string("puW"),
    pileupConf = cms.string("2012"),
)

process.LeptSfBranchAdder = cms.PSet(
    treesToProcess = cms.vstring('zz4lTree'),
    newBranchName = cms.string("lepW"),
    datasetConf = cms.string("2012"),
)



### HERE YOU HAVE TO ADD THE CONFIGURATION FOR ANY NEW BRANCHADDER
#process.myBranchAdder = cms.PSet(
#    treesToProcess = cms.vstring('zz4lTree'),
#    newBranchName = cms.string("puW"),
#)
