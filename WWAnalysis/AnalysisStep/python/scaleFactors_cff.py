import FWCore.ParameterSet.Config as cms
from glob import glob


data2012A = {
    '100': ['SingleElectron2012A'],
    '101': ['SingleMuon2012A'],
    '102': ['DoubleElectron2012A'],
    '103': ['DoubleMuon2012A'],
    '104': ['MuEG2012A'],
}
data2012B = {
    '200': ['SingleElectron2012B'],
    '201': ['SingleMuon2012B'],
    '202': ['DoubleElectron2012B'],
    '203': ['DoubleMuon2012B'],
    '204': ['MuEG2012B'],
    '210': ['SingleElectron2012B'],
    '211': ['SingleMuon2012B'],
    '212': ['DoubleElectron2012B'],
    '213': ['DoubleMuon2012B'],
    '214': ['MuEG2012B'],
}

smallSamples={
        '070':['WZFull',0.0022456186],
        '077':['ZZTo2L2Nu',0.0007327844],
        '011':['TtWFullDR',0.0224598017],
        '012':['TbartWFullDR',0.0226508734],
        '036':['DY10toLLMad',0.1223652890],
        '000':['WWJets2LMad',0.0030065150],
        '001':['GluGluToWWTo4L',0.0016624874],
	'080':['WJetsToLNuMad',2.0392984539],
}

WJets={
        '080':['WJetsToLNuMad',2.0392984539],
}
TTJets={

        '010c1':['TTJetsMad',0.0348202597],
        '010c2':['TTJetsMad',0.0348202597],
}

ZJets={
        '037c0':['DY50toLLMad',0.1160585354],
        '037c1':['DY50toLLMad',0.1160585354],
        '037c2':['DY50toLLMad',0.1160585354],
        '037c3':['DY50toLLMad',0.1160585354],
        '037c4':['DY50toLLMad',0.1160585354],
        '037c5':['DY50toLLMad',0.1160585354],
        '037c6':['DY50toLLMad',0.1160585354],
        '037c7':['DY50toLLMad',0.1160585354],
        '037c8':['DY50toLLMad',0.1160585354],
        '037c9':['DY50toLLMad',0.1160585354],

}
mcSummer12={
        '070':['WZFull',0.0032327020],
        '077':['ZZTo2L2Nu',0.0007327844],
        '011':['TtWFullDR',0.0224598017],
        '012':['TbartWFullDR',0.0226508734],
#        '037':['DY50toLLMad',1.1935135135],
        '036':['DY10toLLMad',0.1223652890],
        '001':['GluGluToWWTo4L',0.0016624874],
#        '010':['TTJetsMad',0.0334311887],
        '000':['WWJets2LMad',0.0030065150],
}

hww={
	'1110':['ggToH110toWWTo2LAndTau2Nu',0.0004262067],
	'1115':['ggToH115toWWTo2LAndTau2Nu',0.0006966061],
	'1120':['ggToH120toWWTo2LAndTau2Nu',0.0011154932],
	'1125':['ggToH125toWWTo2LAndTau2Nu',0.0014792771],
	'1130':['ggToH130toWWTo2LAndTau2Nu',0.0019329875],
	'1135':['ggToH135toWWTo2LAndTau2Nu',0.0023721164],
	'1140':['ggToH140toWWTo2LAndTau2Nu',0.0027582844],
	'1145':['ggToH145toWWTo2LAndTau2Nu',0.0030737075],
	'1150':['ggToH150toWWTo2LAndTau2Nu',0.0033341662],
	'1155':['ggToH155toWWTo2LAndTau2Nu',0.0035425465],
	'1160':['ggToH160toWWTo2LAndTau2Nu',0.0037636891],
	'1170':['ggToH170toWWTo2LAndTau2Nu',0.0034241582],
	'1180':['ggToH180toWWTo2LAndTau2Nu',0.0028858716],
	'1190':['ggToH190toWWTo2LAndTau2Nu',0.0021391306],
	'1200':['ggToH200toWWTo2LAndTau2Nu',0.0018039530],
	'1250':['ggToH250toWWTo2LAndTau2Nu',0.0010968517],
	'1300':['ggToH300toWWTo2LAndTau2Nu',0.0008073467],
	'1350':['ggToH350toWWTo2LAndTau2Nu',0.0007480515],
	'1400':['ggToH400toWWTo2LAndTau2Nu',0.0005960259],
	'1450':['ggToH450toWWTo2LAndTau2Nu',0.0003802122],
	'1500':['ggToH500toWWTo2LAndTau2Nu',0.0002471469],
	'1550':['ggToH550toWWTo2LAndTau2Nu',0.0001567124],
	'1600':['ggToH600toWWTo2LAndTau2Nu',0.0001031688],
	'1700':['ggToH700toWWTo2LAndTau2Nu',0.0000462051],
	'1800':['ggToH800toWWTo2LAndTau2Nu',0.0000227653],
	'1900':['ggToH900toWWTo2LAndTau2Nu',0.0000121168],
	'2000':['ggToH1000toWWTo2LAndTau2Nu',0.0000068825],
}






