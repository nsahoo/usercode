import FWCore.ParameterSet.Config as cms
from glob import glob

data2012A = {
    '100': ['SingleElectron2012A'],
    '101': ['SingleMuon2012A'],
    '102': ['DoubleElectron2012A'],
    '103': ['DoubleMuon2012A'],
    '104': ['MuEG2012A'],
    '110': ['SingleElectron2012A'],
    '111': ['SingleMuon2012A'],
    '112': ['DoubleElectron2012A'],
    '113': ['DoubleMuon2012A'],
    '114': ['MuEG2012A']
}
data2012B = {
    '200': ['SingleElectron2012B'],
    '201': ['SingleMuon2012B'],
    '202': ['DoubleElectron2012B'],
    '203': ['DoubleMuon2012B'],
    '204': ['MuEG2012B'],
}

data2012C = {
    '305': ['SingleElectron2012C'],
    '306': ['SingleMuon2012C'],
    '307': ['DoubleElectron2012C'],
    '308': ['DoubleMuon2012C'],
    '309': ['MuEG2012C'],
    '310': ['SingleElectron2012C'],
    '311': ['SingleMuon2012C'],
    '312': ['DoubleElectron2012C'],
    '313': ['DoubleMuon2012C'],
    '314': ['MuEG2012C'],
    '320': ['SingleElectron2012C'],
    '321': ['SingleMuon2012C'],
    '322': ['DoubleElectron2012C'],
    '323': ['DoubleMuon2012C'],
    '324': ['MuEG2012C'],

}
kk={
	'019':['TTTo2L2Nu2B',0.0021922623],
}
data30X={

    '305': ['SingleElectron2012C'],
    '306': ['SingleMuon2012C'],
    '307': ['DoubleElectron2012C'],
    '308': ['DoubleMuon2012C'],
    '309': ['MuEG2012C'],

}
data33X={

    '330': ['SingleElectron2012C'],
    '331': ['SingleMuon2012C'],
    '332': ['DoubleElectron2012C'],
    '333': ['DoubleMuon2012C'],
    '334': ['MuEG2012C'],

}
ggWW={
	'001': ['GluGluToWWTo4L', 0.0016624874]
}
ZZ={
	'079':['ZZTo2L2QMad',0.0007887541],
	'075':['ZZJetsMad',0.0003902703]
}
WZ={
     '074': ['WZJetsMad',0.0005240456],
     '078': ['WZTo2L2QMad',0.0020455233],   

}
firstSamples={
        '080':['WJetsToLNuMad',2.0392984539],
        '011':['TtWFullDR',0.0224598017],
        '012':['TbartWFullDR',0.0226508734],
        '010':['TTJetsMad',0.0325252934],
     '082': ['WGstarToElNuMad',0.0186650056],
     '083': ['WGstarToMuNuMad',0.0063805743],
     '084': ['WGstarToTauNuMad',0.0067204032],
     '085': ['WgammaToLNuG',0.1350168515],   
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
        '080':['WJetsToLNuMad',0.4928715355],
}
TTJets={

        '010':['TTJetsMad',0.0348202597],
}

DYJets={
	'036':['DY10toLLMad',0.1206496207]
}
ZJets={
        '037':['DY50toLLMad',0.1159835077],

}
WWJets={
        '000':['WWJets2LMad',0.0030065150],
}
mcSummer12={
        '070':['WZFull',0.0032327020],
        '077':['ZZTo2L2Nu',0.0007327844],
        '011':['TtWFullDR',0.0224598017],
        '012':['TbartWFullDR',0.0226508734],
#        '037':['DY50toLLMad',1.1935135135],
#        '036':['DY10toLLMad',0.1223652890],
        '001':['GluGluToWWTo4L',0.0016624874],
#        '010':['TTJetsMad',0.0334311887],
        '000':['WWJets2LMad',0.0030065150],
}
hww125={
        '1125':['ggToH125toWWTo2LAndTau2Nu',0.0014792771]
}

hww145={
        '1145':['ggToH145toWWTo2LAndTau2Nu',0.0030737075],
}

hww155={
        '1155':['ggToH155toWWTo2LAndTau2Nu',0.0035425465],
}
hww={
	'1110':['ggToH110toWWTo2LAndTau2Nu',0.0004262067],
	'1115':['ggToH115toWWTo2LAndTau2Nu',0.0006966061],
	'1120':['ggToH120toWWTo2LAndTau2Nu',0.0011154932],
	'1125':['ggToH125toWWTo2LAndTau2Nu',0.0014792771],
	'1130':['ggToH130toWWTo2LAndTau2Nu',0.0019329875],
	'1135':['ggToH135toWWTo2LAndTau2Nu',0.0024106907],
	'1140':['ggToH140toWWTo2LAndTau2Nu',0.0027582844],
	'1145':['ggToH145toWWTo2LAndTau2Nu',0.0030737075],
	'1150':['ggToH150toWWTo2LAndTau2Nu',0.0033341662],
	'1155':['ggToH155toWWTo2LAndTau2Nu',0.0035425465],
	'1160':['ggToH160toWWTo2LAndTau2Nu',0.0038093945],
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



wztt = {
     '3110': ['wzttH110ToWW',0.0004379237],
     '3115': ['wzttH115ToWW',0.0006895429],
#     '3118': ['wzttH118ToWW',1],
     '3120': ['wzttH120ToWW',0.0010049087],
#     '3122': ['wzttH122ToWW',1],
#     '3124': ['wzttH124ToWW',1],
     '3125': ['wzttH125ToWW',0.0013161051],
#     '3126': ['wzttH126ToWW',1],
#     '3128': ['wzttH128ToWW',1],
     '3130': ['wzttH130ToWW',0.0016506222],
     '3135': ['wzttH135ToWW',0.0019024457],
     '3140': ['wzttH140ToWW',0.0021218560],
     '3145': ['wzttH145ToWW',0.0022285293],
     '3150': ['wzttH150ToWW',0.0022945156],
     '3155': ['wzttH155ToWW',0.0023166058],
     '3160': ['wzttH160ToWW',0.0023124900],
     '3170': ['wzttH170ToWW',0.0020375967],
     '3180': ['wzttH180ToWW',0.0015964283],
     '3190': ['wzttH190ToWW',0.0011153262],
     '3200': ['wzttH200ToWW',0.0008722442],
     '3250': ['wzttH250ToWW',0.0003549525],
     '3300': ['wzttH300ToWW',0.0001708895],
#     '3350': ['wzttH350ToWW',1],
#     '3400': ['wzttH400ToWW',1],
#     '3450': ['wzttH450ToWW',1],
#     '3500': ['wzttH500ToWW',1],
#     '3550': ['wzttH550ToWW',1],
#     '3600': ['wzttH600ToWW',1],
}




vbf = {
     '2110': ['vbfToH110toWWTo2LAndTau2Nu',0.0000305978],
     '2115': ['vbfToH115toWWTo2LAndTau2Nu',0.0000530951],
#     '2118': ['vbfToH118toWWTo2LAndTau2Nu',1],
     '2120': ['vbfToH120toWWTo2LAndTau2Nu',0.0000838699],
#     '2122': ['vbfToH122toWWTo2LAndTau2Nu',1],
#     '2124': ['vbfToH124toWWTo2LAndTau2Nu',1],
     '2125': ['vbfToH125toWWTo2LAndTau2Nu',0.0001227694],
#     '2126': ['vbfToH126toWWTo2LAndTau2Nu',1],
#     '2128': ['vbfToH128toWWTo2LAndTau2Nu',1],
     '2130': ['vbfToH130toWWTo2LAndTau2Nu',0.0002485194],
     '2135': ['vbfToH135toWWTo2LAndTau2Nu',0.0002111442],
     '2140': ['vbfToH140toWWTo2LAndTau2Nu',0.0002532907],
     '2145': ['vbfToH145toWWTo2LAndTau2Nu',0.0004337029],
     '2150': ['vbfToH150toWWTo2LAndTau2Nu',0.0004799856],
     '2155': ['vbfToH155toWWTo2LAndTau2Nu',0.0003474926],
     '2160': ['vbfToH160toWWTo2LAndTau2Nu',0.0003743910],
     '2170': ['vbfToH170toWWTo2LAndTau2Nu',0.0003678044],
     '2180': ['vbfToH180toWWTo2LAndTau2Nu',0.0003310326],
     '2190': ['vbfToH190toWWTo2LAndTau2Nu',0.0002799921],
     '2200': ['vbfToH200toWWTo2LAndTau2Nu',0.0002504096],
     '2250': ['vbfToH250toWWTo2LAndTau2Nu',0.0001728683],
     '2300': ['vbfToH300toWWTo2LAndTau2Nu',0.0001260598],
     '2350': ['vbfToH350toWWTo2LAndTau2Nu',0.0000897655],
     '2400': ['vbfToH400toWWTo2LAndTau2Nu',0.0000616117],
     '2450': ['vbfToH450toWWTo2LAndTau2Nu',0.0000469738],
     '2500': ['vbfToH500toWWTo2LAndTau2Nu',0.0000357958],
     '2550': ['vbfToH550toWWTo2LAndTau2Nu',0.0000277184],
     '2600': ['vbfToH600toWWTo2LAndTau2Nu',0.0000228322],
     '2700': ['vbfToH700toWWTo2LAndTau2Nu',0.0000155016],
     '2800': ['vbfToH800toWWTo2LAndTau2Nu',0.0000110652],
     '2900': ['vbfToH900toWWTo2LAndTau2Nu',0.0000082625],
     '3000': ['vbfToH1000toWWTo2LAndTau2Nu',0.0000064149],
}



VV = {
     '074': ['WZJetsMad',0.0008025050],
     '075': ['ZZJetsMad',0.364868],          # before was ,0.0003875613],
     '078': ['WZTo2L2QMad',0.0009493967],   # new!
     '079': ['ZZTo2L2QMad',0.0019032836],
     '082': ['WGstarToElNuMad',0.0186650056],
     '083': ['WGstarToMuNuMad',0.0063805743],
     '084': ['WGstarToTauNuMad',0.0067204032],
     '085': ['WgammaToLNuG',0.1135745261],   # 0.3180158457],
}


