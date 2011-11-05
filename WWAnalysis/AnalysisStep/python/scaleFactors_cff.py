import FWCore.ParameterSet.Config as cms
# s/^ *\([0-9abc]*\) *\(.*\)/    '\1': ['\2'],/

four1XSamples = {
    '010': ['sTtoBLNu' , 0.0030111098],
    '011': ['tTtoBLNu' , 0.0432392679],
    '012': ['tWTtoBLNu', 0.0216584222],
}

dataSamplesFall = {
    '100': ['SingleElectron2011Av4'],
    '101': ['SingleMuon2011Av4'],
    '102': ['DoubleElectron2011Av4'],
    '103': ['DoubleMuon2011Av4'],
    '104': ['MuEG2011Av4'],

    '120': ['SingleElectron2011Av6'],
    '121': ['SingleMuon2011Av6'],
    '122': ['DoubleElectron2011Av6'],
    '123': ['DoubleMuon2011Av6'],
    '124': ['MuEG2011Av6'],

    '140': ['SingleElectron2011Bv1a'],
    '141': ['SingleMuon2011Bv1a'],
    '142': ['DoubleElectron2011Bv1a'],
    '143': ['DoubleMuon2011Bv1a'],
    '144': ['MuEG2011Bv1a'],

    '150': ['SingleElectron2011AMay10'],
    '151': ['SingleMuon2011AMay10'],
    '152': ['DoubleMuon2011AMay10'],
    '153': ['DoubleElectron2011AMay10'],
    '154': ['MuEG2011AMay10'],

    '160': ['SingleElectron2011AAug05'],
    '161': ['SingleMuon2011AAug05'],
    '162': ['DoubleElectron2011AAug05'],
    '163': ['DoubleMuon2011AAug05'],
    '164': ['MuEG2011AAug05'],
}

dataSamplesB = {
    '140a': ['SingleElectron2011Bv1a'],
    '141a': ['SingleMuon2011Bv1a'],
    '142a': ['DoubleElectron2011Bv1a'],
    '143a': ['DoubleMuon2011Bv1a'],
    '144a': ['MuEG2011Bv1a'],
    '140b': ['SingleElectron2011Bv1b'],
    '141b': ['SingleMuon2011Bv1b'],
    '142b': ['DoubleElectron2011Bv1b'],
    '143b': ['DoubleMuon2011Bv1b'],
    '144b': ['MuEG2011Bv1b'],
    '140c': ['SingleElectron2011Bv1c'],
    '141c': ['SingleMuon2011Bv1c'],
    '142c': ['DoubleElectron2011Bv1c'],
    '143c': ['DoubleMuon2011Bv1c'],
    '144c': ['MuEG2011Bv1c'],
    '140d': ['SingleElectron2011Bv1d'],
    '141d': ['SingleMuon2011Bv1d'],
    '142d': ['DoubleElectron2011Bv1d'],
    '143d': ['DoubleMuon2011Bv1d'],
    '144d': ['MuEG2011Bv1d'],
}


dataSamplesA = {
    '100': ['SingleElectron2011Av4'],
    '101': ['SingleMuon2011Av4'],
    '102': ['DoubleElectron2011Av4'],
    '103': ['DoubleMuon2011Av4'],
    '104': ['MuEG2011Av4'],
    '120': ['SingleElectron2011Av6'],
    '121': ['SingleMuon2011Av6'],
    '122': ['DoubleElectron2011Av6'],
    '123': ['DoubleMuon2011Av6'],
    '124': ['MuEG2011Av6'],
    '130': ['SingleElectron2011Av6v2'],
    '131': ['SingleMuon2011Av6v2'],
    '132': ['DoubleElectron2011Av6v2'],
    '133': ['DoubleMuon2011Av6v2'],
    '134': ['MuEG2011Av6v2'],
    '150': ['SingleElectron2011AMay10'],
    '151': ['SingleMuon2011AMay10'],
    '152': ['DoubleMuon2011AMay10'],
    '153': ['DoubleElectron2011AMay10'],
    '154': ['MuEG2011AMay10'],
    '160': ['SingleElectron2011AAug05'],
    '161': ['SingleMuon2011AAug05'],
    '162': ['DoubleElectron2011AAug05'],
    '163': ['DoubleMuon2011AAug05'],
    '164': ['MuEG2011AAug05'],
}

dataSamples = {
    '100': ['SingleElectron2011Av4'],
    '101': ['SingleMuon2011Av4'],
    '102': ['DoubleElectron2011Av4'],
    '103': ['DoubleMuon2011Av4'],
    '104': ['MuEG2011Av4'],
    '120': ['SingleElectron2011Av6'],
    '121': ['SingleMuon2011Av6'],
    '122': ['DoubleElectron2011Av6'],
    '123': ['DoubleMuon2011Av6'],
    '124': ['MuEG2011Av6'],
    '130': ['SingleElectron2011Av6v2'],
    '131': ['SingleMuon2011Av6v2'],
    '132': ['DoubleElectron2011Av6v2'],
    '133': ['DoubleMuon2011Av6v2'],
    '134': ['MuEG2011Av6v2'],
    '140a': ['SingleElectron2011Bv1a'],
    '141a': ['SingleMuon2011Bv1a'],
    '142a': ['DoubleElectron2011Bv1a'],
    '143a': ['DoubleMuon2011Bv1a'],
    '144a': ['MuEG2011Bv1a'],
    '140b': ['SingleElectron2011Bv1b'],
    '141b': ['SingleMuon2011Bv1b'],
    '142b': ['DoubleElectron2011Bv1b'],
    '143b': ['DoubleMuon2011Bv1b'],
    '144b': ['MuEG2011Bv1b'],
    '140c': ['SingleElectron2011Bv1c'],
    '141c': ['SingleMuon2011Bv1c'],
    '142c': ['DoubleElectron2011Bv1c'],
    '143c': ['DoubleMuon2011Bv1c'],
    '144c': ['MuEG2011Bv1c'],
    '140d': ['SingleElectron2011Bv1d'],
    '141d': ['SingleMuon2011Bv1d'],
    '142d': ['DoubleElectron2011Bv1d'],
    '143d': ['DoubleMuon2011Bv1d'],
    '144d': ['MuEG2011Bv1d'],
    '150': ['SingleElectron2011AMay10'],
    '151': ['SingleMuon2011AMay10'],
    '152': ['DoubleMuon2011AMay10'],
    '153': ['DoubleElectron2011AMay10'],
    '154': ['MuEG2011AMay10'],
    '160': ['SingleElectron2011AAug05'],
    '161': ['SingleMuon2011AAug05'],
    '162': ['DoubleElectron2011AAug05'],
    '163': ['DoubleMuon2011AAug05'],
    '164': ['MuEG2011AAug05'],
}



allBackgroundSamples = {
    '000': ['WWJets2LMad',0.0039976852],
    '001': ['GluGluToWWTo4L',0.0013910103],
    '011': ['TtWFullDR',0.0096636747],
    '012': ['TbartWFullDR',0.0097162413],
    '013': ['TtFull',0.0107482470],
    '014': ['TbartFull',0.0116462861],
    '015': ['TsFull',0.0122705994],
    '016': ['TbarsFull',0.0104362951],
    '019': ['TTTo2L2Nu2B',0.0016549443],
    '030': ['DYtoElEl',0.0567686049],
    '031': ['DYtoMuMu',0.0565481502],
    '032': ['DYtoTauTau',0.0835612165],
    '033': ['DY10toElEl',0.3565305927],
    '034': ['DY10toMuMu',0.3394064576],
    '035': ['DY10toTauTau',1.5089136364],
    '074': ['WZJetsMad',0.0004875988],
    '075': ['ZZJetsMad',0.0002712239],
    '080': ['WJetsToLNuMad',0.3849171055],
    '082': ['WgammaToElNuMad',0.0283125170],
    '083': ['WgammaToMuNuMad',0.0284605980],
    '084': ['WgammaToTauNuMad',0.0286567790],
}

allOtherBackgrounds = {
    '002': ['WWto2L2NuMCatNLO',0.0088698385],
    '003': ['WWto2L2NuMCatNLOUp',0.0088696742],
    '004': ['WWto2L2NuMCatNLODown',0.0088696413],
    '010': ['TTJetsMad',0.0447562801],
    '017': ['TtWFullDS',0.0098946540],
    '018': ['TbartWFullDS',0.0099920140],
    '037': ['DY50toLLMadD6T',0.0840179524],
    '038': ['ZgammaToElElMad',0.0042919330],
    '039': ['ZgammaToMuMuMad',0.0043232610],
    '040': ['ZgammaToTauTauMad',0.0042588620],
    '041': ['ZgammaToNuNuMad',0.0182223700],
    '070': ['WZFull',0.0042670488],
    '071': ['ZZFull',0.0017693896],
    '072': ['GluGluZZ2L2L',0.0000173000],
    '073': ['GluGluZZ4L',0.0000123014],
    '081': ['VGamma',0.1545118876],
}

allSignalSamples = {
    '8110': ['vbfToH110toWWTo2LAndTau2Nu', 0.0001353976],
    '8115': ['vbfToH115toWWTo2LAndTau2Nu', 0.0001226625],
    '9110': ['ggToH110toWWTo2LAndTau2Nu',  0.0003604647],
    '9115': ['ggToH115toWWTo2LAndTau2Nu',  0.0005802742],
    '1120': ['ggToH120toWWto2L2Nu',0.001009155],
    '2120': ['ggToH120toWWtoLNuTauNu',0.0010091643],
    '3120': ['ggToH120toWWto2Tau2Nu',0.0005036730],
    '4120': ['vbfToH120toWWto2L2Nu',0.0000764791],
    '5120': ['vbfToH120toWWtoLNuTauNu',0.0000764763],
    '6120': ['vbfToH120toWWto2Tau2Nu',0.0000397705],
    '1130': ['ggToH130toWWto2L2Nu',0.001828331],
    '2130': ['ggToH130toWWtoLNuTauNu',0.0018442773],
    '3130': ['ggToH130toWWto2Tau2Nu',0.0010668080],
    '4130': ['vbfToH130toWWto2L2Nu',0.0001483693],
    '5130': ['vbfToH130toWWtoLNuTauNu',0.0001483923],
    '6130': ['vbfToH130toWWto2Tau2Nu',0.0000776251],
    '1140': ['ggToH140toWWto2L2Nu',0.002613488],
    '2140': ['ggToH140toWWtoLNuTauNu',0.0025956197],
    '3140': ['ggToH140toWWto2Tau2Nu',0.0014836060],
    '4140': ['vbfToH140toWWto2L2Nu',0.0002230781],
    '5140': ['vbfToH140toWWtoLNuTauNu',0.0002230598],
    '6140': ['vbfToH140toWWto2Tau2Nu',0.0001110787],
    '1150': ['ggToH150toWWto2L2Nu',0.003119261],
    '2150': ['ggToH150toWWtoLNuTauNu',0.0031193178],
    '3150': ['ggToH150toWWto2Tau2Nu',0.0015602270],
    '4150': ['vbfToH150toWWto2L2Nu',0.0002917482],
    '5150': ['vbfToH150toWWtoLNuTauNu',0.0002822004],
    '6150': ['vbfToH150toWWto2Tau2Nu',0.0001419808],
    '1160': ['ggToH160toWWto2L2Nu',0.003619803],
    '2160': ['ggToH160toWWtoLNuTauNu',0.0035011319],
    '3160': ['ggToH160toWWto2Tau2Nu',0.0017511000],
    '4160': ['vbfToH160toWWto2L2Nu',0.0003359096],
    '5160': ['vbfToH160toWWtoLNuTauNu',0.0003359279],
    '6160': ['vbfToH160toWWto2Tau2Nu',0.0001675194],
    '1170': ['ggToH170toWWto2L2Nu',0.003284082],
    '2170': ['ggToH170toWWtoLNuTauNu',0.0031764791],
    '3170': ['ggToH170toWWto2Tau2Nu',0.0017227090],
    '4170': ['vbfToH170toWWto2L2Nu',0.0004050701],
    '5170': ['vbfToH170toWWtoLNuTauNu',0.0003313941],
    '6170': ['vbfToH170toWWto2Tau2Nu',0.0001785504],
    '1180': ['ggToH180toWWto2L2Nu',0.002672061],
    '2180': ['ggToH180toWWtoLNuTauNu',0.0026719396],
    '3180': ['ggToH180toWWto2Tau2Nu',0.0013541897],
    '4180': ['vbfToH180toWWto2L2Nu',0.0002939399],
    '5180': ['vbfToH180toWWtoLNuTauNu',0.0002921980],
    '6180': ['vbfToH180toWWto2Tau2Nu',0.0001568966],
    '1190': ['ggToH190toWWto2L2Nu',0.002013373],
    '2190': ['ggToH190toWWtoLNuTauNu',0.0019737799],
    '3190': ['ggToH190toWWto2Tau2Nu',0.0009873086],
    '4190': ['vbfToH190toWWto2L2Nu',0.0002298619],
    '5190': ['vbfToH190toWWtoLNuTauNu',0.0002284872],
    '6190': ['vbfToH190toWWto2Tau2Nu',0.0001146872],
    '1200': ['ggToH200toWWto2L2Nu',0.001656469],
    '2200': ['ggToH200toWWtoLNuTauNu',0.0016564841],
    '3200': ['ggToH200toWWto2Tau2Nu',0.0008688014],
    '4200': ['vbfToH200toWWto2L2Nu',0.0001984308],
    '5200': ['vbfToH200toWWtoLNuTauNu',0.0002034265],
    '6200': ['vbfToH200toWWto2Tau2Nu',0.0000983338],
    '1210': ['ggToH210toWWto2L2Nu',0.002907484],
    '2210': ['ggToH210toWWtoLNuTauNu',0.0029074842],
    '3210': ['ggToH210toWWto2Tau2Nu',0.0007637523],
    '4210': ['vbfToH210toWWto2L2Nu',0.0003604590],
    '5210': ['vbfToH210toWWtoLNuTauNu',0.0003604720],
    '6210': ['vbfToH210toWWto2Tau2Nu',0.0000892126],
    '1220': ['ggToH220toWWto2L2Nu',0.002598513],
    '2220': ['ggToH220toWWtoLNuTauNu',0.0025983708],
    '3220': ['ggToH220toWWto2Tau2Nu',0.0006491617],
    '4220': ['vbfToH220toWWto2L2Nu',0.0003374910],
    '5220': ['vbfToH220toWWtoLNuTauNu',0.0003295280],
    '6220': ['vbfToH220toWWto2Tau2Nu',0.0000897953],
    '1230': ['ggToH230toWWto2L2Nu',0.00235476],
    '2230': ['ggToH230toWWtoLNuTauNu',0.0023548023],
    '3230': ['ggToH230toWWto2Tau2Nu',0.0006186748],
    '4230': ['vbfToH230toWWto2L2Nu',0.0003021150],
    '5230': ['vbfToH230toWWtoLNuTauNu',0.0003022630],
    '6230': ['vbfToH230toWWto2Tau2Nu',0.0000746405],
    '1250': ['ggToH250toWWto2L2Nu',0.001024161],
    '2250': ['ggToH250toWWtoLNuTauNu',0.0010137122],
    '3250': ['ggToH250toWWto2Tau2Nu',0.0005174620],
    '4250': ['vbfToH250toWWto2L2Nu',0.0001265270],
    '5250': ['vbfToH250toWWtoLNuTauNu',0.0001265201],
    '6250': ['vbfToH250toWWto2Tau2Nu',0.0000637035],
    '1300': ['ggToH300toWWto2L2Nu',0.000710071],
    '2300': ['ggToH300toWWtoLNuTauNu',0.0007100516],
    '3300': ['ggToH300toWWto2Tau2Nu',0.0003723506],
    '4300': ['vbfToH300toWWto2L2Nu',0.0000873839],
    '5300': ['vbfToH300toWWtoLNuTauNu',0.0000873967],
    '6300': ['vbfToH300toWWto2Tau2Nu',0.0000472543],
    '1350': ['ggToH350toWWto2L2Nu',0.00065369],
    '2350': ['ggToH350toWWtoLNuTauNu',0.0006537136],
    '3350': ['ggToH350toWWto2Tau2Nu',0.0003273203],
    '4350': ['vbfToH350toWWto2L2Nu',0.0000610556],
    '5350': ['vbfToH350toWWtoLNuTauNu',0.0001520378],
    '6350': ['vbfToH350toWWto2Tau2Nu',0.0000309812],
    '1400': ['ggToH400toWWto2L2Nu',0.000500941],
    '2400': ['ggToH400toWWtoLNuTauNu',0.0005009319],
    '3400': ['ggToH400toWWto2Tau2Nu',0.0002635197],
    '4400': ['vbfToH400toWWto2L2Nu',0.0000401398],
    '5400': ['vbfToH400toWWtoLNuTauNu',0.0000403922],
    '6400': ['vbfToH400toWWto2Tau2Nu',0.0000200624],
    '1450': ['ggToH450toWWto2L2Nu',0.000320067],
    '2450': ['ggToH450toWWtoLNuTauNu',0.0003200553],
    '3450': ['ggToH450toWWto2Tau2Nu',0.0001715166],
    '4450': ['vbfToH450toWWto2L2Nu',0.0000292107],
    '5450': ['vbfToH450toWWtoLNuTauNu',0.0000292032],
    '6450': ['vbfToH450toWWto2Tau2Nu',0.0000146140],
    '1500': ['ggToH500toWWto2L2Nu',0.000200049],
    '2500': ['ggToH500toWWtoLNuTauNu',0.0002000418],
    '3500': ['ggToH500toWWto2Tau2Nu',0.0001071895],
    '4500': ['vbfToH500toWWto2L2Nu',0.0000220554],
    '5500': ['vbfToH500toWWtoLNuTauNu',0.0000227661],
    '6500': ['vbfToH500toWWto2Tau2Nu',0.0000111208],
    '1550': ['ggToH550toWWto2L2Nu',0.000125485],
    '2550': ['ggToH550toWWtoLNuTauNu',0.0001254762],
    '3550': ['ggToH550toWWto2Tau2Nu',0.0000636560],
    '4550': ['vbfToH550toWWto2L2Nu',0.0000174782],
    '5550': ['vbfToH550toWWtoLNuTauNu',0.0000187541],
    '6550': ['vbfToH550toWWto2Tau2Nu',0.0000091388],
    '1600': ['ggToH600toWWto2L2Nu',0.000079111],
    '2600': ['ggToH600toWWtoLNuTauNu',0.0000791074],
    '3600': ['ggToH600toWWto2Tau2Nu',0.0000428766],
    '4600': ['vbfToH600toWWto2L2Nu',0.0000139747],
    '5600': ['vbfToH600toWWtoLNuTauNu',0.0000138171],
    '6600': ['vbfToH600toWWto2Tau2Nu',0.0000073233],
}

# im disabling these:
# wwSamples = {
#     '000': ['WWJets2LMad',0.0039976852],
#     '001': ['GluGluToWWTo4L',0.0013910103],
# }
# 
# wwSysSamples = {
#     '002': ['WWto2L2NuMCatNLO',0.0088698385],
#     '003': ['WWto2L2NuMCatNLOUp',0.0088696742],
#     '004': ['WWto2L2NuMCatNLODown',0.0088696413],
# }
# 
# ggWWSample = {
#     '001': ['GluGluToWWTo4L',0.0013910103],
# }
# 
# ttbarSamples = {
#     '010': ['TTJetsMad',0.0447562801],
# }
# 
# ttbarPowhegSample={
#      '019': ['TTTo2L2Nu2B',0.0016549443],
#  }
# 
# wzMadgraphSample = {
#      '074': ['WZJetsMad',0.0004875988],
# }
# 
# zzMadgraphSample = {
#      '075': ['ZZJetsMad',0.0002712239],
#  }
# 
# topSamples = {
#     '010': ['TTJetsMad',0.0447562801],
#     '011': ['TtWFullDR',0.0096636747],
#     '012': ['TbartWFullDR',0.0097162413],
#     '013': ['TtFull',0.0107482470],
#     '014': ['TbartFull',0.0116462861],
#     '015': ['TsFull',0.0122705994],
#     '016': ['TbarsFull',0.0104362951],
# }
# 
# tWSamples = {
#     '011': ['TtWFullDR',0.0096636747],
#     '012': ['TbartWFullDR',0.0097162413],
#     '013': ['TtFull',0.0107482470],
#     '014': ['TbartFull',0.0116462861],
#     '015': ['TsFull',0.0122705994],
#     '016': ['TbarsFull',0.0104362951],
# }
# 
# tWSysSamples = {
#     '017': ['TtWFullDS',0.0098946540],
#     '018': ['TbartWFullDS',0.0099920140],
# }
# 
# dySamples = {
#     '030': ['DYtoElEl',0.0567686049],
#     '031': ['DYtoMuMu',0.0565481502],
#     '032': ['DYtoTauTau',0.0835612165],
#     '033': ['DY10toElEl',0.3565305927],
#     '034': ['DY10toMuMu',0.3394064576],
#     '035': ['DY10toTauTau',1.5089136364],
# }
# 
# zGamSamples = {
#     '038': ['ZgammaToElElMad',0.0042919330],
#     '039': ['ZgammaToMuMuMad',0.0043232610],
#     '040': ['ZgammaToTauTauMad',0.0042588620],
#     '041': ['ZgammaToNuNuMad',0.0182223700],
# }
# 
# dyMadSamples = {
#     '037': ['DY50toLLMadD6T',0.0840179520],
# }
# 
# wzSamples = {
#     '070': ['WZFull',0.0042670488],
# }
# 
# zzSamples = {
#     '071': ['ZZFull',0.0018314734],
# }
# 
# diBosonSamples = {
#     '070': ['WZFull',0.0042670488],
#     '071': ['ZZFull',0.0018314734],
#     '072': ['GluGluZZ2L2L',0.0000173000],
#     '073': ['GluGluZZ4L',0.0000123014],
# }
# 
# wJetSamples = {
#     '080': ['WJetsToLNuMad',0.3849171055],
# }
# 
# wGamSamples = {
#     '082': ['WgammaToElNuMad',0.0283125170],
#     '083': ['WgammaToMuNuMad',0.0284605980],
#     '084': ['WgammaToTauNuMad',0.0286567790],
# }
# 
# vGamSamples = {
#     '081': ['VGamma',0.1545118876],
# }
# 
# 
# h115Samples = {
#     '8110': ['vbfToH110toWWTo2LAndTau2Nu', 0.0001353976],
#     '8115': ['vbfToH115toWWTo2LAndTau2Nu', 0.0001226625],
# }
# 
# h110Samples = {
#     '9110': ['ggToH110toWWTo2LAndTau2Nu',  0.0003604647],
#     '9115': ['ggToH115toWWTo2LAndTau2Nu',  0.0005802742],
# }
# 
# h120Samples = {
#     '1120': ['ggToH120toWWto2L2Nu',0.001009155],
#     '2120': ['ggToH120toWWtoLNuTauNu',0.0010091643],
#     '3120': ['ggToH120toWWto2Tau2Nu',0.0005036730],
#     '4120': ['vbfToH120toWWto2L2Nu',0.0000764791],
#     '5120': ['vbfToH120toWWtoLNuTauNu',0.0000764763],
#     '6120': ['vbfToH120toWWto2Tau2Nu',0.0000397705],
# }
# 
# h130Samples = {
#     '1130': ['ggToH130toWWto2L2Nu',0.001828331],
#     '2130': ['ggToH130toWWtoLNuTauNu',0.0018442773],
#     '3130': ['ggToH130toWWto2Tau2Nu',0.0010668080],
#     '4130': ['vbfToH130toWWto2L2Nu',0.0001483693],
#     '5130': ['vbfToH130toWWtoLNuTauNu',0.0001483923],
#     '6130': ['vbfToH130toWWto2Tau2Nu',0.0000776251],
# }
# 
# h140Samples = {
#     '1140': ['ggToH140toWWto2L2Nu',0.002613488],
#     '2140': ['ggToH140toWWtoLNuTauNu',0.0025956197],
#     '3140': ['ggToH140toWWto2Tau2Nu',0.0014836060],
#     '4140': ['vbfToH140toWWto2L2Nu',0.0002230781],
#     '5140': ['vbfToH140toWWtoLNuTauNu',0.0002230598],
#     '6140': ['vbfToH140toWWto2Tau2Nu',0.0001110787],
# }
# 
# h150Samples = {
#     '1150': ['ggToH150toWWto2L2Nu',0.003119261],
#     '2150': ['ggToH150toWWtoLNuTauNu',0.0031193178],
#     '3150': ['ggToH150toWWto2Tau2Nu',0.0015602270],
#     '4150': ['vbfToH150toWWto2L2Nu',0.0002917482],
#     '5150': ['vbfToH150toWWtoLNuTauNu',0.0002822004],
#     '6150': ['vbfToH150toWWto2Tau2Nu',0.0001419808],
# }
# 
# h160Samples = {
#     '1160': ['ggToH160toWWto2L2Nu',0.003619803],
#     '2160': ['ggToH160toWWtoLNuTauNu',0.0035011319],
#     '3160': ['ggToH160toWWto2Tau2Nu',0.0017511000],
#     '4160': ['vbfToH160toWWto2L2Nu',0.0003359096],
#     '5160': ['vbfToH160toWWtoLNuTauNu',0.0003359279],
#     '6160': ['vbfToH160toWWto2Tau2Nu',0.0001675194],
# }
# 
# h170Samples = {
#     '1170': ['ggToH170toWWto2L2Nu',0.003284082],
#     '2170': ['ggToH170toWWtoLNuTauNu',0.0031764791],
#     '3170': ['ggToH170toWWto2Tau2Nu',0.0017227090],
#     '4170': ['vbfToH170toWWto2L2Nu',0.0004050701],
#     '5170': ['vbfToH170toWWtoLNuTauNu',0.0003313941],
#     '6170': ['vbfToH170toWWto2Tau2Nu',0.0001785504],
# }
# 
# h180Samples = {
#     '1180': ['ggToH180toWWto2L2Nu',0.002672061],
#     '2180': ['ggToH180toWWtoLNuTauNu',0.0026719396],
#     '3180': ['ggToH180toWWto2Tau2Nu',0.0013541897],
#     '4180': ['vbfToH180toWWto2L2Nu',0.0002939399],
#     '5180': ['vbfToH180toWWtoLNuTauNu',0.0002921980],
#     '6180': ['vbfToH180toWWto2Tau2Nu',0.0001568966],
# }
# 
# h190Samples = {
#     '1190': ['ggToH190toWWto2L2Nu',0.002013373],
#     '2190': ['ggToH190toWWtoLNuTauNu',0.0019737799],
#     '3190': ['ggToH190toWWto2Tau2Nu',0.0009873086],
#     '4190': ['vbfToH190toWWto2L2Nu',0.0002298619],
#     '5190': ['vbfToH190toWWtoLNuTauNu',0.0002284872],
#     '6190': ['vbfToH190toWWto2Tau2Nu',0.0001146872],
# }
# 
# h200Samples = {
#     '1200': ['ggToH200toWWto2L2Nu',0.001656469],
#     '2200': ['ggToH200toWWtoLNuTauNu',0.0016564841],
#     '3200': ['ggToH200toWWto2Tau2Nu',0.0008688014],
#     '4200': ['vbfToH200toWWto2L2Nu',0.0001984308],
#     '5200': ['vbfToH200toWWtoLNuTauNu',0.0002034265],
#     '6200': ['vbfToH200toWWto2Tau2Nu',0.0000983338],
# }
# 
# h210Samples = {
#     '1210': ['ggToH210toWWto2L2Nu',0.002907484],
#     '2210': ['ggToH210toWWtoLNuTauNu',0.0029074842],
#     '3210': ['ggToH210toWWto2Tau2Nu',0.0007637523],
#     '4210': ['vbfToH210toWWto2L2Nu',0.0003604590],
#     '5210': ['vbfToH210toWWtoLNuTauNu',0.0003604720],
#     '6210': ['vbfToH210toWWto2Tau2Nu',0.0000892126],
# }
# 
# h220Samples = {
#     '1220': ['ggToH220toWWto2L2Nu',0.002598513],
#     '2220': ['ggToH220toWWtoLNuTauNu',0.0025983708],
#     '3220': ['ggToH220toWWto2Tau2Nu',0.0006491617],
#     '4220': ['vbfToH220toWWto2L2Nu',0.0003374910],
#     '5220': ['vbfToH220toWWtoLNuTauNu',0.0003295280],
#     '6220': ['vbfToH220toWWto2Tau2Nu',0.0000897953],
# }
# 
# h230Samples = {
#     '1230': ['ggToH230toWWto2L2Nu',0.00235476],
#     '2230': ['ggToH230toWWtoLNuTauNu',0.0023548023],
#     '3230': ['ggToH230toWWto2Tau2Nu',0.0006186748],
#     '4230': ['vbfToH230toWWto2L2Nu',0.0003021150],
#     '5230': ['vbfToH230toWWtoLNuTauNu',0.0003022630],
#     '6230': ['vbfToH230toWWto2Tau2Nu',0.0000746405],
# }
# 
# h250Samples = {
#     '1250': ['ggToH250toWWto2L2Nu',0.001024161],
#     '2250': ['ggToH250toWWtoLNuTauNu',0.0010137122],
#     '3250': ['ggToH250toWWto2Tau2Nu',0.0005174620],
#     '4250': ['vbfToH250toWWto2L2Nu',0.0001265270],
#     '5250': ['vbfToH250toWWtoLNuTauNu',0.0001265201],
#     '6250': ['vbfToH250toWWto2Tau2Nu',0.0000637035],
# }
# 
# h300Samples = {
#     '1300': ['ggToH300toWWto2L2Nu',0.000710071],
#     '2300': ['ggToH300toWWtoLNuTauNu',0.0007100516],
#     '3300': ['ggToH300toWWto2Tau2Nu',0.0003723506],
#     '4300': ['vbfToH300toWWto2L2Nu',0.0000873839],
#     '5300': ['vbfToH300toWWtoLNuTauNu',0.0000873967],
#     '6300': ['vbfToH300toWWto2Tau2Nu',0.0000472543],
# }
# 
# h350Samples = {
#     '1350': ['ggToH350toWWto2L2Nu',0.00065369],
#     '2350': ['ggToH350toWWtoLNuTauNu',0.0006537136],
#     '3350': ['ggToH350toWWto2Tau2Nu',0.0003273203],
#     '4350': ['vbfToH350toWWto2L2Nu',0.0000610556],
#     '5350': ['vbfToH350toWWtoLNuTauNu',0.0001520378],
#     '6350': ['vbfToH350toWWto2Tau2Nu',0.0000309812],
# }
# 
# h400Samples = {
#     '1400': ['ggToH400toWWto2L2Nu',0.000500941],
#     '2400': ['ggToH400toWWtoLNuTauNu',0.0005009319],
#     '3400': ['ggToH400toWWto2Tau2Nu',0.0002635197],
#     '4400': ['vbfToH400toWWto2L2Nu',0.0000401398],
#     '5400': ['vbfToH400toWWtoLNuTauNu',0.0000403922],
#     '6400': ['vbfToH400toWWto2Tau2Nu',0.0000200624],
# }
# 
# h450Samples = {
#     '1450': ['ggToH450toWWto2L2Nu',0.000320067],
#     '2450': ['ggToH450toWWtoLNuTauNu',0.0003200553],
#     '3450': ['ggToH450toWWto2Tau2Nu',0.0001715166],
#     '4450': ['vbfToH450toWWto2L2Nu',0.0000292107],
#     '5450': ['vbfToH450toWWtoLNuTauNu',0.0000292032],
#     '6450': ['vbfToH450toWWto2Tau2Nu',0.0000146140],
# }
# 
# h500Samples = {
#     '1500': ['ggToH500toWWto2L2Nu',0.000200049],
#     '2500': ['ggToH500toWWtoLNuTauNu',0.0002000418],
#     '3500': ['ggToH500toWWto2Tau2Nu',0.0001071895],
#     '4500': ['vbfToH500toWWto2L2Nu',0.0000220554],
#     '5500': ['vbfToH500toWWtoLNuTauNu',0.0000227661],
#     '6500': ['vbfToH500toWWto2Tau2Nu',0.0000111208],
# }
# 
# h550Samples = {
#     '1550': ['ggToH550toWWto2L2Nu',0.000125485],
#     '2550': ['ggToH550toWWtoLNuTauNu',0.0001254762],
#     '3550': ['ggToH550toWWto2Tau2Nu',0.0000636560],
#     '4550': ['vbfToH550toWWto2L2Nu',0.0000174782],
#     '5550': ['vbfToH550toWWtoLNuTauNu',0.0000187541],
#     '6550': ['vbfToH550toWWto2Tau2Nu',0.0000091388],
# }
# 
# h600Samples = {
#     '1600': ['ggToH600toWWto2L2Nu',0.000079111],
#     '2600': ['ggToH600toWWtoLNuTauNu',0.0000791074],
#     '3600': ['ggToH600toWWto2Tau2Nu',0.0000428766],
#     '4600': ['vbfToH600toWWto2L2Nu',0.0000139747],
#     '5600': ['vbfToH600toWWtoLNuTauNu',0.0000138171],
#     '6600': ['vbfToH600toWWto2Tau2Nu',0.0000073233],
# }
# 
# temp7000s = {
#     '7120': ['zttH120ToWW', 0.0050611818],
#     '7130': ['zttH130ToWW', 0.0038991451],
#     '7140': ['zttH140ToWW', 0.0030164545],
#     '7150': ['zttH150ToWW', 0.0023640455],
#     '7160': ['zttH160ToWW', 0.0018269091],
#     '7170': ['zttH170ToWW', 0.0015049545],
#     '7180': ['zttH180ToWW', 0.0012172273],
#     '7190': ['zttH190ToWW', 0.0010046364],
#     '7200': ['zttH200ToWW', 0.0006146818],
#     '7250': ['zttH250ToWW', 0.0003503182],
#     '7300': ['zttH300ToWW', 0.0001667814],
# }
# 
