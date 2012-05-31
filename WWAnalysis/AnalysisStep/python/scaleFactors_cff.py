import FWCore.ParameterSet.Config as cms
from glob import glob



testSample = {
     '001': ['GluGluToWWTo4L',0.0016624874],


}


vbf = {
#     '2110': ['vbfToH110toWWTo2LAndTau2Nu',1],
#     '2115': ['vbfToH115toWWTo2LAndTau2Nu',1],
#     '2118': ['vbfToH118toWWTo2LAndTau2Nu',1],
#     '2120': ['vbfToH120toWWTo2LAndTau2Nu',1],
#     '2122': ['vbfToH122toWWTo2LAndTau2Nu',1],
#     '2124': ['vbfToH124toWWTo2LAndTau2Nu',1],
#     '2125': ['vbfToH125toWWTo2LAndTau2Nu',1],
#     '2126': ['vbfToH126toWWTo2LAndTau2Nu',1],
#     '2128': ['vbfToH128toWWTo2LAndTau2Nu',1],
     '2130': ['vbfToH130toWWTo2LAndTau2Nu',0.0001660782],
#     '2135': ['vbfToH135toWWTo2LAndTau2Nu',1],
#     '2140': ['vbfToH140toWWTo2LAndTau2Nu',1],
#     '2145': ['vbfToH145toWWTo2LAndTau2Nu',1],
#     '2150': ['vbfToH150toWWTo2LAndTau2Nu',1],
#     '2155': ['vbfToH155toWWTo2LAndTau2Nu',1],
#     '2160': ['vbfToH160toWWTo2LAndTau2Nu',1],
#     '2170': ['vbfToH170toWWTo2LAndTau2Nu',1],
#     '2180': ['vbfToH180toWWTo2LAndTau2Nu',1],
#     '2190': ['vbfToH190toWWTo2LAndTau2Nu',1],
#     '2200': ['vbfToH200toWWTo2LAndTau2Nu',1],
#     '2300': ['vbfToH300toWWTo2LAndTau2Nu',1],
#     '2400': ['vbfToH400toWWTo2LAndTau2Nu',1],
#     '2500': ['vbfToH500toWWTo2LAndTau2Nu',1],
#     '2600': ['vbfToH600toWWTo2LAndTau2Nu',1],
#     '2700': ['vbfToH700toWWTo2LAndTau2Nu',1],
#     '2800': ['vbfToH800toWWTo2LAndTau2Nu',1],
#     '2900': ['vbfToH900toWWTo2LAndTau2Nu',1],
#     '3000': ['vbfToH1000toWWTo2LAndTau2Nu',1],
}



wztt = {
#     '3110': ['wzttH110ToWW',1],
     '3115': ['wzttH115ToWW',1],
#     '3118': ['wzttH118ToWW',1],
     '3120': ['wzttH120ToWW',1],
#     '3122': ['wzttH122ToWW',1],
#     '3124': ['wzttH124ToWW',1],
#     '3125': ['wzttH125ToWW',1],
#     '3126': ['wzttH126ToWW',1],
#     '3128': ['wzttH128ToWW',1],
#     '3130': ['wzttH130ToWW',1],
     '3135': ['wzttH135ToWW',1],
     '3140': ['wzttH140ToWW',1],
#     '3145': ['wzttH145ToWW',1],
     '3150': ['wzttH150ToWW',1],
     '3155': ['wzttH155ToWW',1],
#     '3160': ['wzttH160ToWW',1],
     '3170': ['wzttH170ToWW',1],
#     '3180': ['wzttH180ToWW',1],
     '3190': ['wzttH190ToWW',1],
#     '3200': ['wzttH200ToWW',1],
     '3250': ['wzttH250ToWW',1],
#     '3300': ['wzttH300ToWW',1],
#     '3350': ['wzttH350ToWW',1],
#     '3400': ['wzttH400ToWW',1],
#     '3450': ['wzttH450ToWW',1],
#     '3500': ['wzttH500ToWW',1],
#     '3550': ['wzttH550ToWW',1],
#     '3600': ['wzttH600ToWW',1],
}


dyttdata = {
     '103TT': ['DoubleMu2012ADYttem'],
     '113TT': ['DoubleMu2012BDYttem'],
}


dy = {
     '037a': ['DY50toLLMadCH0',0.1160585354],    
     '037b': ['DY50toLLMadCH1',0.1160585354],    
     '037c': ['DY50toLLMadCH2',0.1160585354],    
     '037d': ['DY50toLLMadCH3',0.1160585354],    
     '037e': ['DY50toLLMadCH4',0.1160585354],    
     '037f': ['DY50toLLMadCH5',0.1160585354],    
     '037g': ['DY50toLLMadCH6',0.1160585354],    
     '037h': ['DY50toLLMadCH7',0.1160585354],    
     '037i': ['DY50toLLMadCH8',0.1160585354],    
}


