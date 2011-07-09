import FWCore.ParameterSet.Config as cms
from glob import glob

dataSamples = {
    '073':  ['SingleElectron2011A', ],
    '074':  ['SingleMuon2011A', ],
    '075':  ['DoubleElectron2011A', ],
    '076':  ['DoubleMuon2011A', ],
    '077':  ['MuEG2011A', ],
    '078':  ['SingleElectron2011Av2', ],
    '079':  ['SingleMuon2011Av2', ],
    '080':  ['DoubleElectron2011Av2', ],
    '081':  ['DoubleMuon2011Av2', ],
    '082':  ['MuEG2011Av2', ],
}

data42xSamples = {
    '083':  ['SingleElectron2011Av4', ],
    '084':  ['SingleMuon2011Av4', ],
    '085':  ['DoubleElectron2011Av4', ],
    '086':  ['DoubleMuon2011Av4', ],
    '087':  ['MuEG2011Av4', ],
    '089':  ['SingleElectron2011AMay10', ],
    '090':  ['SingleMuon2011AMay10', ],
    '091':  ['DoubleMuon2011AMay10', ],
    '092':  ['DoubleElectron2011AMay10', ],
    '093':  ['MuEG2011AMay10', ],
}

allBackgroundSamples = {
    '001':  ['WZtoAny', 0.0086320726],
    '002':  ['ZZtoAny', 0.0035141667],
    '003':  ['DYtoElEl', 0.8362295184],
    '004':  ['DYtoMuMu', 0.8396525673],
    '005':  ['DYtoTauTau', 0.8349332880],
    '017B':  ['DY10toElElZ2', 1.7884117951],
    '018B':  ['DY10toMuMuZ2', 1.5517272336],
    '023':  ['TTJetsMad', 0.1352851037],
    '010':  ['sTtoBLNu', 0.0030111098],
    '011':  ['tTtoBLNu', 0.0432392679],
    '012':  ['tWTtoBLNu', 0.0216584222],
    '014':  ['VVJetsTo4L', 0.0066194202],
    '022':  ['GluGluToWWTo4L', 0.0013909850],
    '021':  ['PhotonVJets', 0.1497015495],
    '026':  ['WJetsToLNuMad', 2.0722688028],
}

diBosonSamples = {
    '001':  ['WZtoAny', 0.0086320726],
    '002':  ['ZZtoAny', 0.0035141667],
}

diBosonSamplesSplit = {
    '001P':  ['WZtoAnyP', 0.0086320726 ],
    '002P':  ['ZZtoAnyP', 0.0035141667 ],
    '001NP': ['WZtoAnyNP',0.0086320726 ],
    '002NP': ['ZZtoAnyNP',0.0035141667 ],
}

wzSamples = {
    '001':  ['WZtoAny', 0.0086320726],
}

zzSamples = {
    '002':  ['ZZtoAny', 0.0035141667],
}

dySamples = {
    '003':  ['DYtoElEl', 0.8362295184],
    '004':  ['DYtoMuMu', 0.8396525673],
    '005':  ['DYtoTauTau', 0.8349332880],
    '017B':  ['DY10toElElZ2', 1.7884117951],
    '018B':  ['DY10toMuMuZ2', 1.5517272336],
}

topSamples = {
    '023':  ['TTJetsMad', 0.1352851037],
    '010':  ['sTtoBLNu', 0.0030111098],
    '011':  ['tTtoBLNu', 0.0432392679],
    '012':  ['tWTtoBLNu', 0.0216584222],
}

ttbarSamples = {
    '023':  ['TTJetsMad', 0.1352851037],
}

tWSamples = {
    '010':  ['sTtoBLNu', 0.0030111098],
    '011':  ['tTtoBLNu', 0.0432392679],
    '012':  ['tWTtoBLNu', 0.0216584222],
}

wwSamples = {
    '014':  ['VVJetsTo4L', 0.0066194202],
    '022':  ['GluGluToWWTo4L', 0.0013909850],
}

wwSample = {
    '014':  ['VVJetsTo4L', 0.0066194202],
}

ggWWSample = {
    '022':  ['GluGluToWWTo4L', 0.0013909850],
}

wJetSamples = {
    '026':  ['WJetsToLNuMad', 2.0722688028],
}

wGamSamples = {
    '021':  ['PhotonVJets', 0.1497015495],
}

dyMadSamples = {
    '027':  ['DY10toLLMadD6T', 4.8437201469],
    '028':  ['DY50toLLMadD6T', 1.2222772051],
}

otherBDTSamples = {
    '007':  ['DYtoEEPythia', 0.8087429678],
    '008':  ['DYtoMuMuPythia', 0.7421818096],
    '009':  ['DYtoTauTauPythia', 0.8097417867],
    '015':  ['WWto2L', 0.0410360727],
    '016':  ['WWto2LMcnlo', 0.050165233],
    '016A':  ['WWto2LMcnloUp', 0.0501786167],
    '016B':  ['WWto2LMcnloDown', 0.0501669056],
    '029':  ['WtoENu', 2.0559571743],
    '030':  ['WtoMuNu', 2.0715907467],
    '031':  ['WtoTauNu', 2.1865409793],
    '032':  ['WPlusENu', 3.0791515116],
    '033':  ['WMinusENu', 2.2016361763],
    '034':  ['WPlusTauNu', 3.0823635395],
    '035':  ['WMinusTauNu', 2.156279447],
    '036':  ['WPlusMuNu', 3.1024702764],
    '037':  ['WMinusMuNu', 2.1661890538],
}

h120Samples = {
    '101120':  ['ggToH120toWWto2L2Nu', 0.0010517638],
    '102120':  ['ggToH120toWWtoLNuTauNu', 0.0010091735],
    '103120':  ['ggToH120toWWto2Tau2Nu', 0.0004197224],
    '104120':  ['vbfToH120toWWto2L2Nu', 0.0000764742],
    '105120':  ['vbfToH120toWWtoLNuTauNu', 0.0000764777],
    '106120':  ['vbfToH120toWWto2Tau2Nu', 0.0000318728],
}

h130Samples = {
    '101130':  ['ggToH130toWWto2L2Nu', 0.0018282649],
    '102130':  ['ggToH130toWWtoLNuTauNu', 0.0018281984],
    '103130':  ['ggToH130toWWto2Tau2Nu', 0.0007621790],
    '104130':  ['vbfToH130toWWto2L2Nu', 0.0001483666],
    '105130':  ['vbfToH130toWWtoLNuTauNu', 0.0001483761],
    '106130':  ['vbfToH130toWWto2Tau2Nu', 0.0000621994],
}

h140Samples = {
    '101140':  ['ggToH140toWWto2L2Nu', 0.0026542584],
    '102140':  ['ggToH140toWWtoLNuTauNu', 0.0025956905],
    '103140':  ['ggToH140toWWto2Tau2Nu', 0.0010819165],
    '104140':  ['vbfToH140toWWto2L2Nu', 0.0002287988],
    '105140':  ['vbfToH140toWWtoLNuTauNu', 0.0002288394],
    '106140':  ['vbfToH140toWWto2Tau2Nu', 0.0000925687],
}

h150Samples = {
    '101150':  ['ggToH150toWWto2L2Nu', 0.0031658885],
    '102150':  ['ggToH150toWWtoLNuTauNu', 0.0031194028],
    '103150':  ['ggToH150toWWto2Tau2Nu', 0.0013000788],
    '104150':  ['vbfToH150toWWto2L2Nu', 0.0002821953],
    '105150':  ['vbfToH150toWWtoLNuTauNu', 0.0002822364],
    '106150':  ['vbfToH150toWWto2Tau2Nu', 0.0001183163],
}

h160Samples = {
    '101160':  ['ggToH160toWWto2L2Nu', 0.0036210625],
    '102160':  ['ggToH160toWWtoLNuTauNu', 0.0035009727],
    '103160':  ['ggToH160toWWto2Tau2Nu', 0.0014591572],
    '104160':  ['vbfToH160toWWto2L2Nu', 0.0003359004],
    '105160':  ['vbfToH160toWWtoLNuTauNu', 0.0003358881],
    '106160':  ['vbfToH160toWWto2Tau2Nu', 0.0001395822],
}

h170Samples = {
    '101170':  ['ggToH170toWWto2L2Nu', 0.0032730986],
    '102170':  ['ggToH170toWWtoLNuTauNu', 0.0031766813],
    '103170':  ['ggToH170toWWto2Tau2Nu', 0.0013228275],
    '104170':  ['vbfToH170toWWto2L2Nu', 0.0003344758],
    '105170':  ['vbfToH170toWWtoLNuTauNu', 0.0003489966],
    '106170':  ['vbfToH170toWWto2Tau2Nu', 0.0001380943],
}

h180Samples = {
    '101180':  ['ggToH180toWWto2L2Nu', 0.0026918110],
    '102180':  ['ggToH180toWWtoLNuTauNu', 0.0026722069],
    '103180':  ['ggToH180toWWto2Tau2Nu', 0.0011137207],
    '104180':  ['vbfToH180toWWto2L2Nu', 0.0002949392],
    '105180':  ['vbfToH180toWWtoLNuTauNu', 0.0003007589],
    '106180':  ['vbfToH180toWWto2Tau2Nu', 0.0001213574],
}

h190Samples = {
    '101190':  ['ggToH190toWWto2L2Nu', 0.0020183708],
    '102190':  ['ggToH190toWWtoLNuTauNu', 0.0019737620],
    '103190':  ['ggToH190toWWto2Tau2Nu', 0.0008227896],
    '104190':  ['vbfToH190toWWto2L2Nu', 0.0002291986],
    '105190':  ['vbfToH190toWWtoLNuTauNu', 0.0002343976],
    '106190':  ['vbfToH190toWWto2Tau2Nu', 0.0000955878],
}

h200Samples = {
    '101200':  ['ggToH200toWWto2L2Nu', 0.0016565444],
    '102200':  ['ggToH200toWWtoLNuTauNu', 0.0027607315],
    '103200':  ['ggToH200toWWto2Tau2Nu', 0.0006894566],
    '104200':  ['vbfToH200toWWto2L2Nu', 0.0001984470],
    '105200':  ['vbfToH200toWWtoLNuTauNu', 0.0001986242],
    '106200':  ['vbfToH200toWWto2Tau2Nu', 0.0000819461],
}

h250Samples = {
    '101250':  ['ggToH250toWWto2L2Nu', 0.0009872997],
    '102250':  ['ggToH250toWWtoLNuTauNu', 0.0009873445],
    '103250':  ['ggToH250toWWto2Tau2Nu', 0.0004106745],
    '104250':  ['vbfToH250toWWto2L2Nu', 0.0001265270],
    '105250':  ['vbfToH250toWWtoLNuTauNu', 0.0001267473],
    '106250':  ['vbfToH250toWWto2Tau2Nu', 0.0000531051],
}

h300Samples = {
    '101300':  ['ggToH300toWWto2L2Nu', 0.0007100387],
    '102300':  ['ggToH300toWWtoLNuTauNu', 0.0008154699],
    '103300':  ['ggToH300toWWto2Tau2Nu', 0.0002954814],
    '104300':  ['vbfToH300toWWto2L2Nu', 0.0000877618],
    '105300':  ['vbfToH300toWWtoLNuTauNu', 0.0000876088],
    '106300':  ['vbfToH300toWWto2Tau2Nu', 0.0000364515],
}

h350Samples = {
    '101350':  ['ggToH350toWWto2L2Nu', 0.0006537017],
    '102350':  ['ggToH350toWWtoLNuTauNu', 0.0006536958],
    '103350':  ['ggToH350toWWto2Tau2Nu', 0.0002727686],
    '104350':  ['vbfToH350toWWto2L2Nu', 0.0000613519],
    '105350':  ['vbfToH350toWWtoLNuTauNu', 0.0000612616],
    '106350':  ['vbfToH350toWWto2Tau2Nu', 0.0000258374],
}

h400Samples = {
    '101400':  ['ggToH400toWWto2L2Nu', 0.0005010275],
    '102400':  ['ggToH400toWWtoLNuTauNu', 0.0005010093],
    '103400':  ['ggToH400toWWto2Tau2Nu', 0.0002091194],
    '104400':  ['vbfToH400toWWto2L2Nu', 0.0000403211],
    '105400':  ['vbfToH400toWWtoLNuTauNu', 0.0000439065],
    '106400':  ['vbfToH400toWWto2Tau2Nu', 0.0000167418],
}

h450Samples = {
    '101450':  ['ggToH450toWWto2L2Nu', 0.0003200698],
    '102450':  ['ggToH450toWWtoLNuTauNu', 0.0003200582],
    '103450':  ['ggToH450toWWto2Tau2Nu', 0.0001333616],
    '104450':  ['vbfToH450toWWto2L2Nu', 0.0000292633],
    '105450':  ['vbfToH450toWWtoLNuTauNu', 0.0000292593],
    '106450':  ['vbfToH450toWWto2Tau2Nu', 0.0000121895],
}

h500Samples = {
    '101500':  ['ggToH500toWWto2L2Nu', 0.0002000491],
    '102500':  ['ggToH500toWWtoLNuTauNu', 0.0002000327],
    '103500':  ['ggToH500toWWto2Tau2Nu', 0.0000833498],
    '104500':  ['vbfToH500toWWto2L2Nu', 0.0000225270],
    '105500':  ['vbfToH500toWWtoLNuTauNu', 0.0000239490],
    '106500':  ['vbfToH500toWWto2Tau2Nu', 0.0000091398],
}

h550Samples = {
    '101550':  ['ggToH550toWWto2L2Nu', 0.0001254819],
    '102550':  ['ggToH550toWWtoLNuTauNu', 0.0001282468],
    '103550':  ['ggToH550toWWto2Tau2Nu', 0.0000530440],
    '104550':  ['vbfToH550toWWto2L2Nu', 0.0000178338],
    '105550':  ['vbfToH550toWWtoLNuTauNu', 0.0000189596],
    '106550':  ['vbfToH550toWWto2Tau2Nu', 0.0000076165],
}

h600Samples = {
    '101600':  ['ggToH600toWWto2L2Nu', 0.0000791103],
    '102600':  ['ggToH600toWWtoLNuTauNu', 0.0000808513],
    '103600':  ['ggToH600toWWto2Tau2Nu', 0.0000333450],
    '104600':  ['vbfToH600toWWto2L2Nu', 0.0000140794],
    '105600':  ['vbfToH600toWWtoLNuTauNu', 0.0000149681],
    '106600':  ['vbfToH600toWWto2Tau2Nu', 0.0000061034],
}

allSignalSamples = {
    '101120':  ['ggToH120toWWto2L2Nu', 0.0010517638],
    '102120':  ['ggToH120toWWtoLNuTauNu', 0.0010091735],
    '103120':  ['ggToH120toWWto2Tau2Nu', 0.0004197224],
    '104120':  ['vbfToH120toWWto2L2Nu', 0.0000764742],
    '105120':  ['vbfToH120toWWtoLNuTauNu', 0.0000764777],
    '106120':  ['vbfToH120toWWto2Tau2Nu', 0.0000318728],
    '101130':  ['ggToH130toWWto2L2Nu', 0.0018282649],
    '102130':  ['ggToH130toWWtoLNuTauNu', 0.0018281984],
    '103130':  ['ggToH130toWWto2Tau2Nu', 0.0007621790],
    '104130':  ['vbfToH130toWWto2L2Nu', 0.0001483666],
    '105130':  ['vbfToH130toWWtoLNuTauNu', 0.0001483761],
    '106130':  ['vbfToH130toWWto2Tau2Nu', 0.0000621994],
    '101140':  ['ggToH140toWWto2L2Nu', 0.0026542584],
    '102140':  ['ggToH140toWWtoLNuTauNu', 0.0025956905],
    '103140':  ['ggToH140toWWto2Tau2Nu', 0.0010819165],
    '104140':  ['vbfToH140toWWto2L2Nu', 0.0002287988],
    '105140':  ['vbfToH140toWWtoLNuTauNu', 0.0002288394],
    '106140':  ['vbfToH140toWWto2Tau2Nu', 0.0000925687],
    '101150':  ['ggToH150toWWto2L2Nu', 0.0031658885],
    '102150':  ['ggToH150toWWtoLNuTauNu', 0.0031194028],
    '103150':  ['ggToH150toWWto2Tau2Nu', 0.0013000788],
    '104150':  ['vbfToH150toWWto2L2Nu', 0.0002821953],
    '105150':  ['vbfToH150toWWtoLNuTauNu', 0.0002822364],
    '106150':  ['vbfToH150toWWto2Tau2Nu', 0.0001183163],
    '101160':  ['ggToH160toWWto2L2Nu', 0.0036210625],
    '102160':  ['ggToH160toWWtoLNuTauNu', 0.0035009727],
    '103160':  ['ggToH160toWWto2Tau2Nu', 0.0014591572],
    '104160':  ['vbfToH160toWWto2L2Nu', 0.0003359004],
    '105160':  ['vbfToH160toWWtoLNuTauNu', 0.0003358881],
    '106160':  ['vbfToH160toWWto2Tau2Nu', 0.0001395822],
    '101170':  ['ggToH170toWWto2L2Nu', 0.0032730986],
    '102170':  ['ggToH170toWWtoLNuTauNu', 0.0031766813],
    '103170':  ['ggToH170toWWto2Tau2Nu', 0.0013228275],
    '104170':  ['vbfToH170toWWto2L2Nu', 0.0003344758],
    '105170':  ['vbfToH170toWWtoLNuTauNu', 0.0003489966],
    '106170':  ['vbfToH170toWWto2Tau2Nu', 0.0001380943],
    '101180':  ['ggToH180toWWto2L2Nu', 0.0026918110],
    '102180':  ['ggToH180toWWtoLNuTauNu', 0.0026722069],
    '103180':  ['ggToH180toWWto2Tau2Nu', 0.0011137207],
    '104180':  ['vbfToH180toWWto2L2Nu', 0.0002949392],
    '105180':  ['vbfToH180toWWtoLNuTauNu', 0.0003007589],
    '106180':  ['vbfToH180toWWto2Tau2Nu', 0.0001213574],
    '101190':  ['ggToH190toWWto2L2Nu', 0.0020183708],
    '102190':  ['ggToH190toWWtoLNuTauNu', 0.0019737620],
    '103190':  ['ggToH190toWWto2Tau2Nu', 0.0008227896],
    '104190':  ['vbfToH190toWWto2L2Nu', 0.0002291986],
    '105190':  ['vbfToH190toWWtoLNuTauNu', 0.0002343976],
    '106190':  ['vbfToH190toWWto2Tau2Nu', 0.0000955878],
    '101200':  ['ggToH200toWWto2L2Nu', 0.0016565444],
    '102200':  ['ggToH200toWWtoLNuTauNu', 0.0027607315],
    '103200':  ['ggToH200toWWto2Tau2Nu', 0.0006894566],
    '104200':  ['vbfToH200toWWto2L2Nu', 0.0001984470],
    '105200':  ['vbfToH200toWWtoLNuTauNu', 0.0001986242],
    '106200':  ['vbfToH200toWWto2Tau2Nu', 0.0000819461],
    '101250':  ['ggToH250toWWto2L2Nu', 0.0009872997],
    '102250':  ['ggToH250toWWtoLNuTauNu', 0.0009873445],
    '103250':  ['ggToH250toWWto2Tau2Nu', 0.0004106745],
    '104250':  ['vbfToH250toWWto2L2Nu', 0.0001265270],
    '105250':  ['vbfToH250toWWtoLNuTauNu', 0.0001267473],
    '106250':  ['vbfToH250toWWto2Tau2Nu', 0.0000531051],
    '101300':  ['ggToH300toWWto2L2Nu', 0.0007100387],
    '102300':  ['ggToH300toWWtoLNuTauNu', 0.0008154699],
    '103300':  ['ggToH300toWWto2Tau2Nu', 0.0002954814],
    '104300':  ['vbfToH300toWWto2L2Nu', 0.0000877618],
    '105300':  ['vbfToH300toWWtoLNuTauNu', 0.0000876088],
    '106300':  ['vbfToH300toWWto2Tau2Nu', 0.0000364515],
    '101350':  ['ggToH350toWWto2L2Nu', 0.0006537017],
    '102350':  ['ggToH350toWWtoLNuTauNu', 0.0006536958],
    '103350':  ['ggToH350toWWto2Tau2Nu', 0.0002727686],
    '104350':  ['vbfToH350toWWto2L2Nu', 0.0000613519],
    '105350':  ['vbfToH350toWWtoLNuTauNu', 0.0000612616],
    '106350':  ['vbfToH350toWWto2Tau2Nu', 0.0000258374],
    '101400':  ['ggToH400toWWto2L2Nu', 0.0005010275],
    '102400':  ['ggToH400toWWtoLNuTauNu', 0.0005010093],
    '103400':  ['ggToH400toWWto2Tau2Nu', 0.0002091194],
    '104400':  ['vbfToH400toWWto2L2Nu', 0.0000403211],
    '105400':  ['vbfToH400toWWtoLNuTauNu', 0.0000439065],
    '106400':  ['vbfToH400toWWto2Tau2Nu', 0.0000167418],
    '101450':  ['ggToH450toWWto2L2Nu', 0.0003200698],
    '102450':  ['ggToH450toWWtoLNuTauNu', 0.0003200582],
    '103450':  ['ggToH450toWWto2Tau2Nu', 0.0001333616],
    '104450':  ['vbfToH450toWWto2L2Nu', 0.0000292633],
    '105450':  ['vbfToH450toWWtoLNuTauNu', 0.0000292593],
    '106450':  ['vbfToH450toWWto2Tau2Nu', 0.0000121895],
    '101500':  ['ggToH500toWWto2L2Nu', 0.0002000491],
    '102500':  ['ggToH500toWWtoLNuTauNu', 0.0002000327],
    '103500':  ['ggToH500toWWto2Tau2Nu', 0.0000833498],
    '104500':  ['vbfToH500toWWto2L2Nu', 0.0000225270],
    '105500':  ['vbfToH500toWWtoLNuTauNu', 0.0000239490],
    '106500':  ['vbfToH500toWWto2Tau2Nu', 0.0000091398],
    '101550':  ['ggToH550toWWto2L2Nu', 0.0001254819],
    '102550':  ['ggToH550toWWtoLNuTauNu', 0.0001282468],
    '103550':  ['ggToH550toWWto2Tau2Nu', 0.0000530440],
    '104550':  ['vbfToH550toWWto2L2Nu', 0.0000178338],
    '105550':  ['vbfToH550toWWtoLNuTauNu', 0.0000189596],
    '106550':  ['vbfToH550toWWto2Tau2Nu', 0.0000076165],
    '101600':  ['ggToH600toWWto2L2Nu', 0.0000791103],
    '102600':  ['ggToH600toWWtoLNuTauNu', 0.0000808513],
    '103600':  ['ggToH600toWWto2Tau2Nu', 0.0000333450],
    '104600':  ['vbfToH600toWWto2L2Nu', 0.0000140794],
    '105600':  ['vbfToH600toWWtoLNuTauNu', 0.0000149681],
    '106600':  ['vbfToH600toWWto2Tau2Nu', 0.0000061034],
}

