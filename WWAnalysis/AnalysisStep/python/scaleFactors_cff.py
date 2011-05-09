import FWCore.ParameterSet.Config as cms
from glob import glob

dataSamples = {
#     '073':          ['SingleElectron2011A'      ], 
    '074':          ['SingleMuon2011A'          ], 
    '075':          ['DoubleElectron2011A'      ], 
    '076':          ['DoubleMuon2011A'          ], 
    '077':          ['MuEG2011A'                ], 
#     '078':          ['SingleElectron2011Av2'    ], 
    '079':          ['SingleMuon2011Av2'        ], 
    '080':          ['DoubleElectron2011Av2'    ], 
    '081':          ['DoubleMuon2011Av2'        ], 
    '082':          ['MuEG2011Av2'              ], 
}

diBosonSamples = {
    '001':             ['WZtoAny',                         0.0086320726 ],
    '002':             ['ZZtoAny',                         0.0036374708 ],
}

dySamples = {
    '003':             ['DYtoElEl',                        0.8362295184 ],
    '004':             ['DYtoMuMu',                        0.8396525673 ],
    '005':             ['DYtoTauTau',                      0.8349332880 ],
    '017B':            ['DY10toElElZ2',                    1.7884117951 ],
    '018B':            ['DY10toMuMuZ2',                    1.5517272336 ],
    '019':             ['DY10toTauTau',                    1.7359419833 ],
}




ttbarSamples = {
    '023':             ['TTJetsMad',                       0.1352851037 ],
    '010':             ['sTtoBLNu',                        0.0030111098 ],
    '011':             ['tTtoBLNu',                        0.0432392679 ],
    '012':             ['tWTtoBLNu',                       0.0216584222 ],
}

wwSamples = {
    '014':             ['VVJetsTo4L',                      0.0066194202 ],
    '022':             ['GluGluToWWTo4L',                  0.0013909850 ],
}

wJetSamples = {
    '026':             ['WJetsToLNuMad',                   2.0722688028 ],
}

otherSamples = {
    '021':             ['PhotonVJets',                     0.1497015495 ],
}


h120Samples = {
    '101120':          ['ggToH120toWWto2L2Nu',             0.0010517638 ],
    '102120':          ['ggToH120toWWtoLNuTauNu',          0.0010092010 ],
    '103120':          ['ggToH120toWWto2Tau2Nu',           0.0004197224 ],
    '104120':          ['vbfToH120toWWto2L2Nu',            0.0000764742 ],
    '105120':          ['vbfToH120toWWtoLNuTauNu',         0.0000764777 ],
    '106120':          ['vbfToH120toWWto2Tau2Nu',          0.0000318728 ],
}

h130Samples = {
    '101130':          ['ggToH130toWWto2L2Nu',             0.0018282649 ],
    '102130':          ['ggToH130toWWtoLNuTauNu',          0.0018283314 ],
    '103130':          ['ggToH130toWWto2Tau2Nu',           0.0007621790 ],
    '104130':          ['vbfToH130toWWto2L2Nu',            0.0001483666 ],
    '105130':          ['vbfToH130toWWtoLNuTauNu',         0.0001573359 ],
    '106130':          ['vbfToH130toWWto2Tau2Nu',          0.0000621994 ],
}

h140Samples = {
    '101140':          ['ggToH140toWWto2L2Nu',             0.0026542584 ],
    '102140':          ['ggToH140toWWtoLNuTauNu',          0.0025954781 ],
    '103140':          ['ggToH140toWWto2Tau2Nu',           0.0010819165 ],
    '104140':          ['vbfToH140toWWto2L2Nu',            0.0002287988 ],
    '105140':          ['vbfToH140toWWtoLNuTauNu',         0.0002230578 ],
    '106140':          ['vbfToH140toWWto2Tau2Nu',          0.0000925687 ],
}

h150Samples = {
    '101150':          ['ggToH150toWWto2L2Nu',             0.0031658885 ],
    '102150':          ['ggToH150toWWtoLNuTauNu',          0.0031193745 ],
    '103150':          ['ggToH150toWWto2Tau2Nu',           0.0013000788 ],
    '104150':          ['vbfToH150toWWto2L2Nu',            0.0002821953 ],
    '105150':          ['vbfToH150toWWtoLNuTauNu',         0.0002821876 ],
    '106150':          ['vbfToH150toWWto2Tau2Nu',          0.0001183163 ],
}

h160Samples = {
    '101160':          ['ggToH160toWWto2L2Nu',             0.0036210625 ],
    '102160':          ['ggToH160toWWtoLNuTauNu',          0.0035012592 ],
    '103160':          ['ggToH160toWWto2Tau2Nu',           0.0014591572 ],
    '104160':          ['vbfToH160toWWto2L2Nu',            0.0003359004 ],
    '105160':          ['vbfToH160toWWtoLNuTauNu',         0.0003446601 ],
    '106160':          ['vbfToH160toWWto2Tau2Nu',          0.0001395822 ],
}

h170Samples = {
    '101170':          ['ggToH170toWWto2L2Nu',             0.0032730986 ],
    '102170':          ['ggToH170toWWtoLNuTauNu',          0.0031764214 ],
    '103170':          ['ggToH170toWWto2Tau2Nu',           0.0013228275 ],
    '104170':          ['vbfToH170toWWto2L2Nu',            0.0003344758 ],
    '105170':          ['vbfToH170toWWtoLNuTauNu',         0.0003314002 ],
    '106170':          ['vbfToH170toWWto2Tau2Nu',          0.0001380943 ],
}

h180Samples = {
    '101180':          ['ggToH180toWWto2L2Nu',             0.0026918110 ],
    '102180':          ['ggToH180toWWtoLNuTauNu',          0.0026720854 ],
    '103180':          ['ggToH180toWWto2Tau2Nu',           0.0011137207 ],
    '104180':          ['vbfToH180toWWto2L2Nu',            0.0002949392 ],
    '105180':          ['vbfToH180toWWtoLNuTauNu',         0.0002921954 ],
    '106180':          ['vbfToH180toWWto2Tau2Nu',          0.0001213574 ],
}

h190Samples = {
    '101190':          ['ggToH190toWWto2L2Nu',             0.0020183708 ],
    '102190':          ['ggToH190toWWtoLNuTauNu',          0.0019739235 ],
    '103190':          ['ggToH190toWWto2Tau2Nu',           0.0008227896 ],
    '104190':          ['vbfToH190toWWto2L2Nu',            0.0002291986 ],
    '105190':          ['vbfToH190toWWtoLNuTauNu',         0.0002406543 ],
    '106190':          ['vbfToH190toWWto2Tau2Nu',          0.0000955878 ],
}

h200Samples = {
    '101200':          ['ggToH200toWWto2L2Nu',             0.0016565444 ],
    '102200':          ['ggToH200toWWtoLNuTauNu',          0.0016564690 ],
    '103200':          ['ggToH200toWWto2Tau2Nu',           0.0006894566 ],
    '104200':          ['vbfToH200toWWto2L2Nu',            0.0001990559 ],
    '105200':          ['vbfToH200toWWtoLNuTauNu',         0.0002042537 ],
    '106200':          ['vbfToH200toWWto2Tau2Nu',          0.0000819461 ],
}

h210Samples = {
    '101210':          ['ggToH210toWWto2L2Nu',             0.0032314781 ],
    '102210':          ['ggToH210toWWtoLNuTauNu',          0.0024228374 ],
    '103210':          ['ggToH210toWWto2Tau2Nu',           0.0022047070 ],
}

h220Samples = {
    '101220':          ['ggToH220toWWto2L2Nu',             0.0021654140 ],
    '102220':          ['ggToH220toWWtoLNuTauNu',          0.0021653812 ],
    '103220':          ['ggToH220toWWto2Tau2Nu',           0.0005409501 ],
}

h230Samples = {
    '101230':          ['ggToH230toWWto2L2Nu',             0.0027710018 ],
    '102230':          ['ggToH230toWWtoLNuTauNu',          0.0019623293 ],
    '103230':          ['ggToH230toWWto2Tau2Nu',           0.0004909388 ],
}

h250Samples = {
    '101250':          ['ggToH250toWWto2L2Nu',             0.0009875061 ],
    '102250':          ['ggToH250toWWtoLNuTauNu',          0.0009873445 ],
    '103250':          ['ggToH250toWWto2Tau2Nu',           0.0004106745 ],
    '104250':          ['vbfToH250toWWto2L2Nu',            0.0001265327 ],
    '105250':          ['vbfToH250toWWtoLNuTauNu',         0.0001298059 ],
    '106250':          ['vbfToH250toWWto2Tau2Nu',          0.0000530930 ],
}

h300Samples = {
    '101300':          ['ggToH300toWWto2L2Nu',             0.0007100387 ],
    '102300':          ['ggToH300toWWtoLNuTauNu',          0.0007100516 ],
    '103300':          ['ggToH300toWWto2Tau2Nu',           0.0002954814 ],
    '104300':          ['vbfToH300toWWto2L2Nu',            0.0000873855 ],
    '105300':          ['vbfToH300toWWtoLNuTauNu',         0.0000873983 ],
    '106300':          ['vbfToH300toWWto2Tau2Nu',          0.0000364149 ],
}

h350Samples = {
    '101350':          ['ggToH350toWWto2L2Nu',             0.0006537017 ],
    '102350':          ['ggToH350toWWtoLNuTauNu',          0.0007507335 ],
    '103350':          ['ggToH350toWWto2Tau2Nu',           0.0002727686 ],
    '104350':          ['vbfToH350toWWto2L2Nu',            0.0000612504 ],
    '105350':          ['vbfToH350toWWtoLNuTauNu',         0.0000610451 ],
    '106350':          ['vbfToH350toWWto2Tau2Nu',          0.0000258198 ],
}

h400Samples = {
    '101400':          ['ggToH400toWWto2L2Nu',             0.0005028106 ],
    '102400':          ['ggToH400toWWtoLNuTauNu',          0.0005009546 ],
    '103400':          ['ggToH400toWWto2Tau2Nu',           0.0002091194 ],
    '104400':          ['vbfToH400toWWto2L2Nu',            0.0000413033 ],
    '105400':          ['vbfToH400toWWtoLNuTauNu',         0.0000401215 ],
    '106400':          ['vbfToH400toWWto2Tau2Nu',          0.0000167183 ],
}

h450Samples = {
    '101450':          ['ggToH450toWWto2L2Nu',             0.0003200611 ],
    '102450':          ['ggToH450toWWtoLNuTauNu',          0.0003200640 ],
    '103450':          ['ggToH450toWWto2Tau2Nu',           0.0001333616 ],
    '104450':          ['vbfToH450toWWto2L2Nu',            0.0000293024 ],
    '105450':          ['vbfToH450toWWtoLNuTauNu',         0.0000292029 ],
    '106450':          ['vbfToH450toWWto2Tau2Nu',          0.0000121758 ],
}

h500Samples = {
    '101500':          ['ggToH500toWWto2L2Nu',             0.0002000473 ],
    '102500':          ['ggToH500toWWtoLNuTauNu',          0.0002000364 ],
    '103500':          ['ggToH500toWWto2Tau2Nu',           0.0000833498 ],
    '104500':          ['vbfToH500toWWto2L2Nu',            0.0000219933 ],
    '105500':          ['vbfToH500toWWtoLNuTauNu',         0.0000219278 ],
    '106500':          ['vbfToH500toWWto2Tau2Nu',          0.0000091422 ],
}

h550Samples = {
    '101550':          ['ggToH550toWWto2L2Nu',             0.0001254819 ],
    '102550':          ['ggToH550toWWtoLNuTauNu',          0.0001254751 ],
    '103550':          ['ggToH550toWWto2Tau2Nu',           0.0000530440 ],
    '104550':          ['vbfToH550toWWto2L2Nu',            0.0000173751 ],
    '105550':          ['vbfToH550toWWtoLNuTauNu',         0.0000173727 ],
    '106550':          ['vbfToH550toWWto2Tau2Nu',          0.0000076165 ],
}

h600Samples = {
    '101600':          ['ggToH600toWWto2L2Nu',             0.0000791103 ],
    '102600':          ['ggToH600toWWtoLNuTauNu',          0.0000808513 ],
    '103600':          ['ggToH600toWWto2Tau2Nu',           0.0000333450 ],
    '104600':          ['vbfToH600toWWto2L2Nu',            0.0000140794 ],
    '105600':          ['vbfToH600toWWtoLNuTauNu',         0.0000149681 ],
    '106600':          ['vbfToH600toWWto2Tau2Nu',          0.0000061034 ],
}
