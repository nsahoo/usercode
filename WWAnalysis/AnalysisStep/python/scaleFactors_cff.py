import FWCore.ParameterSet.Config as cms
from glob import glob

dataSamples = {
    '073': ['SingleElectron2011A','/SingleElectron/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID073_SingleElectron2011A_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '074': ['SingleMuon2011A','/SingleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID074_SingleMuon2011A_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '075': ['DoubleElectron2011A','/DoubleElectron/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID075_DoubleElectron2011A_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '076': ['DoubleMuon2011A','/DoubleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID076_DoubleMuon2011A_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '077': ['MuEG2011A','/MuEG/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID077_MuEG2011A_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '078': ['SingleElectron2011Av2','/SingleElectron/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID078_SingleElectron2011Av2_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '079': ['SingleMuon2011Av2','/SingleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID079_SingleMuon2011Av2_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '080': ['DoubleElectron2011Av2','/DoubleElectron/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID080_DoubleElectron2011Av2_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '081': ['DoubleMuon2011Av2','/DoubleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID081_DoubleMuon2011Av2_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
    '082': ['MuEG2011Av2','/MuEG/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID082_MuEG2011Av2_S2-e29d23d33eb931e4b45d722c65b39f05/USER'],
}

dy42XSamples = {
    '030_42X': ['DYtoElEl',    0.0567686049, '/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID030_DYtoElEl-a0fee9021ee0e974e473d874349bbb62/USER'],
    '031_42X': ['DYtoMuMu',    0.0565481502, '/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID031_DYtoMuMu-a0fee9021ee0e974e473d874349bbb62/USER'],
    '032_42X': ['DYtoTauTau',  0.0835612165, '/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID032_DYtoTauTau-a0fee9021ee0e974e473d874349bbb62/USER'],
    '033_42X': ['DY10toElEl',  1.6292217438, '/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID033_DY10toElEl-a0fee9021ee0e974e473d874349bbb62/USER'],
    '034_42X': ['DY10toMuMu',  1.5767956975, '/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID034_DY10toMuMu-a0fee9021ee0e974e473d874349bbb62/USER'],
    '035_42X': ['DY10toTauTau',1.5713636364, '/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID035_DY10toTauTau-a0fee9021ee0e974e473d874349bbb62/USER'],
}


newDataSamples = {
    '100': ['SingleElectron2011Av4',   '/SingleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID100_SingleElectron2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '101': ['SingleMuon2011Av4',   '/SingleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID101_SingleMuon2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '102': ['DoubleElectron2011Av4',   '/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID102_DoubleElectron2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '103': ['DoubleMuon2011Av4',   '/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID103_DoubleMuon2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '104': ['MuEG2011Av4', '/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID104_MuEG2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER'],
                                
    '121': ['SingleMuon2011Av6',   '/SingleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID121_SingleMuon2011Av6-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '122': ['DoubleElectron2011Av6',   '/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID122_DoubleElectron2011Av6-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '123': ['DoubleMuon2011Av6',   '/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID123_DoubleMuon2011Av6-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '124': ['MuEG2011Av6', '/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID124_MuEG2011Av6v3-eff9861cde6b1617ce28065ed0ad5d45/USER'],
                                
    '150': ['SingleElectron2011AMay10',    '/SingleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID150_SingleElectron2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '151': ['SingleMuon2011AMay10',    '/SingleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID151_SingleMuon2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '152': ['DoubleMuon2011AMay10',    '/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID152_DoubleMuon2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '153': ['DoubleElectron2011AMay10',    '/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID153_DoubleElectron2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '154': ['MuEG2011AMay10',  '/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID154_MuEG2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER'],
                                
    '160': ['SingleElectron2011AAug05',    '/SingleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID160_SingleElectron2011AAug05v2-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '161': ['SingleMuon2011AAug05',  '/SingleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID161_SingleMu2011AAug05v3-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '162': ['DoubleElectron2011AAug05',    '/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID162_DoubleElectron2011AAug05v2-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '163': ['DoubleMuon2011AAug05',  '/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID163_DoubleMu2011AAug05-eff9861cde6b1617ce28065ed0ad5d45/USER'],
    '164': ['MuEG2011AAug05',  '/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID164_MuEG2011AAug05v4-eff9861cde6b1617ce28065ed0ad5d45/USER'],
}

data42xSamples = {
    '083_old': ['SingleElectron2011Av4','/SingleElectron/mangano-R42X_S1_V05_S2_V02_S3_V02_ID083_SingleElectron2011Av4_S2-v3-271b442b03cdae3729dc1e24bf59e029/USER'],
    '084_old': ['SingleMuon2011Av4','/SingleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID084_SingleMuon2011Av4_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],
    '085_old': ['DoubleElectron2011Av4','/DoubleElectron/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID085_DoubleElectron2011Av4_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],
    '086_old': ['DoubleMuon2011Av4','/DoubleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID086_DoubleMuon2011Av4_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],
    '087_old': ['MuEG2011Av4','/MuEG/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID087_MuEG2011Av4_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],

    '083b_old': ['SingleElectron2011Av4','/SingleElectron/mangano-R42X_S1_V05_S2_V02_S3_V02_ID083b_SingleElectron2011Av4_S2-v3-271b442b03cdae3729dc1e24bf59e029/USER'],
    '084b_old': ['SingleMuon2011Av4','/SingleMu/mangano-R42X_S1_V05_S2_V02_S3_V02_ID084b_SingleMuon2011Av4_S2-bb1d7d67c682eb6954dd2da78c9b0e4a/USER'],
    '085b_old': ['DoubleElectron2011Av4','/DoubleElectron/mangano-R42X_S1_V05_S2_V02_S3_V02_ID085b_DoubleElectron2011Av4_S2-bb1d7d67c682eb6954dd2da78c9b0e4a/USER'],
    '086b_old': ['DoubleMuon2011Av4','/DoubleMu/mangano-R42X_S1_V05_S2_V02_S3_V02_ID086b_DoubleMuon2011Av4_S2-bb1d7d67c682eb6954dd2da78c9b0e4a/USER'],
    '087b_old': ['MuEG2011Av4','/MuEG/mangano-R42X_S1_V05_S2_V02_S3_V02_ID087b_MuEG2011Av4_S2-bb1d7d67c682eb6954dd2da78c9b0e4a/USER'],

    '089_old': ['SingleElectron2011AMay10','/SingleElectron/mangano-R42X_S1_V05_S2_V02_S3_V02_ID089_SingleElectron2011AMay10_S2-v3-271b442b03cdae3729dc1e24bf59e029/USER'],
    '090_old': ['SingleMuon2011AMay10','/SingleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID090_SingleMu2011AReRecoMay10_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],
    '091_old': ['DoubleMuon2011AMay10','/DoubleMu/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID091_DoubleMu2011AReRecoMay10_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],
    '092_old': ['DoubleElectron2011AMay10','/DoubleElectron/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID092_DoubleElectron2011AReRecoMay10_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],
    '093_old': ['MuEG2011AMay10','/MuEG/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID093_MuEG2011AReRecoMay10_S2-b7e00ff49f7a03f5c5361f6cf358a673/USER'],
}

allBackgroundSamples = {
    '001': ['WZtoAny',0.0086320730,'/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID001_WZtoAny_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '002': ['ZZtoAny',0.0035141670,'/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID002_ZZtoAny_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '003': ['DYtoElEl',0.8362295180,'/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID003_DYtoElEl_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '004': ['DYtoMuMu',0.8396525670,'/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID004_DYtoMuMu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '005': ['DYtoTauTau',0.8349332880,'/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID005_DYtoTauTau_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '017B': ['DY10toElElZ2',1.7884117950,'/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID017B_DY10toElElZ2_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '018B': ['DY10toMuMuZ2',1.5517272340,'/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID018B_DY10toMuMuZ2_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '023': ['TTJetsMad',0.1352851040,'/TTJets_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID023_TTJetsMad_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '010': ['sTtoBLNu',0.0030111100,'/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID010_sTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '011': ['tTtoBLNu',0.0432392680,'/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID011_tTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '012': ['tWTtoBLNu',0.0216584220,'/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID012_tWTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '014': ['VVJetsTo4L',0.0066194200,'/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID014_VVJetsTo4L_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '022': ['GluGluToWWTo4L',0.0013909850,'/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID022_GluGluToWWTo4L_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '021': ['PhotonVJets',0.1497015500,'/PhotonVJets_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID021_PhotonVJets_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '026': ['WJetsToLNuMad',2.0722688030,'/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID026_WJetsToLNuMad_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

diBosonSamples = {
    '001': ['WZtoAny',0.0086320730,'/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID001_WZtoAny_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '002': ['ZZtoAny',0.0035141670,'/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID002_ZZtoAny_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

wzSamples = {
    '001': ['WZtoAny',0.0086320730,'/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID001_WZtoAny_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

zzSamples = {
    '002': ['ZZtoAny',0.0035141670,'/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID002_ZZtoAny_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

dySamples = {
    '003': ['DYtoElEl',0.8362295180,'/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID003_DYtoElEl_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '004': ['DYtoMuMu',0.8396525670,'/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID004_DYtoMuMu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '005': ['DYtoTauTau',0.8349332880,'/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID005_DYtoTauTau_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '017B': ['DY10toElElZ2',1.7884117950,'/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID017B_DY10toElElZ2_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '018B': ['DY10toMuMuZ2',1.5517272340,'/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID018B_DY10toMuMuZ2_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

topSamples = {
    '023': ['TTJetsMad',0.1352851040,'/TTJets_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID023_TTJetsMad_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '010': ['sTtoBLNu',0.0030111100,'/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID010_sTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '011': ['tTtoBLNu',0.0432392680,'/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID011_tTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '012': ['tWTtoBLNu',0.0216584220,'/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID012_tWTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

ttbarSamples = {
    '023': ['TTJetsMad',0.1352851040,'/TTJets_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID023_TTJetsMad_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

tWSamples = {
    '010': ['sTtoBLNu',0.0030111100,'/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID010_sTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '011': ['tTtoBLNu',0.0432392680,'/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID011_tTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '012': ['tWTtoBLNu',0.0216584220,'/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID012_tWTtoBLNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

wwSamples = {
    '014': ['VVJetsTo4L',0.0066194200,'/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID014_VVJetsTo4L_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '022': ['GluGluToWWTo4L',0.0013909850,'/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID022_GluGluToWWTo4L_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

wwSample = {
    '014': ['VVJetsTo4L',0.0066194200,'/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID014_VVJetsTo4L_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

ggWWSample = {
    '022': ['GluGluToWWTo4L',0.0013909850,'/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID022_GluGluToWWTo4L_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

wJetSamples = {
    '026': ['WJetsToLNuMad',2.0722688030,'/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID026_WJetsToLNuMad_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

wGamSamples = {
    '021': ['PhotonVJets',0.1497015500,'/PhotonVJets_7TeV-madgraph/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID021_PhotonVJets_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

dyMadSamples = {
    '027': ['DY10toLLMadD6T',67.4720304000,'/DYJetsToLL_TuneD6T_M-10To50_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID027_DY10toLLMadD6T_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '028': ['DY50toLLMadD6T',1.1982516850,'/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID028_DY50toLLMadD6T_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

otherBDTSamples = {
    '007': ['DYtoEEPythia','0.808742968'],
    '008': ['DYtoMuMuPythia','0.74218181'],
    '009': ['DYtoTauTauPythia','0.809741787'],
    '015': ['WWto2L','0.041036073'],
    '016': ['WWto2LMcnlo','0.050165233'],
    '016A': ['WWto2LMcnloUp','0.050178617'],
    '016B': ['WWto2LMcnloDown','0.050166906'],
    '029': ['WtoENu','2.055957174'],
    '030': ['WtoMuNu','2.071590747'],
    '031': ['WtoTauNu','2.186540979'],
    '032': ['WPlusENu','3.079151512'],
    '033': ['WMinusENu','2.201636176'],
    '034': ['WPlusTauNu','3.08236354'],
    '035': ['WMinusTauNu','2.156279447'],
    '036': ['WPlusMuNu','3.102470276'],
    '037': ['WMinusMuNu','2.166189054'],
}

h120Samples = {
    '101120': ['ggToH120toWWto2L2Nu',0.0010517640,'/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101120_ggToH120toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102120': ['ggToH120toWWtoLNuTauNu',0.0010091740,'/GluGluToHToWWToLNuTauNu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102120_ggToH120toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103120': ['ggToH120toWWto2Tau2Nu',0.0004197220,'/GluGluToHToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103120_ggToH120toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104120': ['vbfToH120toWWto2L2Nu',0.0000764742,'/VBF_HToWWTo2L2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104120_vbfToH120toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105120': ['vbfToH120toWWtoLNuTauNu',0.0000764777,'/VBF_HToWWToLNuTauNu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105120_vbfToH120toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106120': ['vbfToH120toWWto2Tau2Nu',0.0000318728,'/VBF_HToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106120_vbfToH120toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h130Samples = {
    '101130': ['ggToH130toWWto2L2Nu',0.0018282650,'/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101130_ggToH130toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102130': ['ggToH130toWWtoLNuTauNu',0.0018281980,'/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102130_ggToH130toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103130': ['ggToH130toWWto2Tau2Nu',0.0007621790,'/GluGluToHToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103130_ggToH130toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104130': ['vbfToH130toWWto2L2Nu',0.0001483670,'/VBF_HToWWTo2L2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104130_vbfToH130toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105130': ['vbfToH130toWWtoLNuTauNu',0.0001483760,'/VBF_HToWWToLNuTauNu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105130_vbfToH130toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106130': ['vbfToH130toWWto2Tau2Nu',0.0000621994,'/VBF_HToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106130_vbfToH130toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h140Samples = {
    '101140': ['ggToH140toWWto2L2Nu',0.0026542580,'/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101140_ggToH140toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102140': ['ggToH140toWWtoLNuTauNu',0.0025956910,'/GluGluToHToWWToLNuTauNu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102140_ggToH140toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103140': ['ggToH140toWWto2Tau2Nu',0.0010819170,'/GluGluToHToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103140_ggToH140toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104140': ['vbfToH140toWWto2L2Nu',0.0002287990,'/VBF_HToWWTo2L2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104140_vbfToH140toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105140': ['vbfToH140toWWtoLNuTauNu',0.0002288390,'/VBF_HToWWToLNuTauNu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105140_vbfToH140toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106140': ['vbfToH140toWWto2Tau2Nu',0.0000925687,'/VBF_HToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106140_vbfToH140toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h150Samples = {
    '101150': ['ggToH150toWWto2L2Nu',0.0031658890,'/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101150_ggToH150toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102150': ['ggToH150toWWtoLNuTauNu',0.0031194030,'/GluGluToHToWWToLNuTauNu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102150_ggToH150toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103150': ['ggToH150toWWto2Tau2Nu',0.0013000790,'/GluGluToHToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103150_ggToH150toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104150': ['vbfToH150toWWto2L2Nu',0.0002821950,'/VBF_HToWWTo2L2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104150_vbfToH150toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105150': ['vbfToH150toWWtoLNuTauNu',0.0002822360,'/VBF_HToWWToLNuTauNu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105150_vbfToH150toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106150': ['vbfToH150toWWto2Tau2Nu',0.0001183160,'/VBF_HToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106150_vbfToH150toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h160Samples = {
    '101160': ['ggToH160toWWto2L2Nu',0.0036210630,'/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101160_ggToH160toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102160': ['ggToH160toWWtoLNuTauNu',0.0035009730,'/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102160_ggToH160toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103160': ['ggToH160toWWto2Tau2Nu',0.0014591570,'/GluGluToHToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103160_ggToH160toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104160': ['vbfToH160toWWto2L2Nu',0.0003359000,'/VBF_HToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104160_vbfToH160toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105160': ['vbfToH160toWWtoLNuTauNu',0.0003358880,'/VBF_HToWWToLNuTauNu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105160_vbfToH160toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106160': ['vbfToH160toWWto2Tau2Nu',0.0001395820,'/VBF_HToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106160_vbfToH160toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h170Samples = {
    '101170': ['ggToH170toWWto2L2Nu',0.0032730990,'/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101170_ggToH170toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102170': ['ggToH170toWWtoLNuTauNu',0.0031766810,'/GluGluToHToWWToLNuTauNu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102170_ggToH170toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103170': ['ggToH170toWWto2Tau2Nu',0.0013228280,'/GluGluToHToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103170_ggToH170toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104170': ['vbfToH170toWWto2L2Nu',0.0003344760,'/VBF_HToWWTo2L2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104170_vbfToH170toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105170': ['vbfToH170toWWtoLNuTauNu',0.0003489970,'/VBF_HToWWToLNuTauNu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105170_vbfToH170toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106170': ['vbfToH170toWWto2Tau2Nu',0.0001380940,'/VBF_HToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106170_vbfToH170toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h180Samples = {
    '101180': ['ggToH180toWWto2L2Nu',0.0026918110,'/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101180_ggToH180toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102180': ['ggToH180toWWtoLNuTauNu',0.0026722070,'/GluGluToHToWWToLNuTauNu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102180_ggToH180toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103180': ['ggToH180toWWto2Tau2Nu',0.0011137210,'/GluGluToHToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103180_ggToH180toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104180': ['vbfToH180toWWto2L2Nu',0.0002949390,'/VBF_HToWWTo2L2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104180_vbfToH180toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105180': ['vbfToH180toWWtoLNuTauNu',0.0003007590,'/VBF_HToWWToLNuTauNu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105180_vbfToH180toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106180': ['vbfToH180toWWto2Tau2Nu',0.0001213570,'/VBF_HToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106180_vbfToH180toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h190Samples = {
    '101190': ['ggToH190toWWto2L2Nu',0.0020183710,'/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101190_ggToH190toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102190': ['ggToH190toWWtoLNuTauNu',0.0019737620,'/GluGluToHToWWToLNuTauNu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102190_ggToH190toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103190': ['ggToH190toWWto2Tau2Nu',0.0008227900,'/GluGluToHToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103190_ggToH190toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104190': ['vbfToH190toWWto2L2Nu',0.0002291990,'/VBF_HToWWTo2L2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104190_vbfToH190toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105190': ['vbfToH190toWWtoLNuTauNu',0.0002343980,'/VBF_HToWWToLNuTauNu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105190_vbfToH190toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106190': ['vbfToH190toWWto2Tau2Nu',0.0000955878,'/VBF_HToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106190_vbfToH190toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h200Samples = {
    '101200': ['ggToH200toWWto2L2Nu',0.0016565440,'/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101200_ggToH200toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102200': ['ggToH200toWWtoLNuTauNu',0.0027607320,'/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102200_ggToH200toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103200': ['ggToH200toWWto2Tau2Nu',0.0006894570,'/GluGluToHToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103200_ggToH200toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104200': ['vbfToH200toWWto2L2Nu',0.0001984470,'/VBF_HToWWTo2L2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104200_vbfToH200toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105200': ['vbfToH200toWWtoLNuTauNu',0.0001986240,'/VBF_HToWWToLNuTauNu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105200_vbfToH200toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106200': ['vbfToH200toWWto2Tau2Nu',0.0000819461,'/VBF_HToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106200_vbfToH200toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h250Samples = {
    '101250': ['ggToH250toWWto2L2Nu',0.0009873000,'/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101250_ggToH250toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102250': ['ggToH250toWWtoLNuTauNu',0.0009873450,'/GluGluToHToWWToLNuTauNu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102250_ggToH250toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103250': ['ggToH250toWWto2Tau2Nu',0.0004106750,'/GluGluToHToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103250_ggToH250toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104250': ['vbfToH250toWWto2L2Nu',0.0001265270,'/VBF_HToWWTo2L2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104250_vbfToH250toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105250': ['vbfToH250toWWtoLNuTauNu',0.0001267470,'/VBF_HToWWToLNuTauNu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105250_vbfToH250toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106250': ['vbfToH250toWWto2Tau2Nu',0.0000531051,'/VBF_HToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106250_vbfToH250toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h300Samples = {
    '101300': ['ggToH300toWWto2L2Nu',0.0007100390,'/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101300_ggToH300toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102300': ['ggToH300toWWtoLNuTauNu',0.0008154700,'/GluGluToHToWWToLNuTauNu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102300_ggToH300toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103300': ['ggToH300toWWto2Tau2Nu',0.0002954810,'/GluGluToHToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103300_ggToH300toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104300': ['vbfToH300toWWto2L2Nu',0.0000877618,'/VBF_HToWWTo2L2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104300_vbfToH300toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105300': ['vbfToH300toWWtoLNuTauNu',0.0000876088,'/VBF_HToWWToLNuTauNu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105300_vbfToH300toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106300': ['vbfToH300toWWto2Tau2Nu',0.0000364515,'/VBF_HToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106300_vbfToH300toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h350Samples = {
    '101350': ['ggToH350toWWto2L2Nu',0.0006537020,'/GluGluToHToWWTo2L2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101350_ggToH350toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102350': ['ggToH350toWWtoLNuTauNu',0.0006536960,'/GluGluToHToWWToLNuTauNu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102350_ggToH350toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103350': ['ggToH350toWWto2Tau2Nu',0.0002727690,'/GluGluToHToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103350_ggToH350toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104350': ['vbfToH350toWWto2L2Nu',0.0000613519,'/VBF_HToWWTo2L2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104350_vbfToH350toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105350': ['vbfToH350toWWtoLNuTauNu',0.0000612616,'/VBF_HToWWToLNuTauNu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105350_vbfToH350toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106350': ['vbfToH350toWWto2Tau2Nu',0.0000258374,'/VBF_HToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106350_vbfToH350toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h400Samples = {
    '101400': ['ggToH400toWWto2L2Nu',0.0005010280,'/GluGluToHToWWTo2L2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101400_ggToH400toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102400': ['ggToH400toWWtoLNuTauNu',0.0005010090,'/GluGluToHToWWToLNuTauNu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102400_ggToH400toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103400': ['ggToH400toWWto2Tau2Nu',0.0002091190,'/GluGluToHToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103400_ggToH400toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104400': ['vbfToH400toWWto2L2Nu',0.0000403211,'/VBF_HToWWTo2L2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104400_vbfToH400toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105400': ['vbfToH400toWWtoLNuTauNu',0.0000439065,'/VBF_HToWWToLNuTauNu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105400_vbfToH400toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106400': ['vbfToH400toWWto2Tau2Nu',0.0000167418,'/VBF_HToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106400_vbfToH400toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h450Samples = {
    '101450': ['ggToH450toWWto2L2Nu',0.0003200700,'/GluGluToHToWWTo2L2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101450_ggToH450toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102450': ['ggToH450toWWtoLNuTauNu',0.0003200580,'/GluGluToHToWWToLNuTauNu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102450_ggToH450toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103450': ['ggToH450toWWto2Tau2Nu',0.0001333620,'/GluGluToHToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103450_ggToH450toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104450': ['vbfToH450toWWto2L2Nu',0.0000292633,'/VBF_HToWWTo2L2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104450_vbfToH450toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105450': ['vbfToH450toWWtoLNuTauNu',0.0000292593,'/VBF_HToWWToLNuTauNu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105450_vbfToH450toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106450': ['vbfToH450toWWto2Tau2Nu',0.0000121895,'/VBF_HToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106450_vbfToH450toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h500Samples = {
    '101500': ['ggToH500toWWto2L2Nu',0.0002000490,'/GluGluToHToWWTo2L2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101500_ggToH500toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102500': ['ggToH500toWWtoLNuTauNu',0.0002000330,'/GluGluToHToWWToLNuTauNu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102500_ggToH500toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103500': ['ggToH500toWWto2Tau2Nu',0.0000833498,'/GluGluToHToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103500_ggToH500toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104500': ['vbfToH500toWWto2L2Nu',0.0000225270,'/VBF_HToWWTo2L2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104500_vbfToH500toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105500': ['vbfToH500toWWtoLNuTauNu',0.0000239490,'/VBF_HToWWToLNuTauNu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105500_vbfToH500toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106500': ['vbfToH500toWWto2Tau2Nu',0.0000091398,'/VBF_HToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106500_vbfToH500toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h550Samples = {
    '101550': ['ggToH550toWWto2L2Nu',0.0001254820,'/GluGluToHToWWTo2L2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101550_ggToH550toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102550': ['ggToH550toWWtoLNuTauNu',0.0001282470,'/GluGluToHToWWToLNuTauNu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102550_ggToH550toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103550': ['ggToH550toWWto2Tau2Nu',0.0000530440,'/GluGluToHToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103550_ggToH550toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104550': ['vbfToH550toWWto2L2Nu',0.0000178338,'/VBF_HToWWTo2L2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104550_vbfToH550toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105550': ['vbfToH550toWWtoLNuTauNu',0.0000189596,'/VBF_HToWWToLNuTauNu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105550_vbfToH550toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106550': ['vbfToH550toWWto2Tau2Nu',0.0000076165,'/VBF_HToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106550_vbfToH550toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

h600Samples = {
    '101600': ['ggToH600toWWto2L2Nu',0.0000791103,'/GluGluToHToWWTo2L2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101600_ggToH600toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102600': ['ggToH600toWWtoLNuTauNu',0.0000808513,'/GluGluToHToWWToLNuTauNu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102600_ggToH600toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103600': ['ggToH600toWWto2Tau2Nu',0.0000333450,'/GluGluToHToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103600_ggToH600toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104600': ['vbfToH600toWWto2L2Nu',0.0000140794,'/VBF_HToWWTo2L2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104600_vbfToH600toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105600': ['vbfToH600toWWtoLNuTauNu',0.0000149681,'/VBF_HToWWToLNuTauNu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105600_vbfToH600toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106600': ['vbfToH600toWWto2Tau2Nu',0.0000061034,'/VBF_HToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106600_vbfToH600toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

allSignalSamples = {
    '101120': ['ggToH120toWWto2L2Nu',0.0010517640,'/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101120_ggToH120toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102120': ['ggToH120toWWtoLNuTauNu',0.0010091740,'/GluGluToHToWWToLNuTauNu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102120_ggToH120toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103120': ['ggToH120toWWto2Tau2Nu',0.0004197220,'/GluGluToHToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103120_ggToH120toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104120': ['vbfToH120toWWto2L2Nu',0.0000764742,'/VBF_HToWWTo2L2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104120_vbfToH120toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105120': ['vbfToH120toWWtoLNuTauNu',0.0000764777,'/VBF_HToWWToLNuTauNu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105120_vbfToH120toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106120': ['vbfToH120toWWto2Tau2Nu',0.0000318728,'/VBF_HToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106120_vbfToH120toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101130': ['ggToH130toWWto2L2Nu',0.0018282650,'/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101130_ggToH130toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102130': ['ggToH130toWWtoLNuTauNu',0.0018281980,'/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102130_ggToH130toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103130': ['ggToH130toWWto2Tau2Nu',0.0007621790,'/GluGluToHToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103130_ggToH130toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104130': ['vbfToH130toWWto2L2Nu',0.0001483670,'/VBF_HToWWTo2L2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104130_vbfToH130toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105130': ['vbfToH130toWWtoLNuTauNu',0.0001483760,'/VBF_HToWWToLNuTauNu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105130_vbfToH130toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106130': ['vbfToH130toWWto2Tau2Nu',0.0000621994,'/VBF_HToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106130_vbfToH130toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101140': ['ggToH140toWWto2L2Nu',0.0026542580,'/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101140_ggToH140toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102140': ['ggToH140toWWtoLNuTauNu',0.0025956910,'/GluGluToHToWWToLNuTauNu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102140_ggToH140toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103140': ['ggToH140toWWto2Tau2Nu',0.0010819170,'/GluGluToHToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103140_ggToH140toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104140': ['vbfToH140toWWto2L2Nu',0.0002287990,'/VBF_HToWWTo2L2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104140_vbfToH140toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105140': ['vbfToH140toWWtoLNuTauNu',0.0002288390,'/VBF_HToWWToLNuTauNu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105140_vbfToH140toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106140': ['vbfToH140toWWto2Tau2Nu',0.0000925687,'/VBF_HToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106140_vbfToH140toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101150': ['ggToH150toWWto2L2Nu',0.0031658890,'/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101150_ggToH150toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102150': ['ggToH150toWWtoLNuTauNu',0.0031194030,'/GluGluToHToWWToLNuTauNu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102150_ggToH150toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103150': ['ggToH150toWWto2Tau2Nu',0.0013000790,'/GluGluToHToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103150_ggToH150toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104150': ['vbfToH150toWWto2L2Nu',0.0002821950,'/VBF_HToWWTo2L2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104150_vbfToH150toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105150': ['vbfToH150toWWtoLNuTauNu',0.0002822360,'/VBF_HToWWToLNuTauNu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105150_vbfToH150toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106150': ['vbfToH150toWWto2Tau2Nu',0.0001183160,'/VBF_HToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106150_vbfToH150toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101160': ['ggToH160toWWto2L2Nu',0.0036210630,'/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101160_ggToH160toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102160': ['ggToH160toWWtoLNuTauNu',0.0035009730,'/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102160_ggToH160toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103160': ['ggToH160toWWto2Tau2Nu',0.0014591570,'/GluGluToHToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103160_ggToH160toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104160': ['vbfToH160toWWto2L2Nu',0.0003359000,'/VBF_HToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104160_vbfToH160toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105160': ['vbfToH160toWWtoLNuTauNu',0.0003358880,'/VBF_HToWWToLNuTauNu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105160_vbfToH160toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106160': ['vbfToH160toWWto2Tau2Nu',0.0001395820,'/VBF_HToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106160_vbfToH160toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101170': ['ggToH170toWWto2L2Nu',0.0032730990,'/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101170_ggToH170toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102170': ['ggToH170toWWtoLNuTauNu',0.0031766810,'/GluGluToHToWWToLNuTauNu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102170_ggToH170toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103170': ['ggToH170toWWto2Tau2Nu',0.0013228280,'/GluGluToHToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103170_ggToH170toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104170': ['vbfToH170toWWto2L2Nu',0.0003344760,'/VBF_HToWWTo2L2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104170_vbfToH170toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105170': ['vbfToH170toWWtoLNuTauNu',0.0003489970,'/VBF_HToWWToLNuTauNu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105170_vbfToH170toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106170': ['vbfToH170toWWto2Tau2Nu',0.0001380940,'/VBF_HToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106170_vbfToH170toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101180': ['ggToH180toWWto2L2Nu',0.0026918110,'/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101180_ggToH180toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102180': ['ggToH180toWWtoLNuTauNu',0.0026722070,'/GluGluToHToWWToLNuTauNu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102180_ggToH180toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103180': ['ggToH180toWWto2Tau2Nu',0.0011137210,'/GluGluToHToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103180_ggToH180toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104180': ['vbfToH180toWWto2L2Nu',0.0002949390,'/VBF_HToWWTo2L2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104180_vbfToH180toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105180': ['vbfToH180toWWtoLNuTauNu',0.0003007590,'/VBF_HToWWToLNuTauNu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105180_vbfToH180toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106180': ['vbfToH180toWWto2Tau2Nu',0.0001213570,'/VBF_HToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106180_vbfToH180toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101190': ['ggToH190toWWto2L2Nu',0.0020183710,'/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101190_ggToH190toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102190': ['ggToH190toWWtoLNuTauNu',0.0019737620,'/GluGluToHToWWToLNuTauNu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102190_ggToH190toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103190': ['ggToH190toWWto2Tau2Nu',0.0008227900,'/GluGluToHToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103190_ggToH190toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104190': ['vbfToH190toWWto2L2Nu',0.0002291990,'/VBF_HToWWTo2L2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104190_vbfToH190toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105190': ['vbfToH190toWWtoLNuTauNu',0.0002343980,'/VBF_HToWWToLNuTauNu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105190_vbfToH190toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106190': ['vbfToH190toWWto2Tau2Nu',0.0000955878,'/VBF_HToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106190_vbfToH190toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101200': ['ggToH200toWWto2L2Nu',0.0016565440,'/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101200_ggToH200toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102200': ['ggToH200toWWtoLNuTauNu',0.0027607320,'/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102200_ggToH200toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103200': ['ggToH200toWWto2Tau2Nu',0.0006894570,'/GluGluToHToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103200_ggToH200toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104200': ['vbfToH200toWWto2L2Nu',0.0001984470,'/VBF_HToWWTo2L2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104200_vbfToH200toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105200': ['vbfToH200toWWtoLNuTauNu',0.0001986240,'/VBF_HToWWToLNuTauNu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105200_vbfToH200toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106200': ['vbfToH200toWWto2Tau2Nu',0.0000819461,'/VBF_HToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106200_vbfToH200toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101250': ['ggToH250toWWto2L2Nu',0.0009873000,'/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101250_ggToH250toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102250': ['ggToH250toWWtoLNuTauNu',0.0009873450,'/GluGluToHToWWToLNuTauNu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102250_ggToH250toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103250': ['ggToH250toWWto2Tau2Nu',0.0004106750,'/GluGluToHToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103250_ggToH250toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104250': ['vbfToH250toWWto2L2Nu',0.0001265270,'/VBF_HToWWTo2L2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104250_vbfToH250toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105250': ['vbfToH250toWWtoLNuTauNu',0.0001267470,'/VBF_HToWWToLNuTauNu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105250_vbfToH250toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106250': ['vbfToH250toWWto2Tau2Nu',0.0000531051,'/VBF_HToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106250_vbfToH250toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101300': ['ggToH300toWWto2L2Nu',0.0007100390,'/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101300_ggToH300toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102300': ['ggToH300toWWtoLNuTauNu',0.0008154700,'/GluGluToHToWWToLNuTauNu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102300_ggToH300toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103300': ['ggToH300toWWto2Tau2Nu',0.0002954810,'/GluGluToHToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103300_ggToH300toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104300': ['vbfToH300toWWto2L2Nu',0.0000877618,'/VBF_HToWWTo2L2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104300_vbfToH300toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105300': ['vbfToH300toWWtoLNuTauNu',0.0000876088,'/VBF_HToWWToLNuTauNu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105300_vbfToH300toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106300': ['vbfToH300toWWto2Tau2Nu',0.0000364515,'/VBF_HToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106300_vbfToH300toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101350': ['ggToH350toWWto2L2Nu',0.0006537020,'/GluGluToHToWWTo2L2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101350_ggToH350toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102350': ['ggToH350toWWtoLNuTauNu',0.0006536960,'/GluGluToHToWWToLNuTauNu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102350_ggToH350toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103350': ['ggToH350toWWto2Tau2Nu',0.0002727690,'/GluGluToHToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103350_ggToH350toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104350': ['vbfToH350toWWto2L2Nu',0.0000613519,'/VBF_HToWWTo2L2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104350_vbfToH350toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105350': ['vbfToH350toWWtoLNuTauNu',0.0000612616,'/VBF_HToWWToLNuTauNu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105350_vbfToH350toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106350': ['vbfToH350toWWto2Tau2Nu',0.0000258374,'/VBF_HToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106350_vbfToH350toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101400': ['ggToH400toWWto2L2Nu',0.0005010280,'/GluGluToHToWWTo2L2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101400_ggToH400toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102400': ['ggToH400toWWtoLNuTauNu',0.0005010090,'/GluGluToHToWWToLNuTauNu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102400_ggToH400toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103400': ['ggToH400toWWto2Tau2Nu',0.0002091190,'/GluGluToHToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103400_ggToH400toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104400': ['vbfToH400toWWto2L2Nu',0.0000403211,'/VBF_HToWWTo2L2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104400_vbfToH400toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105400': ['vbfToH400toWWtoLNuTauNu',0.0000439065,'/VBF_HToWWToLNuTauNu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105400_vbfToH400toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106400': ['vbfToH400toWWto2Tau2Nu',0.0000167418,'/VBF_HToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106400_vbfToH400toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101450': ['ggToH450toWWto2L2Nu',0.0003200700,'/GluGluToHToWWTo2L2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101450_ggToH450toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102450': ['ggToH450toWWtoLNuTauNu',0.0003200580,'/GluGluToHToWWToLNuTauNu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102450_ggToH450toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103450': ['ggToH450toWWto2Tau2Nu',0.0001333620,'/GluGluToHToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103450_ggToH450toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104450': ['vbfToH450toWWto2L2Nu',0.0000292633,'/VBF_HToWWTo2L2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104450_vbfToH450toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105450': ['vbfToH450toWWtoLNuTauNu',0.0000292593,'/VBF_HToWWToLNuTauNu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105450_vbfToH450toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106450': ['vbfToH450toWWto2Tau2Nu',0.0000121895,'/VBF_HToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106450_vbfToH450toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101500': ['ggToH500toWWto2L2Nu',0.0002000490,'/GluGluToHToWWTo2L2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101500_ggToH500toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102500': ['ggToH500toWWtoLNuTauNu',0.0002000330,'/GluGluToHToWWToLNuTauNu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102500_ggToH500toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103500': ['ggToH500toWWto2Tau2Nu',0.0000833498,'/GluGluToHToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103500_ggToH500toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104500': ['vbfToH500toWWto2L2Nu',0.0000225270,'/VBF_HToWWTo2L2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104500_vbfToH500toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105500': ['vbfToH500toWWtoLNuTauNu',0.0000239490,'/VBF_HToWWToLNuTauNu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105500_vbfToH500toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106500': ['vbfToH500toWWto2Tau2Nu',0.0000091398,'/VBF_HToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106500_vbfToH500toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101550': ['ggToH550toWWto2L2Nu',0.0001254820,'/GluGluToHToWWTo2L2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101550_ggToH550toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102550': ['ggToH550toWWtoLNuTauNu',0.0001282470,'/GluGluToHToWWToLNuTauNu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102550_ggToH550toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103550': ['ggToH550toWWto2Tau2Nu',0.0000530440,'/GluGluToHToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103550_ggToH550toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104550': ['vbfToH550toWWto2L2Nu',0.0000178338,'/VBF_HToWWTo2L2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104550_vbfToH550toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105550': ['vbfToH550toWWtoLNuTauNu',0.0000189596,'/VBF_HToWWToLNuTauNu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105550_vbfToH550toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106550': ['vbfToH550toWWto2Tau2Nu',0.0000076165,'/VBF_HToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106550_vbfToH550toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '101600': ['ggToH600toWWto2L2Nu',0.0000791103,'/GluGluToHToWWTo2L2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID101600_ggToH600toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '102600': ['ggToH600toWWtoLNuTauNu',0.0000808513,'/GluGluToHToWWToLNuTauNu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID102600_ggToH600toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '103600': ['ggToH600toWWto2Tau2Nu',0.0000333450,'/GluGluToHToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID103600_ggToH600toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '104600': ['vbfToH600toWWto2L2Nu',0.0000140794,'/VBF_HToWWTo2L2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID104600_vbfToH600toWWto2L2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '105600': ['vbfToH600toWWtoLNuTauNu',0.0000149681,'/VBF_HToWWToLNuTauNu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID105600_vbfToH600toWWtoLNuTauNu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
    '106600': ['vbfToH600toWWto2Tau2Nu',0.0000061034,'/VBF_HToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V05_S2_V02_S3_V02_ID106600_vbfToH600toWWto2Tau2Nu_S2-e08cc3ddcace102ccbd78f8ef4522d96/USER'],
}

