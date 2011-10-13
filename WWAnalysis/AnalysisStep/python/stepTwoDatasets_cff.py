import FWCore.ParameterSet.Config as cms
# to remove lines with no dataset:
# %s/^[0-9]\+ *$//
# convert rest of the lines to dictionaries
# :%s/^\([0-9]\+\) *\(.*\)$/    '\1':'\2',/


stepTwoDatasets = {
    '000':'/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID000_WWJets2LMad-a0fee9021ee0e974e473d874349bbb62/USER',
    '001':'/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID001_GluGluToWWTo4L-a0fee9021ee0e974e473d874349bbb62/USER',
    '002':'/WWTo2L2Nu_CT10_7TeV-mcatnlo/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID002_WWto2L2NuMCatNLO-27316856f5d7d57332f8edbc6c095e34/USER',
    '003':'/WWTo2L2Nu_scaleup_CT10_7TeV-mcatnlo/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID003_WWto2L2NuMCatNLOUp-27316856f5d7d57332f8edbc6c095e34/USER',
    '004':'/WWTo2L2Nu_scaledown_CT10_7TeV-mcatnlo/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID004_WWto2L2NuMCatNLODown-27316856f5d7d57332f8edbc6c095e34/USER',
            
    '010':'/TTJets_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID010_TTJetsMad-27316856f5d7d57332f8edbc6c095e34/USER',
    '011':'/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID011_TtWFullDR-a0fee9021ee0e974e473d874349bbb62/USER',
    '012':'/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID012_TbartWFullDR-a0fee9021ee0e974e473d874349bbb62/USER',
    '013':'/T_TuneZ2_t-channel_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID013_TtFull-27316856f5d7d57332f8edbc6c095e34/USER',
    '014':'/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID014_TbartFull-a0fee9021ee0e974e473d874349bbb62/USER',
    '015':'/T_TuneZ2_s-channel_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID015_TsFull-a0fee9021ee0e974e473d874349bbb62/USER',
    '016':'/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID016_TbarsFull-a0fee9021ee0e974e473d874349bbb62/USER',
    '017':'/T_TuneZ2_tW-channel-DS_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID017_TtWFullDS-27316856f5d7d57332f8edbc6c095e34/USER',
    '018':'/Tbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauola/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID018_TbartWFullDS-27316856f5d7d57332f8edbc6c095e34/USER',
    '019':'/TTTo2L2Nu2B_7TeV-powheg-pythia6/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID019_TTTo2L2Nu2B-3458aa0bd471577b09aa6b73035cd433/USER',
            
    '030':'/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID030_DYtoElEl-a0fee9021ee0e974e473d874349bbb62/USER',
    '031':'/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID031_DYtoMuMu-a0fee9021ee0e974e473d874349bbb62/USER',
    '032':'/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID032_DYtoTauTau-a0fee9021ee0e974e473d874349bbb62/USER',
    '033':'/DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID0033_DY10toElEl-27316856f5d7d57332f8edbc6c095e34/USER',
    '034':'/DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID0034_DY10toMuMu-27316856f5d7d57332f8edbc6c095e34/USER',
    '035':'/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID035_DY10toTauTau-a0fee9021ee0e974e473d874349bbb62/USER',

    '037':'/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID037_DY50toLLMadD6T-27316856f5d7d57332f8edbc6c095e34/USER',
    '038':'/ZGToEEG_TuneZ2_7TeV-madgraph/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID038_ZgammaToElElMad-27316856f5d7d57332f8edbc6c095e34/USER',
    '039':'/ZGToMuMuG_TuneZ2_7TeV-madgraph/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID039_ZgammaToMuMuMad-27316856f5d7d57332f8edbc6c095e34/USER',
    '040':'/ZGToTauTauG_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID040_ZgammaToTauTauMad-27316856f5d7d57332f8edbc6c095e34/USER',
    '041':'/ZGToNuNuG_TuneZ2_7TeV-madgraph/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID041_ZgammaToNuNuMad-27316856f5d7d57332f8edbc6c095e34/USER',
            
    '070':'/WZ_TuneZ2_7TeV_pythia6_tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID070_WZFull-a0fee9021ee0e974e473d874349bbb62/USER',
    '071':'/ZZ_TuneZ2_7TeV_pythia6_tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID071_ZZFull-a0fee9021ee0e974e473d874349bbb62/USER',
    '072':'/GluGluToZZTo2L2L_7TeV-gg2zz-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID072_GluGluZZ2L2L-27316856f5d7d57332f8edbc6c095e34/USER',
    '073':'/GluGluToZZTo4L_7TeV-gg2zz-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID073_GluGluZZ4L-a0fee9021ee0e974e473d874349bbb62/USER',
    '074':'/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID074_WZJetsMad-a90b596e1a39c995fc1403cf7cc2b14b/USER',
    '075':'/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID075_ZZJetsMad-a90b596e1a39c995fc1403cf7cc2b14b/USER',
        
    '080':'/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID080_WJetsToLNuMad-a0fee9021ee0e974e473d874349bbb62/USER',
    '081':'/GVJets_7TeV-madgraph/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID081_VGamma-a0fee9021ee0e974e473d874349bbb62/USER',
    '082':'/WGToENuG_TuneZ2_7TeV-madgraph/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID082_WgammaToElNuMad-27316856f5d7d57332f8edbc6c095e34/USER',
    '083':'/WGToMuNuG_TuneZ2_7TeV-madgraph/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID083_WgammaToMuNuMad-27316856f5d7d57332f8edbc6c095e34/USER',
    '084':'/WGToTauNuG_TuneZ2_7TeV-madgraph-tauola/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID084_WgammaToTauNuMad-27316856f5d7d57332f8edbc6c095e34/USER',
            
    '100':'/SingleElectron/mangano-R42X_S1_V06_S2_V04_S3_V06_ID100_SingleElectron2011Av4_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '101':'/SingleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID101_SingleMuon2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '102':'/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID102_DoubleElectron2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '103':'/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID103_DoubleMuon2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '104':'/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID104_MuEG2011Av4-eff9861cde6b1617ce28065ed0ad5d45/USER',

            





            
    '120':'/SingleElectron/mangano-R42X_S1_V06_S2_V04_S3_V06_ID120_SingleElectron2011Av6_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '121':'/SingleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID121_SingleMuon2011Av6-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '122':'/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID122_DoubleElectron2011Av6-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '123':'/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID123_DoubleMuon2011Av6-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '124':'/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID124_MuEG2011Av6v3-eff9861cde6b1617ce28065ed0ad5d45/USER',
            
    '130':'/SingleElectron/mangano-R42X_S1_V06_S2_V04_S3_V06_ID130_SingleElectron2011Av6v2_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '131':'/SingleMu/mangano-R42X_S1_V06_S2_V04_S3_V06_ID131_SingleMuon2011Av6v2_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '132':'/DoubleElectron/mangano-R42X_S1_V06_S2_V04_S3_V06_ID132_DoubleElectron2011Av6v2_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '133':'/DoubleMu/mangano-R42X_S1_V06_S2_V04_S3_V06_ID133_DoubleMuon2011Av6v2_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '134':'/MuEG/mangano-R42X_S1_V06_S2_V04_S3_V06_ID134_MuEG2011Av6v2_S2-9aa2bd39036dba1ed10c175bf492df72/USER',


    '140a':'/SingleElectron/mangano-R42X_S1_V06_S2_V02_S3_V05_ID130_SingleElectron2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '141a':'/SingleMu/mangano-R42X_S1_V06_S2_V02_S3_V05_ID131_SingleMuon2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '142a':'/DoubleElectron/mangano-R42X_S1_V06_S2_V02_S3_V05_ID132_DoubleElectron2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '143a':'/DoubleMu/mangano-R42X_S1_V06_S2_V02_S3_V05_ID133_DoubleMuon2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '144a':'/MuEG/mangano-R42X_S1_V06_S2_V02_S3_V05_ID144_MuEG2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',

    '140b':'/SingleElectron/mangano-R42X_S1_V06_S2_V02_S3_V05_ID140b_SingleElectron2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '141b':'/SingleMu/mangano-R42X_S1_V06_S2_V02_S3_V05_ID141b_SingleMuon2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '142b':'/DoubleElectron/mangano-R42X_S1_V06_S2_V02_S3_V05_ID142b_DoubleElectron2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '143b':'/DoubleMu/mangano-R42X_S1_V06_S2_V02_S3_V05_ID143b_DoubleMuon2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',
    '144b':'/MuEG/mangano-R42X_S1_V06_S2_V02_S3_V05_ID144b_MuEG2011Bv1a-ddde5080c9ce1099b0b1ba7bc18443f1/USER',

    '140c':'/SingleElectron/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID140c_SingleElectron2011Bv1a-29f699e44c6b7f711ce3481cf319381b/USER',
    '141c':'/SingleMu/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID141c_SingleMuon2011Bv1a-29f699e44c6b7f711ce3481cf319381b/USER',
    '142c':'/DoubleElectron/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID142c_DoubleElectron2011Bv1a-29f699e44c6b7f711ce3481cf319381b/USER',
    '143c':'/DoubleMu/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID143c_DoubleMuon2011Bv1a-29f699e44c6b7f711ce3481cf319381b/USER',
    '144c':'/MuEG/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID144c_MuEG2011Bv1a-29f699e44c6b7f711ce3481cf319381b/USER',

    '140d':'/SingleElectron/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID140d_SingleElectron2011Bv1d-29f699e44c6b7f711ce3481cf319381b/USER',
    '141d':'/SingleMu/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID141d_SingleMuon2011Bv1d-29f699e44c6b7f711ce3481cf319381b/USER',
    '142d':'/DoubleElectron/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID142d_DoubleElectron2011Bv1d2-29f699e44c6b7f711ce3481cf319381b/USER',
    '143d':'/DoubleMu/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID143d_DoubleMuon2011Bv1d-29f699e44c6b7f711ce3481cf319381b/USER',
    '144d':'/MuEG/jfernan2-R42X_S1_V06_S2_V02_S3_V05_ID144d_MuEG2011Bv1d-29f699e44c6b7f711ce3481cf319381b/USER',            
            
    '150':'/SingleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID150_SingleElectron2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '151':'/SingleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID151_SingleMuon2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '152':'/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID152_DoubleMuon2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '153':'/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID153_DoubleElectron2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '154':'/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID154_MuEG2011AMay10-eff9861cde6b1617ce28065ed0ad5d45/USER',

            
    '160':'/SingleElectron/mangano-R42X_S1_V06_S2_V04_S3_V06_ID160_SingleElectron2011AAug05_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '161':'/SingleMu/mangano-R42X_S1_V06_S2_V04_S3_V06_ID161_SingleMuon2011AAug05_S2-9aa2bd39036dba1ed10c175bf492df72/USER',
    '162':'/DoubleElectron/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID162_DoubleElectron2011AAug05v2-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '163':'/DoubleMu/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID163_DoubleMu2011AAug05-eff9861cde6b1617ce28065ed0ad5d45/USER',
    '164':'/MuEG/mwlebour-R42X_S1_V06_S2_V03_S3_V06_S2_ID164_MuEG2011AAug05v4-eff9861cde6b1617ce28065ed0ad5d45/USER',
            
            
            
            
            
            
            
            
            

    '1120':'/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1120_ggToH120toWWto2L2Nuv2-27316856f5d7d57332f8edbc6c095e34/USER',
    '1130':'/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1130_ggToH130toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1140':'/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1140_ggToH140toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1150':'/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1150_ggToH150toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1160':'/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1160_ggToH160toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1170':'/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1170_ggToH170toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1180':'/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1180_ggToH180toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1190':'/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1190_ggToH190toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1200':'/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID1200_ggToH200toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '1210':'/GluGluToHToWWTo2L2Nu_M-210_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1210_ggToH210toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1220':'/GluGluToHToWWTo2L2Nu_M-220_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1220_ggToH220toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1230':'/GluGluToHToWWTo2L2Nu_M-230_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1230_ggToH230toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1250':'/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1250_ggToH250toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1300':'/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1300_ggToH300toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1350':'/GluGluToHToWWTo2L2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1350_ggToH350toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1400':'/GluGluToHToWWTo2L2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1400_ggToH400toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1450':'/GluGluToHToWWTo2L2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1450_ggToH450toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1500':'/GluGluToHToWWTo2L2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1500_ggToH500toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1550':'/GluGluToHToWWTo2L2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1550_ggToH550toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '1600':'/GluGluToHToWWTo2L2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID1600_ggToH600toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
            

    '2120':'/GluGluToHToWWToLNuTauNu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2120_ggToH120toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2130':'/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2130_ggToH130toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2140':'/GluGluToHToWWToLNuTauNu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2140_ggToH140toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2150':'/GluGluToHToWWToLNuTauNu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2150_ggToH150toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2160':'/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2160_ggToH160toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2170':'/GluGluToHToWWToLNuTauNu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2170_ggToH170toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2180':'/GluGluToHToWWToLNuTauNu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2180_ggToH180toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2190':'/GluGluToHToWWToLNuTauNu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2190_ggToH190toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2200':'/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2200_ggToH200toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2210':'/GluGluToHToWWToLNuTauNu_M-210_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2210_ggToH210toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2220':'/GluGluToHToWWToLNuTauNu_M-220_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2220_ggToH220toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2230':'/GluGluToHToWWToLNuTauNu_M-230_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2230_ggToH230toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2250':'/GluGluToHToWWToLNuTauNu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2250_ggToH250toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2300':'/GluGluToHToWWToLNuTauNu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2300_ggToH300toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2350':'/GluGluToHToWWToLNuTauNu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2350_ggToH350toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2400':'/GluGluToHToWWToLNuTauNu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2400_ggToH400toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2450':'/GluGluToHToWWToLNuTauNu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2450_ggToH450toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2500':'/GluGluToHToWWToLNuTauNu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2500_ggToH500toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2550':'/GluGluToHToWWToLNuTauNu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2550_ggToH550toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '2600':'/GluGluToHToWWToLNuTauNu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID2600_ggToH600toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
            

    '3120':'/GluGluToHToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3120_ggToH120toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3130':'/GluGluToHToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID3130_ggToH130toWWto2Tau2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '3140':'/GluGluToHToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID3140_ggToH140toWWto2Tau2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '3150':'/GluGluToHToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID3150_ggToH150toWWto2Tau2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '3160':'/GluGluToHToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID3160_ggToH160toWWto2Tau2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '3170':'/GluGluToHToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3170_ggToH170toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3180':'/GluGluToHToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3180_ggToH180toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3190':'/GluGluToHToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3190_ggToH190toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3200':'/GluGluToHToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3200_ggToH200toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3210':'/GluGluToHToWWTo2Tau2Nu_M-210_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3210_ggToH210toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3220':'/GluGluToHToWWTo2Tau2Nu_M-220_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3220_ggToH220toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3230':'/GluGluToHToWWTo2Tau2Nu_M-230_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3230_ggToH230toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3250':'/GluGluToHToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3250_ggToH250toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3300':'/GluGluToHToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3300_ggToH300toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3350':'/GluGluToHToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3350_ggToH350toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3400':'/GluGluToHToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3400_ggToH400toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3450':'/GluGluToHToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3450_ggToH450toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3500':'/GluGluToHToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3500_ggToH500toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3550':'/GluGluToHToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3550_ggToH550toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '3600':'/GluGluToHToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID3600_ggToH600toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
            
    '4120':'/VBF_HToWWTo2L2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4120_vbfToH120toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4130':'/VBF_HToWWTo2L2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4130_vbfToH130toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4140':'/VBF_HToWWTo2L2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4140_vbfToH140toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4150':'/VBF_HToWWTo2L2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4150_vbfToH150toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4160':'/VBF_HToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4160_vbfToH160toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4170':'/VBF_HToWWTo2L2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4170_vbfToH170toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4180':'/VBF_HToWWTo2L2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4180_vbfToH180toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4190':'/VBF_HToWWTo2L2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4190_vbfToH190toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4200':'/VBF_HToWWTo2L2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4200_vbfToH200toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4210':'/VBF_HToWWTo2L2Nu_M-210_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4210_vbfToH210toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4220':'/VBF_HToWWTo2L2Nu_M-220_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4220_vbfToH220toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4230':'/VBF_HToWWTo2L2Nu_M-230_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4230_vbfToH230toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4250':'/VBF_HToWWTo2L2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4250_vbfToH250toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4300':'/VBF_HToWWTo2L2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID4300_vbfToH300toWWto2L2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '4350':'/VBF_HToWWTo2L2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4350_vbfToH350toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4400':'/VBF_HToWWTo2L2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4400_vbfToH400toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4450':'/VBF_HToWWTo2L2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4450_vbfToH450toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4500':'/VBF_HToWWTo2L2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4500_vbfToH500toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4550':'/VBF_HToWWTo2L2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4550_vbfToH550toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '4600':'/VBF_HToWWTo2L2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID4600_vbfToH600toWWto2L2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
            
    '5120':'/VBF_HToWWToLNuTauNu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5120_vbfToH120toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5130':'/VBF_HToWWToLNuTauNu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5130_vbfToH130toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5140':'/VBF_HToWWToLNuTauNu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5140_vbfToH140toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5150':'/VBF_HToWWToLNuTauNu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5150_vbfToH150toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5160':'/VBF_HToWWToLNuTauNu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5160_vbfToH160toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5170':'/VBF_HToWWToLNuTauNu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5170_vbfToH170toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5180':'/VBF_HToWWToLNuTauNu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5180_vbfToH180toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5190':'/VBF_HToWWToLNuTauNu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5190_vbfToH190toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5200':'/VBF_HToWWToLNuTauNu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5200_vbfToH200toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5210':'/VBF_HToWWToLNuTauNu_M-210_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5210_vbfToH210toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5220':'/VBF_HToWWToLNuTauNu_M-220_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5220_vbfToH220toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5230':'/VBF_HToWWToLNuTauNu_M-230_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5230_vbfToH230toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5250':'/VBF_HToWWToLNuTauNu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5250_vbfToH250toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5300':'/VBF_HToWWToLNuTauNu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5300_vbfToH300toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5350':'/VBF_HToWWToLNuTauNu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5350_vbfToH350toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5400':'/VBF_HToWWToLNuTauNu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5400_vbfToH400toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5450':'/VBF_HToWWToLNuTauNu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5450_vbfToH450toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5500':'/VBF_HToWWToLNuTauNu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5500_vbfToH500toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5550':'/VBF_HToWWToLNuTauNu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5550_vbfToH550toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
    '5600':'/VBF_HToWWToLNuTauNu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID5600_vbfToH600toWWtoLNuTauNu-a0fee9021ee0e974e473d874349bbb62/USER',
            
    '6120':'/VBF_HToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6120_vbfToH120toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6130':'/VBF_HToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6130_vbfToH130toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6140':'/VBF_HToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6140_vbfToH140toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6150':'/VBF_HToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6150_vbfToH150toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6160':'/VBF_HToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6160_vbfToH160toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6170':'/VBF_HToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6170_vbfToH170toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6180':'/VBF_HToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6180_vbfToH180toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6190':'/VBF_HToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6190_vbfToH190toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6200':'/VBF_HToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6200_vbfToH200toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6210':'/VBF_HToWWTo2Tau2Nu_M-210_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID6210_vbfToH210toWWto2Tau2Nu-27316856f5d7d57332f8edbc6c095e34/USER',
    '6220':'/VBF_HToWWTo2Tau2Nu_M-220_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6220_vbfToH220toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6230':'/VBF_HToWWTo2Tau2Nu_M-230_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6230_vbfToH230toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6250':'/VBF_HToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6250_vbfToH250toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6300':'/VBF_HToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6300_vbfToH300toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6350':'/VBF_HToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6350_vbfToH350toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6400':'/VBF_HToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6400_vbfToH400toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6450':'/VBF_HToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6450_vbfToH450toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6500':'/VBF_HToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6500_vbfToH500toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6550':'/VBF_HToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6550_vbfToH550toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',
    '6600':'/VBF_HToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6/mwlebour-R42X_S1_V06_S2_V02_S3_V05_S2_ID6600_vbfToH600toWWto2Tau2Nu-a0fee9021ee0e974e473d874349bbb62/USER',

    '7120':'/WH_ZH_TTH_HToWW_M-120_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7120_wzttH120ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7130':'/WH_ZH_TTH_HToWW_M-130_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7130_wzttH130ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7140':'/WH_ZH_TTH_HToWW_M-140_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7140_wzttH140ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7150':'/WH_ZH_TTH_HToWW_M-150_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7150_wzttH150ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7160':'/WH_ZH_TTH_HToWW_M-160_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7160_wzttH160ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7170':'/WH_ZH_TTH_HToWW_M-170_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7170_wzttH170ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7180':'/WH_ZH_TTH_HToWW_M-180_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7180_wzttH180ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7190':'/WH_ZH_TTH_HToWW_M-190_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7190_wzttH190ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7200':'/WH_ZH_TTH_HToWW_M-200_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7200_wzttH200ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7250':'/WH_ZH_TTH_HToWW_M-250_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7250_wzttH250ToWW-27316856f5d7d57332f8edbc6c095e34/USER',
    '7300':'/WH_ZH_TTH_HToWW_M-300_7TeV-pythia6/mwlebour-R42X_S1_V06_S2_V04_S3_V06_S2_ID7300_wzttH300ToWW-27316856f5d7d57332f8edbc6c095e34/USER',

}

