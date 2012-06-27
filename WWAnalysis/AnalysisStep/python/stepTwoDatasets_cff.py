import FWCore.ParameterSet.Config as cms
# to remove lines with no dataset:
# %s/^[0-9]\+ *$//
# convert rest of the lines to dictionaries
# :%s/^\([0-9]\+\) *\(.*\)$/    '\1':'\2',/

stepTwoDatasets = {
    '001':'/GluGluToWWTo4L_TuneZ2star_8TeV-gg2ww-pythia6/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID001_GluGluToWWTo4L-bc58659808d390d23203e20a37543a0e/USER',
	'000':'/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID000_WWJets2LMad-bc58659808d390d23203e20a37543a0e/USER',

	'100':'/SingleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID100_SingleElectron2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'101':'/SingleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID101_SingleMuon2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'102':'/DoubleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID102_DoubleElectron2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'103':'/DoubleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID103_DoubleMuon2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'104':'/MuEG/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID104_MuEG2012A-5d3127590cd54452f5f3b5b5f9c49e55/USER',

	'200':'/SingleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID200_SingleElectron2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'201':'/SingleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID201_SingleMuon2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'202':'/DoubleElectron/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID202_DoubleElectron2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'203':'/DoubleMu/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID203_DoubleMuon2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',
	'204':'/MuEG/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID204_MuEG2012B-5d3127590cd54452f5f3b5b5f9c49e55/USER',

	'012':'/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID012_TbartWFullDR-bc58659808d390d23203e20a37543a0e/USER',
	'011':'/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID011_TtWFullDR-bc58659808d390d23203e20a37543a0e/USER',
	'010c1':'/TTJets_TuneZ2star_8TeV-madgraph-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID010_TTJetsMadChunck1-bc58659808d390d23203e20a37543a0e/USER',
        '010c2':'/TTJets_TuneZ2star_8TeV-madgraph-tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID010_TTJetsMadChunck2-bc58659808d390d23203e20a37543a0e/USER',

	'080':'/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID080_WJetsToLNuMad-bc58659808d390d23203e20a37543a0e/USER',

	'070':'/WZ_TuneZ2star_8TeV_pythia6_tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID070_WZFull-bc58659808d390d23203e20a37543a0e/USER',
	'077':'/ZZTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID077_ZZTo2L2Nu-bc58659808d390d23203e20a37543a0e/USER',

	'036':'/DYJetsToLL_M-10To50filter_8TeV-madgraph/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID036_DY10toLLMad-bc58659808d390d23203e20a37543a0e/USER',

	'037c0':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck0-bc58659808d390d23203e20a37543a0e/USER',
        '037c1':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck1-bc58659808d390d23203e20a37543a0e/USER',
        '037c2':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck2-bc58659808d390d23203e20a37543a0e/USER',
        '037c3':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck3-bc58659808d390d23203e20a37543a0e/USER',
        '037c4':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck4-bc58659808d390d23203e20a37543a0e/USER',
        '037c5':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck5-bc58659808d390d23203e20a37543a0e/USER',
        '037c6':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck6-bc58659808d390d23203e20a37543a0e/USER',
        '037c7':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck7-bc58659808d390d23203e20a37543a0e/USER',
        '037c8':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck8-bc58659808d390d23203e20a37543a0e/USER',
        '037c9':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck9-bc58659808d390d23203e20a37543a0e/USER',


        '3110':'/WH_ZH_TTH_HToWW_M-110_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3110_wzttH110ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3115':'/WH_ZH_TTH_HToWW_M-115_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3115_wzttH115ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3118':'wzttH118ToWW',
        '3120':'/WH_ZH_TTH_HToWW_M-120_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3120_wzttH120ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3122':'wzttH122ToWW',
        '3124':'wzttH124ToWW',
        '3125':'/WH_ZH_TTH_HToWW_M-125_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3125_wzttH125ToWW_31May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3126':'wzttH126ToWW',
        '3128':'wzttH128ToWW',
        '3130':'/WH_ZH_TTH_HToWW_M-130_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3130_wzttH130ToWW_02Jun2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3135':'/WH_ZH_TTH_HToWW_M-135_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3135_wzttH135ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3140':'/WH_ZH_TTH_HToWW_M-140_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3140_wzttH140ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3145':'/WH_ZH_TTH_HToWW_M-145_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3145_wzttH145ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3150':'/WH_ZH_TTH_HToWW_M-150_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3150_wzttH150ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3155':'/WH_ZH_TTH_HToWW_M-155_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3155_wzttH155ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3160':'/WH_ZH_TTH_HToWW_M-160_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3160_wzttH160ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3170':'/WH_ZH_TTH_HToWW_M-170_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3170_wzttH170ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3180':'/WH_ZH_TTH_HToWW_M-180_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3180_wzttH180ToWW_31May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3190':'/WH_ZH_TTH_HToWW_M-190_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3190_wzttH190ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3200':'/WH_ZH_TTH_HToWW_M-200_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3200_wzttH200ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3250':'/WH_ZH_TTH_HToWW_M-250_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3250_wzttH250ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3300':'/WH_ZH_TTH_HToWW_M-300_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3300_wzttH300ToWW_31May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3350':'wzttH350ToWW',
        '3400':'wzttH400ToWW',
        '3450':'wzttH450ToWW',
        '3500':'wzttH500ToWW',
        '3550':'wzttH550ToWW',
        '3600':'wzttH600ToWW',   



}
