import FWCore.ParameterSet.Config as cms
#  regex: ^\([0-9]\+\) *\(.*\)$
# full line :%s/^\([0-9]\+\) *\(.*\)$/^I'\1':'\2',/


stepOneDatasets = {

	'001':'/GluGluToWWTo4L_TuneZ2star_8TeV-gg2ww-pythia6/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID001_GluGluToWWTo4L-bc58659808d390d23203e20a37543a0e/USER',

        '2110':'/VBF_HToWWTo2LAndTau2Nu_M-110_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2110_vbfToH110toWWTo2LAndTau2Nu_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2115':'vbfToH115toWWTo2LAndTau2Nu',
        '2118':'vbfToH118toWWTo2LAndTau2Nu',
        '2120':'vbfToH120toWWTo2LAndTau2Nu',
        '2122':'vbfToH122toWWTo2LAndTau2Nu',
        '2124':'vbfToH124toWWTo2LAndTau2Nu',
        '2125':'vbfToH125toWWTo2LAndTau2Nu',
        '2126':'vbfToH126toWWTo2LAndTau2Nu',
        '2128':'vbfToH128toWWTo2LAndTau2Nu',
        '2130':'/VBF_HToWWTo2LAndTau2Nu_M-130_8TeV-powheg-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID2130_vbfToH130toWWTo2LAndTau2Nu_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '2135':'vbfToH135toWWTo2LAndTau2Nu',
        '2140':'vbfToH140toWWTo2LAndTau2Nu',
        '2145':'vbfToH145toWWTo2LAndTau2Nu',
        '2150':'vbfToH150toWWTo2LAndTau2Nu',
        '2155':'vbfToH155toWWTo2LAndTau2Nu',
        '2160':'vbfToH160toWWTo2LAndTau2Nu',
        '2170':'vbfToH170toWWTo2LAndTau2Nu',
        '2180':'vbfToH180toWWTo2LAndTau2Nu',
        '2190':'vbfToH190toWWTo2LAndTau2Nu',
        '2200':'vbfToH200toWWTo2LAndTau2Nu',
        '2300':'vbfToH300toWWTo2LAndTau2Nu',
        '2400':'vbfToH400toWWTo2LAndTau2Nu',
        '2500':'vbfToH500toWWTo2LAndTau2Nu',
        '2600':'vbfToH600toWWTo2LAndTau2Nu',
        '2700':'vbfToH700toWWTo2LAndTau2Nu',
        '2800':'vbfToH800toWWTo2LAndTau2Nu',
        '2900':'vbfToH900toWWTo2LAndTau2Nu',
        '3000':'vbfToH1000toWWTo2LAndTau2Nu',



        '3110':'wzttH110ToWW',
        '3115':'/WH_ZH_TTH_HToWW_M-115_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3115_wzttH115ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3118':'wzttH118ToWW',
        '3120':'/WH_ZH_TTH_HToWW_M-120_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3120_wzttH120ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3122':'wzttH122ToWW',
        '3124':'wzttH124ToWW',
        '3125':'wzttH125ToWW',
        '3126':'wzttH126ToWW',
        '3128':'wzttH128ToWW',
        '3130':'wzttH130ToWW',
        '3135':'/WH_ZH_TTH_HToWW_M-135_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3135_wzttH135ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3140':'/WH_ZH_TTH_HToWW_M-140_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3140_wzttH140ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3145':'wzttH145ToWW',
        '3150':'/WH_ZH_TTH_HToWW_M-150_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3150_wzttH150ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3155':'/WH_ZH_TTH_HToWW_M-155_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3155_wzttH155ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3160':'wzttH160ToWW',
        '3170':'/WH_ZH_TTH_HToWW_M-170_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3170_wzttH170ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3180':'wzttH180ToWW',
        '3190':'/WH_ZH_TTH_HToWW_M-190_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3190_wzttH190ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3200':'wzttH200ToWW',
        '3250':'/WH_ZH_TTH_HToWW_M-250_8TeV-pythia6/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID3250_wzttH250ToWW_28May2012-e0faaf2eb7156623b7cf3c14c660b46e/USER',
        '3300':'wzttH300ToWW',
        '3350':'wzttH350ToWW',
        '3400':'wzttH400ToWW',
        '3450':'wzttH450ToWW',
        '3500':'wzttH500ToWW',
        '3550':'wzttH550ToWW',
        '3600':'wzttH600ToWW',   





        '103TT':'/DoubleMu/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID103TT_DoubleMu2012ADYttem_28May2012-cee84bd278931e99af50d92faa75b2af/USER',
        '113TT':'/DoubleMu/amassiro-R52X_S1_V04_S2_V02_S3_V03_ID113TT_DoubleMu2012BDYttem_28May2012-cee84bd278931e99af50d92faa75b2af/USER',




# dy

        '037a':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck0-bc58659808d390d23203e20a37543a0e/USER',    
        '037b':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck1-bc58659808d390d23203e20a37543a0e/USER',    
        '037c':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck2-bc58659808d390d23203e20a37543a0e/USER',    
        '037d':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck3-bc58659808d390d23203e20a37543a0e/USER',    
        '037e':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck4-bc58659808d390d23203e20a37543a0e/USER',    
        '037f':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck5-bc58659808d390d23203e20a37543a0e/USER',    
        '037g':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck6-bc58659808d390d23203e20a37543a0e/USER',    
        '037h':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck7-bc58659808d390d23203e20a37543a0e/USER',    
        '037i':'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID037_DY50toLLMadChunck8-bc58659808d390d23203e20a37543a0e/USER',    


}


