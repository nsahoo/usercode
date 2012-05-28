import FWCore.ParameterSet.Config as cms
# to remove lines with no dataset:
# %s/^[0-9]\+ *$//
# convert rest of the lines to dictionaries
# :%s/^\([0-9]\+\) *\(.*\)$/    '\1':'\2',/

stepTwoDatasets = {
    '001':'/GluGluToWWTo4L_TuneZ2star_8TeV-gg2ww-pythia6/jfernan2-R52X_S1_V04_S2_V02_S3_V03_ID001_GluGluToWWTo4L-bc58659808d390d23203e20a37543a0e/USER',

}
