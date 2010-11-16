# fix global tag
process.GlobalTag.globaltag = 'GR_R_38X_V14::All'

# JSON
import re
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
minRun = 0
for f in "Cert_132440-144114_7TeV_Sep17ReReco_Collisions10_CMSSWConfig.txt", "Cert_132440-149442_7TeV_StreamExpress_Collisions10_CMSSWConfig_v2.txt":
    jsontxt = open(f,"r")
    for l in jsontxt:
        if "process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(": break
    for l in jsontxt:
        m = re.search(r"'((\d+):\d+-\d+:\d+)',",l)
        if m:
            run = int(m.group(2))
            if run < minRun: continue
            minRun = run
            process.source.lumisToProcess.append(m.group(1))


# leave some lines blank at the end please


