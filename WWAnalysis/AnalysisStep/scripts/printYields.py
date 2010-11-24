from ROOT import gROOT, TCanvas, TH1F, TFile
from ROOT import kCyan, kGreen, kMagenta, kBlue

f = TFile('all.root')
keys = [key.GetName() for key in f.GetListOfKeys()]

channels = ['all']
# channels = ['all','wwelel','wwelmu','wwmumu']

nbins = 16

backgrounds = ['we','wm','wt','wz','zz','ttbar','tt','ts','tw','dyee','dymm','dytt']
signals     = ['ww']

counts = {}
for k in keys:
    if k.split('_')[0] in channels:
        counts[k] = {}
        h = f.Get(k)
        for i in range(1,nbins+1):
            counts[k][i] = h.GetBinContent(i)

            
lines = {}
lines[0] = "  Cut  "
for i in range(1,nbins+1):
    lines[i] = "{0:^7}".format(i)

all = backgrounds
all.extend(signals)
for sample in all:
    lines[0] += "{0:^12}".format(sample)
    for cut in counts['all_'+sample]:
        lines[cut] += "{0:^12.3}".format(counts['all_'+sample][cut])

for line in lines:
    print lines[line]
