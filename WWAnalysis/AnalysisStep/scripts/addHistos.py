from ROOT import gROOT, TCanvas, TH1F, TFile
from ROOT import kCyan, kGreen, kMagenta, kBlue

# gROOT.Reset()
c1 = TCanvas('c1','Example with Formula')
 
f = TFile('all.root')
keys = [key.GetName() for key in f.GetListOfKeys()]

channels = ['all']
# channels = ['all','wwelel','wwelmu','wwmumu']

# nbins = f.Get(keys[0]).GetNbinsX()
nbins = 16

# subPlots = {
#     'w' :[kCyan-9,   'W+Jets',    ['we','wm','wt'],        TH1F('w','w',nbins,0,nbins)],
#     'vv':[kGreen-9,  'diBoson',   ['ww','wz','zz'],        TH1F('ww','ww',nbins,0,nbins)],
#     't' :[kMagenta-9,'t#bar{t},t',['ttbar','tt','ts','tw'],TH1F('t','t',nbins,0,nbins)],
#     'dy':[kBlue-9,   'Drell-Yan', ['dyee','dymm','dytt'],  TH1F('dy','dy',nbins,0,nbins)],
# }

backgrounds = ['we','wm','wt','wz','zz','ttbar','tt','ts','tw','dyee','dymm','dytt']
signals     = ['ww']

counts = {}
for k in keys:
#     for channel in channels:
    if k.split('_')[0] in channels:
        counts[k] = {}
        h = f.Get(k)
#        for plot in subPlots:
#            if k in [channel+'_%s' %x for x in subPlots[plot][2]]:
#                subPlots[plot][3].Add(h) 
        for i in range(1,nbins+1):
            counts[k][i] = h.GetBinContent(i)

            
# c1.Divide(2,2)
# for plot in subPlots:
#     for channel in channels:
        

lines = {}
lines[0] = "  Cut  "
for i in range(1,nbins+1):
    lines[i] = "{0:^7}".format(i)

for sample in counts:
    lines[0] += "{0:^12}".format(sample)
    for cut in counts[sample]:
        lines[cut] += "{0:^12.3}".format(counts[sample][cut])

for line in lines:
    print lines[line]
