from WWAnalysis.AnalysisStep.mcAnalysis import *
from optparse import OptionParser
parser = OptionParser(usage="%prog [options] mcanalysis.txt cuts.txt")
addMCAnalysisOptions(parser)
(options, args) = parser.parse_args()
options.nMinusOne = False
options.weight    = True
mca = MCAnalysis(args[0], options)
mca.scaleProcess("Top",  1.0)
cf0j = CutsFile(args[1],options)
cf1j = CutsFile(cf0j).replace('jet veto', 'one jet', 'njet == 1 && dphilljet*sameflav < 165./180.*3.1415926').replace('top veto', 'top veto', 'bveto && nbjet == 0')
cf2j = CutsFile(cf0j).replace('jet veto', 'one jet', 'njet == 2').replace('top veto', 'top veto', 'bveto && nbjet == 0')
selections = { 
    '0j0b' : cf0j,
    '1j0b' : cf1j,
    #'2j'   : cf2j,
    '0j0bip'  : CutsFile(cf0j).replace('top veto','soft ip veto','bveto_ip'),
    '0j1sbip' : CutsFile(cf0j).replace('top veto','soft b ip',   '!bveto_ip'),
    #'0j1sbmu' : CutsFile(cf0j).replace('top veto','soft b mu','!bveto_mu && bveto_ip'),
    '0j1sb'   : CutsFile(cf0j).replace('top veto','soft b','!bveto'),
    '1j0hb'   : CutsFile(cf1j).replace('top veto','hard b veto','nbjet == 0'),
    '1j0hb1sb' : CutsFile(cf1j).replace('top veto','hard b veto','nbjet == 0 && !bveto_ip'),
    '1j0hb0ip' : CutsFile(cf1j).replace('top veto','hard b veto','nbjet == 0 && !bveto_ip'),
    '1j1sbip'  : CutsFile(cf1j).replace('top veto','soft b ip','!bveto_ip'),
    #'1j1sbmu' : CutsFile(cf1j).replace('top veto','soft b mu','!bveto_mu && bveto_ip'),
    #'1j1sb'   : CutsFile(cf1j).replace('top veto','soft b','!bveto'),
    '1j1hb'   : CutsFile(cf1j).replace('top veto','hard b','nbjet == 1'),
    '1j1hb1ip': CutsFile(cf1j).replace('top veto','hard b','nbjet == 1 && !bveto_ip'),
    '1j1hb0sb': CutsFile(cf1j).replace('top veto','hard b but no soft b','bveto && nbjet == 1'),
    '2jinc'   : CutsFile(cf2j).remove('top veto'),
    '2j1hb'   : CutsFile(cf2j).replace('top veto','hard b',    'nbjet == 1'),
    '2j2hb'   : CutsFile(cf2j).replace('top veto','two hard b','nbjet == 2'),
    '2j1hb0sip' : CutsFile(cf2j).replace('top veto','hard b',    'nbjet == 1 && bveto_ip'),
    '2j2hb0sip' : CutsFile(cf2j).replace('top veto','two hard b','nbjet == 2 && bveto_ip'),
    '2j1sbmu' : CutsFile(cf1j).replace('top veto','soft b mu', '!bveto_mu'),
}

#print "\n ==== All Selections (full) ==== "
reports = {}; yields = {}
index = -1 if options.inclusive else EM
for S,C in selections.items():
    reports[S] = mca.getYields(C)
    yields[S] = {
        'top':   reports[S]['Top' ][-1][1][index][1],
        'data':  reports[S]['data'][-1][1][index][1],
        'other': sum( reports[S][p][-1][1][index][1] for p in reports[S].keys() if p not in ['data','Top']),
    }
    yields[S]['all'] = yields[S]['top'] + yields[S]['other']
    #print " ---- %s ---- " % S
    #mca.prettyPrint(reports[S])

def sorted(x):
    all = [i for i in x]
    all.sort()
    return all

print "\n ==== All Selections (overview) ==== "
for S, y in sorted(yields.items()):
    print "Selection %-10s:  n(top, sim.) = %8.2f    n(other, sim.) = %8.2f       n(data) = %8.2f" % ( S, y['top'], y['other'], y['data'] )

print ""

def effND(num,den): 
    if den == 0: return 0
    eff = num/float(den)
    return ( eff, sqrt(eff*(1-eff)/den) )
def eff12(n1,n2): 
    if n1+n2 == 0: return 0
    eff = (2.*n2)/(n1+2.*n2);
    return ( eff, 2*sqrt(n1*n2*(n1+n2))/pow(n1+2.*n2,2) )
def eff2_from_eff1(eff1): return ( eff1[0]*(2-eff1[0]), (1-eff1[0])*eff1[1] )
def eff_over_ineff(eff):  return ( eff[0]/(1.0 - eff[0]), eff[1]/pow(1.0-eff[0], 2) )
def ineff_over_eff(eff):  return ( (1.0 - eff[0])/eff[0], eff[1]/pow(eff[0], 2) )
def N_times_alpha(N,alpha): return ( N*alpha[0], sqrt(N)*alpha[0], N*alpha[1] )
def scalefact(eff,mc): return (eff[0]/mc[0], eff[1]/mc[0])
def scalevec(eff,scale): return (eff[0]*scale, eff[1]*scale)
def mean(*effs): 
    """Do the right thing"""
    ws = [ 1./pow(eff[1], 2.0) for eff in effs ]
    avgeff = sum([eff[0]*w for eff,w in zip(effs,ws)])/sum(ws)
    avgerr = sqrt(1.0/sum(ws))
    return ( avgeff, avgerr )
### Now estimating tagging efficiencies
print "\n ==== All Efficiencies ==== "

# soft ip tag from 1 jet bin with hard jet tagged (no background subtraction)
effsoftip_data = effND(yields['1j1hb1ip']['data'], yields['1j1hb']['data']);
effsoftip_simt = effND(yields['1j1hb1ip']['top'] , yields['1j1hb']['top'] )
effsoftip_sima = effND(yields['1j1hb1ip']['all'] , yields['1j1hb']['all'] )
print "Eff(soft ip) from 1 jet bin:  data %.1f +/- %.1f %%   sim(top) %.1f%%  (s.f. %4.2f)    sim(all) %.1f%% (s.f. %4.2f)" % (
         effsoftip_data[0]*100, effsoftip_data[1]*100, 
         effsoftip_simt[0]*100, effsoftip_data[0]/effsoftip_simt[0], 
         effsoftip_sima[0]*100, effsoftip_data[0]/effsoftip_sima[0]) 

# soft mu tag from 2 jet bin with no tagging requirement (no background subtraction)
effsoftmu2_data = effND(yields['2j1sbmu']['data'], yields['1j1hb']['data']);
effsoftmu2_simt = effND(yields['2j1sbmu']['top'] , yields['1j1hb']['top'] )
effsoftmu2_sima = effND(yields['2j1sbmu']['all'] , yields['1j1hb']['all'] )
print "Eff(soft mu) from 2 jet bin:  data %.1f +/- %.1f %%   sim(top) %.1f%%  (s.f. %4.2f)    sim(all) %.1f%% (s.f. %4.2f)" % (
         effsoftmu2_data[0]*100, effsoftmu2_data[1]*100, 
         effsoftmu2_simt[0]*100, effsoftmu2_data[0]/effsoftmu2_simt[0], 
         effsoftmu2_sima[0]*100, effsoftmu2_data[0]/effsoftmu2_sima[0]) 

# hard ip tag from 1 jet bin soft jet tagged (no background subtraction)
effhard1ip_data = effND(yields['1j1hb1ip']['data'], yields['1j1sbip']['data']);
effhard1ip_simt = effND(yields['1j1hb1ip']['top'] , yields['1j1sbip']['top'] )
effhard1ip_sima = effND(yields['1j1hb1ip']['all'] , yields['1j1sbip']['all'] )
print "Eff(hard) from 1 jet bin with soft ip tag:   data %.1f +/- %.1f %%   sim(top) %.1f%%  (s.f. %4.2f)    sim(all) %.1f%% (s.f. %4.2f)" % (
         effhard1ip_data[0]*100, effhard1ip_data[1]*100, 
         effhard1ip_simt[0]*100, effhard1ip_data[0]/effhard1ip_simt[0], 
         effhard1ip_sima[0]*100, effhard1ip_data[0]/effhard1ip_sima[0]) 

# hard ip tag from 2 jet bin soft jet tagged (no background subtraction)
# N(1tag) = 2*eff*(1-eff);  N(2 tag) = eff*eff
effhard2j1b_data = eff12(yields['2j1hb']['data'], yields['2j2hb']['data']);
effhard2j1b_simt = eff12(yields['2j1hb']['top'] , yields['2j2hb']['top'] )
effhard2j1b_sima = eff12(yields['2j1hb']['all'] , yields['2j2hb']['all'] )
print "Eff(hard) from 2 jet bin (1 tag vs 2 tags):  data %.1f +/- %.1f %%   sim(top) %.1f%%  (s.f. %4.2f)    sim(all) %.1f%% (s.f. %4.2f)" % (
         effhard2j1b_data[0]*100, effhard2j1b_data[1]*100, 
         effhard2j1b_simt[0]*100, effhard2j1b_data[0]/effhard2j1b_simt[0], 
         effhard2j1b_sima[0]*100, effhard2j1b_data[0]/effhard2j1b_sima[0]) 

# hard ip tag from 2 jet bin soft jet tagged (no background subtraction)
# N(1tag) = 2*eff*(1-eff);  N(2 tag) = eff*eff
effhard2j1b_data = eff12(yields['2j1hb0sip']['data'], yields['2j2hb0sip']['data']);
effhard2j1b_simt = eff12(yields['2j1hb0sip']['top'] , yields['2j2hb0sip']['top'] )
effhard2j1b_sima = eff12(yields['2j1hb0sip']['all'] , yields['2j2hb0sip']['all'] )
print "      applying soft ip veto at denominator:  data %.1f +/- %.1f %%   sim(top) %.1f%%  (s.f. %4.2f)    sim(all) %.1f%% (s.f. %4.2f)" % (
         effhard2j1b_data[0]*100, effhard2j1b_data[1]*100, 
         effhard2j1b_simt[0]*100, effhard2j1b_data[0]/effhard2j1b_simt[0], 
         effhard2j1b_sima[0]*100, effhard2j1b_data[0]/effhard2j1b_sima[0]) 

sfhardr2j1b = scalefact(ineff_over_eff(effhard2j1b_data), ineff_over_eff(effhard2j1b_sima))
sfhardr1ip  = scalefact(ineff_over_eff(effhard1ip_data),  ineff_over_eff(effhard1ip_sima))
sfhardr = mean(sfhardr2j1b, sfhardr1ip)
print "Scale factor for [1-eff(hard)]/eff(hard): %4.2f +/- %4.2f from soft ip, %4.2f +/- %4.2f from 2j, %4.2f +/- %4.2f combined" % (
            sfhardr1ip[0], sfhardr1ip[1], sfhardr2j1b[0], sfhardr2j1b[1], sfhardr[0], sfhardr[1] )

### Now estimating yields 
print "\n ==== All Background Yields  ==== "

# top scale factor for 0 jet bin 
eff2bip_data = eff2_from_eff1(effsoftip_data); alpha0jip_data = ineff_over_eff(eff2bip_data)
eff2bip_simt = eff2_from_eff1(effsoftip_simt); alpha0jip_simt = ineff_over_eff(eff2bip_simt)
eff2bip_sima = eff2_from_eff1(effsoftip_sima); alpha0jip_sima = ineff_over_eff(eff2bip_sima)
top0jip_data_nosub = N_times_alpha(yields['0j1sbip']['data'], alpha0jip_data)
top0jip_sima_nosub = N_times_alpha(yields['0j1sbip']['all' ], alpha0jip_sima)
top0jip_data_sub   = N_times_alpha(yields['0j1sbip']['data'] - yields['0j1sbip']['other'], alpha0jip_data)
top0jip_simt_sub   = N_times_alpha(yields['0j1sbip']['all' ] - yields['0j1sbip']['other'], alpha0jip_simt)
print "\nNow using full soft-b tagging (ip + mu)"
print "N(top 0j soft b, Data): %5.2f +/- %5.2f (stat) with no background subtraction        "  % (yields['0j1sb']['data'], sqrt(yields['0j1sb']['data']))
print "N(top 0j soft b, Data): %5.2f +/- %5.2f (stat) +/- %5.2f (sub) with sim. subtraction "  % (yields['0j1sb']['data'] - yields['0j1sb']['other'], sqrt(yields['0j1sb']['data']), yields['0j1sb']['other'])
scaleFac0j_N =  ( (yields['0j1sb']['data'] - yields['0j1sb']['other'])/(yields['0j1sb']['top']),
                  hypot(sqrt(yields['0j1sb']['data']), yields['0j1sb']['other'])/(yields['0j1sb']['top']) )

print "Scale factor up to here: %4.2f +/- %4.2f" % scaleFac0j_N

print ""
eff0jip_true = yields['0j1sbip']['top']/(yields['0j1sb']['top']+yields['0j0b']['top'])
effRatio0j1j = eff0jip_true/eff2bip_simt[0]
print "Eff(soft ip tag, 0j, Data): %5.3f +/- %5.3f (no background sub.)" % eff2bip_data
print "Eff(soft ip tag, 0j, Sim.): %5.3f           (no background sub.)" % eff2bip_sima[0]
print "Eff(soft ip tag, 0j, True): %5.3f           (from 1 jet bin)    " % eff2bip_simt[0]
print "Eff(soft ip tag, 0j, True): %5.3f           (from 0 jet bin)    " % (eff0jip_true) 
print "             [ ratio 0j/1j: %5.3f ]  " % (effRatio0j1j)
eff2bip_data_scaled = scalevec(eff2bip_data, effRatio0j1j)
print "Eff(soft ip tag, 0j, Data): %5.3f +/- %5.3f (scaled to 0 jet bin)" % eff2bip_data_scaled
print ""
eff2muResidual_sim  = ((yields['0j1sb']['top']-yields['0j1sbip']['top'])/(yields['0j1sb']['top']+yields['0j0b']['top']), 0)
eff2muResidual_data = scalevec( scalefact(effsoftmu2_data, effsoftmu2_sima), eff2muResidual_sim[0] )
print "Eff(soft mu tag, 0j, Sim.): %5.3f on top of the soft ip tag" % eff2muResidual_sim[0]
print "Eff(soft mu tag, 0j, Data): %5.3f +/- %5.3f from taking the scale factor in the 2 jet measurement"  % eff2muResidual_data
print ""

effSoftTagTotal_sim =  (yields['0j1sb']['top']/(yields['0j1sb']['top']+yields['0j0b']['top']), 0)
effSoftTagTotal_data = (eff2bip_data_scaled[0] + eff2muResidual_data[0], hypot(eff2bip_data_scaled[1], eff2muResidual_data[1]))
print "Eff(any soft tag, 0j, Sim.):  %5.3f "          % effSoftTagTotal_sim[0]
print "Eff(any soft tag, 0j, Data):  %5.3f +-/ %5.3f" %  effSoftTagTotal_data

alphaTotal_sim  = ineff_over_eff(effSoftTagTotal_sim)
alphaTotal_data = ineff_over_eff(effSoftTagTotal_data)
scaleFac0j_alpha = scalefact(alphaTotal_data, alphaTotal_sim)
print "alpha(top any tag, 0j, Data): %5.3f +/- %5.3f " % alphaTotal_data
print "alpha(top any tag, 0j, Sim.): %5.3f +/- %5.3f " % alphaTotal_sim
print "Scale factor for alpha:  %4.2f +/- %4.2f" % scaleFac0j_alpha

print ""
scaleFac0j_total = ( scaleFac0j_alpha[0]*scaleFac0j_N[0], 
                     hypot(scaleFac0j_alpha[1]*scaleFac0j_N[0], scaleFac0j_alpha[0]*scaleFac0j_N[1]) )
print "Final scale factor in 0 jets:  %4.2f +/- %4.2f    [ dB/B = %4.2f ]" % (scaleFac0j_total[0], scaleFac0j_total[1], scaleFac0j_total[1]/scaleFac0j_total[0])

## Other method, not mainstream
# print "\nNow using only soft ip b tagging"
# print "alpha(top, 0j, Data): %5.3f +/- %5.3f " % alpha0jip_data
# print "alpha(top, 0j, Sim.): %5.3f +/- %5.3f " % alpha0jip_sima
# print "N(top, 0j, Data) from soft ip: %8.2f +/- %6.2f (N) +/- %6.2f (alpha) with no subtraction   " % top0jip_data_nosub
# print "N(top, 0j, Data) from soft ip: %8.2f +/- %6.2f (N) +/- %6.2f (alpha) with sim. subtraction " % top0jip_data_sub
# print "N(top, 0j, Sim.) from soft ip: %8.2f without any background subtraction in control region  " % top0jip_sima_nosub[0]
# print "N(top, 0j, Sim.) from soft ip: %8.2f for perfect background subtraction in control region  " % top0jip_simt_sub[0]
# print "N(top, 0j, Sim.) truth value : %8.2f                                                       " % yields['0j0bip']['top'] 
# print "Scale factor, 0j: %5.3f " % (top0jip_data_sub[0]/top0jip_simt_sub[0])
print ""

alpha1jMC = yields['1j0b']['top']/yields['1j1hb0sb']['top']
topYield1hb0sb    = yields['1j1hb0sb']['data'] - yields['1j1hb0sb']['other']
topYield1hb0sbErr = sqrt(yields['1j1hb0sb']['data'] + pow(yields['1j1hb0sb']['other'],2))
#topYield0hb2E = (topYield1hb0sb*alpha1jMC*sfhardr[0], topYield1hb0sbErr*alpha1jMC*sfhardr[0], topYield1hb0sb*alpha1jMC*sfhardr[1])
topYield0hb3E = (topYield1hb0sb*alpha1jMC*sfhardr[0], 
                 sqrt(yields['1j1hb0sb']['data'])*alpha1jMC*sfhardr[0], 
                 yields['1j1hb0sb']['other']*alpha1jMC*sfhardr[0], 
                 topYield1hb0sb*alpha1jMC*sfhardr[1])
topYield0hb1E = (topYield0hb3E[0], hypot(topYield0hb3E[1], hypot(topYield0hb3E[2], topYield0hb3E[3])))
print "[1-eff(hard)]/eff(hard) in 1 jet bin with soft b veto in sim.: %4.2f  " % (alpha1jMC)
print "        alpha (data) = alpha(mc) times scale factor from data: %4.2f +/- %4.2f " % (alpha1jMC*sfhardr[0], alpha1jMC*sfhardr[1])
print "Events in 1 jet bin with 1 hard b and soft b veto in data: %5.1f  " % yields['1j1hb0sb']['data']
print "                            expected background from sim.: %5.1f  " % yields['1j1hb0sb']['other']
print "              data after background subtraction from sim.: %5.1f +/- %4.1f " % ( topYield1hb0sb, topYield1hb0sbErr )
print "Estimated background in data: %5.1f +/- %4.1f (N) +/- %4.1f (bg.) +/- %4.1f (alpha) " % topYield0hb3E
print "                              %5.1f +/- %4.1f (all)     [ dB/B = %4.2f ] " % (topYield0hb1E[0], topYield0hb1E[1], topYield0hb1E[1]/topYield0hb1E[0])
print "Background from simulation:   %5.1f " % yields['1j0b']['top']
print "    data/sim. scale factor:   %4.2f +/- %4.2f   [ dB/B = %4.2f ] " % (topYield0hb1E[0]/yields['1j0b']['top'], topYield0hb1E[1]/yields['1j0b']['top'], topYield0hb1E[1]/topYield0hb1E[0])
#print "N(top, 1j, Sim.) from hard b: %8.2f without any background subtraction in control region  " % top1j_sima_nosub[0]
#print "N(top, 1j, Sim.) from hard b: %8.2f for perfect background subtraction in control region  " % top1j_simt_sub[0]
#print "N(top, 1j, Sim.) truth value: %8.2f                                                       " % yields['1j0b']['top'] 
#print "Scale factor, 1j: %5.3f " % (top1j_data_sub[0]/top1j_simt_sub[0])
#print ""
#

if len(args) > 2: print "cutfile mm err me err em err ee err"
for cutfile in args[2:]:
    mca.scaleProcess("Top", 1.17)
    report = mca.getYields(CutsFile(cutfile,options), process="Top")
    print cutfile, 
    print report['Top'][-1][1][MM][1], report['Top'][-1][1][MM][2], 
    print report['Top'][-1][1][ME][1], report['Top'][-1][1][ME][2], 
    print report['Top'][-1][1][EM][1], report['Top'][-1][1][EM][2], 
    print report['Top'][-1][1][EE][1], report['Top'][-1][1][EE][2], 
    print ""
    
