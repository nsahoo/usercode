import os, re
from math import *

def file2map(file):
    map = {}
    for line in open(file,"r"):
        fields = [x.strip() for x in line.split()]
        if len(fields) > 1 and re.match("\d+", fields[0]):
            map[float(fields[0])] = [float(y) for y in fields[1:] if re.match(r"-?\d+.*",y) ]
    return map

SYST_PATH = os.getenv('CMSSW_BASE')+'/src/WWAnalysis/AnalysisStep/test/datacard/'
YR_ggH  = file2map(SYST_PATH+'YR-ggH.txt')
YR_vbfH = file2map(SYST_PATH+'YR-vbfH.txt')

ggH_pdfErrYR  = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ggH.items()] )
ggH_scaErrYR  = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ggH.items()] )
vbfH_pdfErrYR = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_vbfH.items()] )
vbfH_scaErrYR = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_vbfH.items()] )

for X in 450, 550:
    ggH_pdfErrYR[X]  = 0.5*(ggH_pdfErrYR[X-10] +ggH_pdfErrYR[X+10])
    ggH_scaErrYR[X]  = 0.5*(ggH_scaErrYR[X-10] +ggH_scaErrYR[X+10])
    vbfH_pdfErrYR[X] = 0.5*(vbfH_pdfErrYR[X-10]+vbfH_pdfErrYR[X+10])
    vbfH_scaErrYR[X] = 0.5*(vbfH_scaErrYR[X-10]+vbfH_scaErrYR[X+10])

ggH_jets = dict([(m, dict(zip(['f0','f1','f2','k1','k2'], vals))) for m,vals in file2map(SYST_PATH+"ggH_jetBins.txt").items()]) 


def getCommonSysts(mass,channel,jets,qqWWfromData):
    nuisances = {} 
    MCPROC = ['ggH', 'vbfH', 'DTT', 'ggWW', 'VV', 'Vg' ]; 
    if channel == 'elmu' or channel == 'muel': MCPROC+=['DYMM','DYEE']
    if not qqWWfromData: MCPROC.append('WW')
    # -- Luminosity ---------------------
    nuisances['lumi'] = [ ['lnN'], dict([(p,1.045) for p in MCPROC])]
    # -- PDF ---------------------
    nuisances['pdf_gg']    = [ ['lnN'], { 'ggH':ggH_pdfErrYR[mass], 'ggWW':1.04 }]
    nuisances['pdf_qqbar'] = [ ['lnN'], { 'vbfH':vbfH_pdfErrYR[mass], 'VV':1.04, 'WW':(1.0 if mass < 200 else 1.04) }]
    # -- Theory ---------------------
    if jets == 0:
        # appendix D of https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=0&confId=135333
        k0 = pow(ggH_scaErrYR[mass],     1/ggH_jets[mass]['f0'])
        k1 = pow(ggH_jets[mass]['k1'], 1-1/ggH_jets[mass]['f0']) # -f1-f2=f0-1
        nuisances['QCDscale_ggH']    = [  ['lnN'], { 'ggH':k0 }]
        nuisances['QCDscale_ggH1in'] = [  ['lnN'], { 'ggH':k1 }]
    elif jets == 1:
        k1 = pow(ggH_jets[mass]['k1'], 1+ggH_jets[mass]['f2']/ggH_jets[mass]['f1'])
        k2 = pow(ggH_jets[mass]['k2'],  -ggH_jets[mass]['f2']/ggH_jets[mass]['f1'])
        nuisances['QCDscale_ggH1in'] = [  ['lnN'], { 'ggH':k1 }]
        nuisances['QCDscale_ggH2in'] = [  ['lnN'], { 'ggH':k2 }]
#     nuisances['QCDscale_qqH']  = [ ['lnN'], { 'vbfH':vbfH_scaErrYR[mass] }]
    nuisances['QCDscale_qqH']  = [ ['lnN'], { 'vbfH':1.01 }]
    nuisances['QCDscale_ggVV'] = [ ['lnN'], { 'ggWW':1.5}]
#     nuisances['QCDscale_VV']   = [ ['lnN'], { 'VV':1.03, 'WW':(1.0 if qqWWfromData else 1.03) }]
    nuisances['QCDscale_VV']   = [ ['lnN'], { 'VV':1.04 }]
    nuisances['QCDscale_Vg'] = [ ['lnN'], {'Vg':1.50}]
    # -- Experimental ---------------------
    if 'm' in channel: nuisances['CMS_eff_m'] = [ ['lnN'], dict([(p,pow(1.02,channel.count('m'))) for p in MCPROC])]
    if 'e' in channel: nuisances['CMS_eff_e'] = [ ['lnN'], dict([(p,pow(1.02,channel.count('e'))) for p in MCPROC if p != 'Vg'])]
    if   channel == 'mumu': nuisances['CMS_p_scale_m'] = [ ['lnN'], dict([(p,1.015) for p in MCPROC if p != 'DTT'] )]
    elif channel == 'elmu': nuisances['CMS_p_scale_m'] = [ ['lnN'], dict([(p,1.015) for p in MCPROC if p != 'DTT'] )]
    elif channel == 'muel': nuisances['CMS_p_scale_e'] = [ ['lnN'], dict([(p,1.020) for p in MCPROC if p != 'DTT'] )]
    elif channel == 'elel': nuisances['CMS_p_scale_e'] = [ ['lnN'], dict([(p,1.020) for p in MCPROC if p != 'DTT'] )]
    nuisances['CMS_met'] = [ ['lnN'], dict([(p,1.02) for p in MCPROC if p != 'DTT'])]
#     if (jets == 0): nuisances['CMS_p_scale_j'] = [ ['lnN'], {'ggH':0.97}]
#     else:           nuisances['CMS_p_scale_j'] = [ ['lnN'], {'ggH':1.01}]
    nuisances['CMS_scale_j'] = [ ['lnN'], dict([(p,1.02) for p in MCPROC if p != 'DTT'])]
    nuisances['QCDscale_ggH_ACEPT'] = [ ['lnN'], {'ggH':1.02}]
    nuisances['QCDscale_qqH_ACEPT'] = [ ['lnN'], {'vbfH':1.02}]
    if (jets == 0): nuisances['UEPS'] = [ ['lnN'], {'ggH':0.94}]
    else:           nuisances['UEPS'] = [ ['lnN'], {'ggH':1.11}]
    if (jets == 0): nuisances['CMS_QCDscale_WW_EXTRAP'] = [ ['lnN'], {'WW':0.954}]
    else:           nuisances['CMS_QCDscale_WW_EXTRAP'] = [ ['lnN'], {'WW':1.206}]
    return nuisances

def addFakeBackgroundSysts(nuisances, mass,channel,jets,errWW=0.2,errDY=1.0,errTop0j=1.0,errTop1j=0.3,errWJ=0.5):
    if errWW:  nuisances['CMS_norm_WW'           ] = [ ['lnN'], { 'WW':(1+errWW)} ] 
    if errTop0j and jets == 0: nuisances['CMS_norm_Top0j'] = [ ['lnN'], { 'Top':(1+errTop0j)}] 
    if errTop1j and jets == 1: nuisances['CMS_norm_Top1j'] = [ ['lnN'], { 'Top':(1+errTop1j)}] 
    if errWJ: nuisances['CMS_fake_%s'%channel[2]] = [ ['lnN'], { 'WJet':(1+errWJ) } ] 
    if errDY and channel in ['mumu', 'elel']:
        proc = "DY{0}{0}".format(channel[0].upper())
        nuisances['CMS_norm_DY_%s' % channel[0]]  = [ ['lnN'], { proc:(1+errDY) } ] 

