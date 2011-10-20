#!/usr/bin/env python

import re
from sys import argv
import os.path
from optparse import OptionParser
from math import sqrt,fabs
parser = OptionParser()
parser.add_option("-s", "--stat",   dest="stat",          default=False, action="store_true")  # ignore systematic uncertainties to consider statistical uncertainties only
parser.add_option("-S", "--force-shape", dest="shape",    default=False, action="store_true")  # ignore systematic uncertainties to consider statistical uncertainties only
parser.add_option("-a", "--asimov", dest="asimov",  default=False, action="store_true")

(options, args) = parser.parse_args()
options.bin = True # fake that is a binary output, so that we parse shape lines

from HiggsAnalysis.CombinedLimit.DatacardParser import *

DC = parseCard(file("datacards/hww-1.55fb.mH120.comb.txt"), options)

nuisToConsider = [ y for y in DC.systs if 'CMS' in y[0] ]

errors = {}
for nuis in nuisToConsider:
    if nuis[2] == 'gmN': gmN = nuis[3][0]
    else               : gmN = 0
    for channel in nuis[4]:
        if channel not in errors.keys(): errors[channel] = {}
        for process in nuis[4][channel]:
            if nuis[4][channel][process] == 0: continue
            if gmN != 0:
#                 newError = nuis[4][channel][process] * sqrt(gmN) / DC.exp[channel][process]
                newError = nuis[4][channel][process] * sqrt(gmN) / DC.exp[channel][process]
            else:
                newError = fabs(1-nuis[4][channel][process])
#             print process,nuis[4][channel][process],gmN,newError
            if process in errors[channel].keys():
                errors[channel][process] += newError*newError
            else:
                errors[channel][process] = newError*newError

for channel in errors:
    for process in errors[channel]:
        errors[channel][process] = sqrt(errors[channel][process])


for x in DC.exp:
    for y in DC.exp[x]:
        print "%10s %10s %10.2f +/- %10.2f (rel = %10.2f)" % (x,y,DC.exp[x][y],DC.exp[x][y]*errors[x][y],errors[x][y])


