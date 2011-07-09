#!/usr/bin/env python

import WWAnalysis.AnalysisStep.scaleFactors_cff
from optparse import OptionParser
from glob import glob
import re, os, os.path

parser = OptionParser(usage="%prog [options] block [base path (to switch to pattern)]")
parser.add_option("-p", "--pattern", dest="pattern", default='%(id)s.%(name)s/*.root', help="Pattern for file names");
parser.add_option("-q", "--query", dest="query",  default='find file where dataset = %(dataset)s', help="Pattern for file names");
parser.add_option("-j", "--json",    dest="json",    default='certifiedUCSD',                   help="JSON file, under WWAnalysis/Misc/Jsons")
parser.add_option("-n", "--num",    dest="num",   type='int',  default=1,                   help="number of file per job")
(options, args) = parser.parse_args()

usePattern = False
if len(args) < 1: raise RuntimeError, "usage: makeStep3.py [options] block [base path (to switch to pattern)]"
if len(args) > 1: usePattern = True
if not hasattr(WWAnalysis.AnalysisStep.scaleFactors_cff, args[0]):
    print "Sample block %s not found in scaleFactors_cff" % args[0]
    for k in dir(WWAnalysis.AnalysisStep.scaleFactors_cff):
        print "\t",k
    raise RuntimeError, "Sample block %s not found in scaleFactors_cff" % args[0]
if usePattern and not os.path.exists(args[1]): 
    raise RuntimeError, "Path %s does not exist" % args[1]
for id,list in getattr(WWAnalysis.AnalysisStep.scaleFactors_cff, args[0]).items(): 
    idn = re.sub('[^0-9]','',id)
    #print "ID %s: name %s" % (id, list[0]); continue
    arg3  = list[1] if len(list) == 3 else options.json
    if usePattern: 
        pattern = args[1]+"/"+(options.pattern % {'name':list[0], 'id':id})
        os.system("cmsSplit.pl step3.py %(dataset)s %(id)s %(arg3)s -a --bash --files=%(pattern)s --label=%(id2)s_%(dataset)s --fj %(num)d" % {
                    'dataset':list[0], 'id':idn, 'arg3':arg3, 'pattern':pattern, 'num':options.num, 'id2':id
                  })
    else : 
        query = options.query % {'dataset':list[-1]}
        os.system('cmsSplit.pl step3.py %(dataset)s %(id)s %(arg3)s -a --bash --usePhys --dbsql="%(query)s" --label=%(id2)s_%(dataset)s --fj %(num)d' % {
                'dataset':list[0], 'id':idn, 'arg3':arg3, 'query':query, 'num':options.num, 'id2':id
              })
