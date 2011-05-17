#!/usr/bin/env python

import WWAnalysis.AnalysisStep.scaleFactors_cff
from optparse import OptionParser
from glob import glob
import re, os, os.path

parser = OptionParser(usage="%prog [options] block path")
parser.add_option("-p", "--pattern", dest="pattern", default='%(id)s.%(name)s/%(name)s_*.root', help="Pattern for file names");
parser.add_option("-j", "--json",    dest="json",    default='certifiedUCSD',                   help="JSON file, under WWAnalysis/Misc/Jsons")
(options, args) = parser.parse_args()

if len(args) < 2: raise RuntimeError, "usage: makeStep3.py [options] block path"
if not hasattr(WWAnalysis.AnalysisStep.scaleFactors_cff, args[0]): 
    raise RuntimeError, "Sample block %s not found in scaleFactors_cff" % args[0]
if not os.path.exists(args[1]): 
    raise RuntimeError, "Path %s does not exist" % args[1]
for id,list in getattr(WWAnalysis.AnalysisStep.scaleFactors_cff, args[0]).items(): 
    idn = re.sub('[^0-9]','',id)
    #print "ID %s: name %s" % (id, list[0]); continue
    pattern = args[1]+"/"+(options.pattern % {'name':list[0], 'id':id})
    arg3  = list[1] if len(list) == 2 else options.json
    #print "Pattern: ",pattern
    #if len(list) == 2: print "\tis MC, scale factor = ", list[1]
    #for f in glob(pattern): print "\t",f
    os.system("cmsSplit.pl step3.py %(dataset)s %(id)s %(arg3)s -a --bash --files=%(pattern)s --label=%(dataset)s --fj 1" % {
                'dataset':list[0], 'id':idn, 'arg3':arg3, 'pattern':pattern
              })
"""
foreach my $f (@files) {
    my $dataset = undef;
    $f =~ m{.*/\d+[AB]?\.(\w+?)/\1_.*.root} and $dataset = $1;
    $f =~ m{.*/\d+[AB]?/(\w+?)_.*.root} and $dataset = $1;
    die "Unrecognised filename '$f'\n" unless defined($dataset);
    $dataset =~ s/2011A(v2)?$//;
    $dataset =~ m/SingleElectron/ and next; ## skip this
    my $mass = '';
    if ($dataset =~ /^ggToH(\d+)to.*/) {
        $dataset = 'ggH'; $mass = $1;
    } elsif ($dataset =~ /^vbfToH(\d+)to.*/) {
        $dataset = 'vbfH'; $mass = $1;
    }
    $pmap{$dataset . $mass} = [] unless defined($pmap{$dataset . $mass});
    push @{$pmap{$dataset . $mass}}, $f;
}

foreach my $p (keys(%pmap)) {
    open FLIST, "> files.$p" or die "Can't write to files.$p";
    print FLIST join("\n", @{$pmap{$p}}, "\n");
    close FLIST;
    system "~/pl/cmsSplit.pl step3.py $p --args --label $p --filelist=files.$p --fj 1 --bash\n";
}
"""
