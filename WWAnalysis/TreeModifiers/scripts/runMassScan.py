#!/usr/bin/env python
# mass MH scan: ./runMassScan.py -c -n 500 -p 10 -m 125.6 --model=Mass -q 1nh
# mass MH mu scan: ./runMassScan.py -c -n 90000 -p 20 -m 125.6 --model=Mass --scanMu -q 8nh  
# width scan: ./runMassScan.py -c -n 200 -p 1 -m 125.6 --model=Width -q 1nd

import sys
import os.path
import optparse
import fnmatch

def main():
    usage = '''usage: %prog [opts] scenario'''
    parser = optparse.OptionParser(usage)
    parser.add_option('-a','--config1D',dest='oneD',help='Run 1D fit (m4l) ', action='store_true', default=False)
    parser.add_option('-b','--config2D',dest='twoD',help='Run 2D fit (m4l-merr) ', action='store_true', default=False)
    parser.add_option('-c','--config3D',dest='threeD',help='Run 1D fit (m4l-merr-KD) ', action='store_true', default=False)
    parser.add_option('--queue','-q',dest='queue',help='run in batch in queue specified as option (default -q 8nh)', default='8nh')
    parser.add_option('-n','--numPoints',dest='numPoints',help='Define the number of points of the grid',default=1000)
    parser.add_option('-p','--pointsPerJob',dest='pointsPerJob',help='Define the number of points to be run for each job',default=100)
    parser.add_option('-m','--mass',dest='mass',help='Define the central mass value',default=125.8)
    parser.add_option('-f', '--fast'     , dest='fastScan'    , help='do a fast Scan'                        , default=False   )
    parser.add_option('--model'          , type="string"      , dest='Model'       , help='default model is Mass, set to Width to run width scan' , default='Mass' )
    parser.add_option('--scanMu'         , dest='scanMu'      , help='do a 2D scan including mu' , default=False, action="store_true" )

    (opt, args) = parser.parse_args()

    dcdir  = os.getcwd()+'/cards_mass/' if opt.Model=='Mass' else os.getcwd()+'/cards_width/'
    outdir = os.getcwd()+'/results/'
    logdir = os.getcwd()+'/log/'
    srcdir = os.getcwd()+'/Src/'
    
    onedext   ='c1DNoMassErr'
    twodext   ='c1DMassErr'
    threedext ='c2DMassErr'

    ext=''
    if(opt.oneD):
        print 'you chose 1D'
        ext='c1DNoMassErr/'
    elif(opt.twoD):
        print 'you chose 2D'
        ext='c1DMassErr/'
    elif(opt.threeD):
        print 'you chose 3D'
        ext='c2DMassErr/'
    else:
        print 'You have to use 1D, 2D or 3D fit. Exiting.' 
        return

    dcdir += ext
    outdir += ext
    logdir += ext
    srcdir += ext

    if opt.Model=='Mass':
        ndscan = 'MH' 
        if opt.scanMu == True: ndscan='MHMu'
    elif opt.Model=='Width':
        ndscan = 'Width'
    else:
        print "ERROR! Model should be Mass or Width"

    channels = ['4mu','4e','2e2mu','comb']
    if opt.scanMu or opt.Model=='Width': channels = ['comb']
    for ch in channels:
        print "submitting toys for channel "+ch+"..."
        command = 'combine -M MultiDimFit '
        command = command+dcdir+'Float'+opt.Model+'_'+ch+'_hzz.root '
        if opt.Model=='Mass':    command += '-P MH '
        elif opt.Model=='Width': command += '-P HiggsDecayWidth '
        if opt.scanMu: command += '-P r --setPhysicsModelParameterRanges r=0,3 '
        command += '-m 125.6 --floatOtherPOI=1 --algo=grid -n SCAN'+opt.Model+'_'+ch+' --points='+str(opt.numPoints)
        if(opt.fastScan): command += ' --fastScan '

        jobs = int(opt.numPoints) / int(opt.pointsPerJob)
        firstPoint=0
        lastPoint=0
        for j in range(jobs):
            lastPoint += int(opt.pointsPerJob)
            runfile = srcdir+'scan-'+ndscan+'-'+ch+'-j'+str(j)+'.src'
            f = open(runfile, 'w')
            f.write('cd ~/workspace/hzz4l/CMSSW_6_1_1/\n')
            f.write('eval `scram ru -sh` \n')
            f.write('cd - \n')
            extraflags = ' --firstPoint='+str(firstPoint)+' --lastPoint='+str(lastPoint)
            f.write(command+extraflags+'\n')
            f.write('mv higgsCombineSCAN'+opt.Model+'_%s.MultiDimFit.mH*.root %sscan%s-%s-%s-j%d.root \n' % (ch,outdir,opt.Model,ndscan,ch,j) )
            firstPoint = lastPoint+1
            f.close()
            extraflags = ' -M 2000000 ' if opt.Model=='Width' else ' '
            bsub = 'bsub -q '+opt.queue+extraflags+' -J '+ch+'-j'+str(j)+' -o '+logdir+'/job-'+ndscan+'-ch'+ch+'-j'+str(j)+'.log source '+runfile
            print '    job # '+str(j)
            #print bsub
            os.system(bsub)

    print 'Done. Wait and bye bye.'



if __name__ == '__main__':
    main()
        
