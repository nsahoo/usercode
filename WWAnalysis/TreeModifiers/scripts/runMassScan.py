#!/usr/bin/env python

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

    (opt, args) = parser.parse_args()

    dcdir  = os.getcwd()+'/cards/'
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

    channels = ['4mu','4e','2e2mu','2mu2e','comb']
    for ch in channels:
        print "submitting toys for channel "+ch+"..."
        command = 'combine -M MultiDimFit '
        command = command+dcdir+'FloatMass_'+ch+'_hzz.root '
        command = command+'-m 125.8 -P MH -P r --floatOtherPOI=1 --algo=grid -n SCAN_'+ch+' --points='+str(opt.numPoints)
        if(opt.fastScan): command += ' --fastScan '

        jobs = int(opt.numPoints) / int(opt.pointsPerJob)
        firstPoint=0
        lastPoint=0
        for j in range(jobs):
            lastPoint += int(opt.pointsPerJob)
            runfile = srcdir+'scan-'+ch+'-j'+str(j)+'.src'
            f = open(runfile, 'w')
            f.write('cd ~/workspace/hzz4l/CMSSW_6_1_1/\n')
            f.write('eval `scram ru -sh` \n')
            f.write('cd - \n')
            extraflags = ' --firstPoint='+str(firstPoint)+' --lastPoint='+str(lastPoint)
            f.write(command+extraflags+'\n')
            f.write('mv higgsCombineSCAN_%s.MultiDimFit.mH*.root %sscan-%s-j%d.root \n' % (ch,outdir,ch,j) )
            firstPoint = lastPoint+1
            f.close()
            bsub = 'bsub -q '+opt.queue+' -J '+ch+'-j'+str(j)+' -o '+logdir+'/job-ch'+ch+'-j'+str(j)+'.log source '+runfile
            print '    job # '+str(j)
            #print bsub
            os.system(bsub)

    print 'Done. Wait and bye bye.'



if __name__ == '__main__':
    main()
        
