#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>

#include "scales2.h"
#include "FakeRateCalculator.h"
#include "CardTemplate.h"
#include "FitMaker.h"
#include "YieldMaker.h"
#include "findAndReplace.h"

using namespace std;
using namespace RooFit;

double wavg(double x, double y, double ex, double ey) {
  return (x/ex/ex + y/ey/ey)/(1/ex/ex + 1/ey/ey);
}

pair<double,double> werror(double x, double y, double ex, double ey) {
  double max = TMath::Max(x+ex,y+ey);
  double min = TMath::Min(x-ex,y-ey);

  double val = wavg(x,y,ex,ey);
  double eu=fabs(max-val);
  double ed=fabs(val-min);
  
  return std::make_pair<double,double>(ed,eu);
}

struct QQZZ4muCardInfo {

  float mean_lan, sigma_lan, k_lnN, f1;

    QQZZ4muCardInfo() {
      mean_lan  = 0.0;
      sigma_lan = 0.0;
      k_lnN     = 0.0;
      f1        = 0.0;
    }


    QQZZ4muCardInfo(ZZEbE4muFitMaker fm) {
      mean_lan  = fm.getVarLanMean();
      sigma_lan = fm.getVarLanSigma();
      k_lnN     = fm.getVarLnNK();
      f1        = fm.getVarF1();
    }


    std::string printInfo() {
        
        stringstream info;
        info << "\t\t\t\tmLan = " << mean_lan  << ";" << std::endl;
        info << "\t\t\t\tsLan = " << sigma_lan << ";" << std::endl;
        info << "\t\t\t\tkLnN = " << k_lnN     << ";" << std::endl;
        info << "\t\t\t\tf1   = " << f1        << ";" << std::endl;
        info << "\n\n";
        return info.str();
    }
};

struct QQZZ4eAnd2e2muCardInfo {

  float mean_lan, sigma_lan, mean_gaus, sigma_gaus, f1;

    QQZZ4eAnd2e2muCardInfo() {
      mean_lan   = 0.0;
      sigma_lan  = 0.0;
      mean_gaus  = 0.0;
      sigma_gaus = 0.0;
      f1         = 0.0;
    }


    QQZZ4eAnd2e2muCardInfo(ZZEbE4eAnd2e2muFitMaker fm) {
      mean_lan   = fm.getVarLanMean();
      sigma_lan  = fm.getVarLanSigma();
      mean_gaus  = fm.getVarGauMean();
      sigma_gaus = fm.getVarGauSigma();
      f1         = fm.getVarF1();
    }


    std::string printInfo() {
        
        stringstream info;
        info << "\t\t\t\tqqmLan = " << mean_lan   << ";" << std::endl;
        info << "\t\t\t\tqqsLan = " << sigma_lan  << ";" << std::endl;
        info << "\t\t\t\tqqmLnN = " << mean_gaus  << ";" << std::endl;
        info << "\t\t\t\tqqkLnN = " << sigma_gaus << ";" << std::endl;
        info << "\t\t\t\tqqf1   = " << f1         << ";" << std::endl;
        info << "\n\n";
        return info.str();
    }
};

struct HiggsMassPointInfo {

    float z1min;
    float z2min;
    float massLow;
    float massHigh;
    float massErrLow;
    float massErrHigh;
    float melacut;
    bool  do7TeV;
    int   njetscut;
    bool  njetveto;
    float ptcut;
    float ptmax;

    ZXYieldMaker   ymaker_zxss;
    ZXYieldMaker   ymaker_zxos;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;

    void analyze(int ch, ofstream& file) {

        float yield_qq  = ymaker_qqzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        float yield_gg  = ymaker_ggzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        float yield_zjss  = ymaker_zxss.getYield(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        //        float yield_zjos  = ymaker_zxos.getYield(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);

        float yielde_qq  = ymaker_qqzz.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        float yielde_gg  = ymaker_ggzz.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        float yielde_zjss  = ymaker_zxss.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        //        float yielde_zjos  = ymaker_zxos.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);

        if (ch == 2) {
            yield_qq += ymaker_qqzz.getYield(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
            yield_gg += ymaker_ggzz.getYield(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
            yield_zjss += ymaker_zxss.getYield(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
            //            yield_zjos += ymaker_zxos.getYield(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);

            yielde_qq = sqrt(yielde_qq + pow(ymaker_qqzz.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax),2));
            yielde_gg = sqrt(yielde_gg + pow(ymaker_ggzz.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax),2));
            yielde_zjss = sqrt(yielde_zjss + pow(ymaker_zxss.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax),2));
            //            yielde_zjos = sqrt(yielde_zjos + pow(ymaker_zxos.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax),2));
        }

        std::string chstr;

        if(ch==0) {
          ZZEbE4muFitMaker  fitmaker_qz ("bkg_qqzz" , massLow, massHigh);
          ZZEbE4muFitMaker  fitmaker_gz ("bkg_ggzz" , massLow, massHigh);
          ZZEbE4muFitMaker  fitmaker_zxss ("bkg_zjets_ss", massLow, massHigh);
          
          fitmaker_gz.add(ymaker_ggzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          fitmaker_qz.add(ymaker_qqzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          fitmaker_zxss.add(ymaker_zxss.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          
          if (ch == 2) {
            fitmaker_gz.add(ymaker_ggzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
            fitmaker_qz.add(ymaker_qqzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
            fitmaker_zxss.add(ymaker_zxss.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          }        
          
          chstr = "4mu";
          
          std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
          /*
            fitmaker_zxss.fit();
            fitmaker_zxss.print("cat0/zjets_ss_"+chstr+"_"+tevstr+".pdf", 15);
            fitmaker_zxos.fit();
            fitmaker_zxos.print("cat0/zjets_os_"+chstr+"_"+tevstr+".pdf", 15);
            ZXCardInfo zci(fitmaker_zxss,fitmaker_zxos);
          */
          
          fitmaker_qz.fit();
          fitmaker_qz.print("cat0/qqZZ_"+chstr+tevstr, 200);
          fitmaker_qz.print("cat0/qqZZ_"+chstr+"_zoom"+tevstr, 100, true);
          QQZZ4muCardInfo qzci(fitmaker_qz);
          
          fitmaker_gz.fit();
          fitmaker_gz.print("cat0/ggZZ_"+chstr+tevstr, 50, true);
          fitmaker_gz.print("cat0/ggZZ_"+chstr+"_zoom"+tevstr, 50, true);
          QQZZ4muCardInfo gzci(fitmaker_gz);

          file << "\n\n\n***********************************************************";
          file << "Channel = " << chstr << "********************************************************************\n\n";
          //file << zci.printInfo() ;
          file << qzci.printInfo();
          file << gzci.printInfo();

        } else {
          if (ch == 1) chstr = "4e";
          if (ch == 2) chstr = "2e2mu";

          ZZEbE4eAnd2e2muFitMaker  fitmaker_qz ("bkg_qqzz" , massLow, massHigh);
          ZZEbE4eAnd2e2muFitMaker  fitmaker_gz ("bkg_ggzz" , massLow, massHigh);
          ZZEbE4eAnd2e2muFitMaker  fitmaker_zxss ("bkg_zjets_ss", massLow, massHigh);

//           float f1val = (ch == 1) ? 0.50 : 0.65;
//           fitmaker_gz.setVarF1(f1val);
//           fitmaker_qz.setVarF1(f1val);
          
          fitmaker_gz.add(ymaker_ggzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          fitmaker_qz.add(ymaker_qqzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          fitmaker_zxss.add(ymaker_zxss.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          
          if (ch == 2) {
            fitmaker_gz.add(ymaker_ggzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
            fitmaker_qz.add(ymaker_qqzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
            fitmaker_zxss.add(ymaker_zxss.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
          }        
          
          std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
          /*
            fitmaker_zxss.fit();
            fitmaker_zxss.print("cat0/zjets_ss_"+chstr+"_"+tevstr+".pdf", 15);
            fitmaker_zxos.fit();
            fitmaker_zxos.print("cat0/zjets_os_"+chstr+"_"+tevstr+".pdf", 15);
            ZXCardInfo zci(fitmaker_zxss,fitmaker_zxos);
          */
          
          fitmaker_qz.fit();
          fitmaker_qz.print("cat0/qqZZ_"+chstr+tevstr, 200);
          fitmaker_qz.print("cat0/qqZZ_"+chstr+"_zoom"+tevstr, 100, true);
          QQZZ4eAnd2e2muCardInfo qzci(fitmaker_qz);
          
          fitmaker_gz.fit();
          fitmaker_gz.print("cat0/ggZZ_"+chstr+tevstr, 50, true);
          fitmaker_gz.print("cat0/ggZZ_"+chstr+"_zoom"+tevstr, 50, true);
          QQZZ4eAnd2e2muCardInfo gzci(fitmaker_gz);

          file << "\n\n\n***********************************************************";
          file << "Channel = " << chstr << "********************************************************************\n\n";
          //file << zci.printInfo() ;
          file << qzci.printInfo();
          file << gzci.printInfo();
          
        }
      


//         file << "\n--------------- Yields ------------------ \n";
//         file << "\t\t\t\tqq->ZZ : " << yield_qq << " +/- " << yielde_qq << "\n";
//         file << "\t\t\t\tgg->ZZ : " << yield_gg << " +/- " << yielde_gg << "\n";
//         file << "\t\t\t\tZ+jets (SS) : " << yield_zjss << " +/- " << yielde_zjss << "\n";
//         file << "\t\t\t\tZ+jets (OS) : " << yield_zjos << " +/- " << yielde_zjos << "\n";
//         file << "\t\t\t\tZ+jets (Avg): " << wavg(yield_zjss,yield_zjos,yielde_zjss,yielde_zjos) 
//              << " + " << werror(yield_zjss,yield_zjos,yielde_zjss,yielde_zjos).second 
//              << " - " << werror(yield_zjss,yield_zjos,yielde_zjss,yielde_zjos).first << "\n";
 
    } 


};

void getParams() {

    float lumi7     = 5.05;
    float lumi8     = 19.6;

    std::string treeFolder7 = "root://pcmssd12//data/hzz4l/step2/HZZ4L_42X_S1_V12_S2_V03/";
    std::string treeFolder8 = "root://pcmssd12//data/hzz4l/step2/HZZ4L_53X_S1_V12_S2_V03/";

    enum zxalgo { ss=0, os=1};
    float Z1min[2] = { 40,   81.186 };
    float Z1max[2] = { 120, 101.186 };

    /*    
    HiggsMassPointInfo hmpi7;
    hmpi7.z1min    = 40.;
    hmpi7.z2min    = 12.;
    hmpi7.massLow  = 100.;
    hmpi7.massHigh = 160;
    hmpi7.melacut  = -1.;
    hmpi7.do7TeV   = true;
    hmpi7.njetscut = -1;
    hmpi7.ptcut    = 0;
    hmpi7.ptmax    = 1000000;
    hmpi7.njetveto = true;

    init(hmpi7.do7TeV);

    FakeRateCalculator FR_7TeV_SS(treeFolder7+"DATA/7TeV/yesRegrYesCalibYesMu/data2011.root", hmpi7.do7TeV, Z1min[ss], Z1max[ss], 0.0, 0.0, true);
    FakeRateCalculator FR_7TeV_OS(treeFolder7+"DATA/7TeV/yesRegrYesCalibYesMu/data2011.root", hmpi7.do7TeV, Z1min[os], Z1max[os], 0.0, 0.0, true);
    
    hmpi7.ymaker_zxss.fill(treeFolder7+"DATA/7TeV/yesRegrYesCalibYesMu/data2011.root"          , 1.0, FR_7TeV_SS, true);
    hmpi7.ymaker_zxos.fill(treeFolder7+"DATA/7TeV/yesRegrYesCalibYesMu/data2011.root"          , 1.0, FR_7TeV_OS, false);

    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id102.root" , getBkgXsec(102)*lumi7/evt_7TeV(102), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id103.root" , getBkgXsec(103)*lumi7/evt_7TeV(103), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id104.root" , getBkgXsec(104)*lumi7/evt_7TeV(104), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id105.root" , getBkgXsec(105)*lumi7/evt_7TeV(105), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id106.root" , getBkgXsec(106)*lumi7/evt_7TeV(106), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id107.root" , getBkgXsec(107)*lumi7/evt_7TeV(107), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id101.root" , getBkgXsec(101)*lumi7/evt_7TeV(101), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id100.root" , getBkgXsec(100)*lumi7/evt_7TeV(100), 0.0, false);

    ofstream file7;
    file7.open("cat0/info7TeV.txt");

    hmpi7.analyze(0, file7);
    hmpi7.analyze(1, file7);
    hmpi7.analyze(2, file7);
    
    file7.close();
    */

    HiggsMassPointInfo hmpi8;
    hmpi8.z1min    = 40.;
    hmpi8.z2min    = 12.;
    hmpi8.massLow  = 100;
    hmpi8.massHigh = 160;
    hmpi8.do7TeV   = false;
    hmpi8.melacut  = -1.;
    hmpi8.njetscut = -1;
    hmpi8.ptcut    = 0;
    hmpi8.ptmax    = 1000000.;
    hmpi8.njetveto = true;


    setmc(false);    
    init(hmpi8.do7TeV);

    FakeRateCalculator FR_8TeV_SS(treeFolder8+"DATA/8TeV/yesRegrYesCalibYesMu/hcp.root", hmpi8.do7TeV, Z1min[ss], Z1max[ss], 0.0, 0.0, true);
    FakeRateCalculator FR_8TeV_OS(treeFolder8+"DATA/8TeV/yesRegrYesCalibYesMu/hcp.root", hmpi8.do7TeV, Z1min[os], Z1max[os], 0.0, 0.0, true);
    
    hmpi8.ymaker_zxss.fill(treeFolder8+"DATA/8TeV/yesRegrYesCalibYesMu/hcp.root"          , 1.0, FR_8TeV_SS, true);
    //    hmpi8.ymaker_zxos.fill(treeFolder8+"DATA/8TeV/yesRegrYesCalibYesMu/hcp.root"          , 1.0, FR_8TeV_OS, false);

    setmc(true);    
    init(hmpi8.do7TeV);

    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id102.root" , getBkgXsec(102)*lumi8/evt_8TeV(102), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id103.root" , getBkgXsec(103)*lumi8/evt_8TeV(103), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id104.root" , getBkgXsec(104)*lumi8/evt_8TeV(104), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id105.root" , getBkgXsec(105)*lumi8/evt_8TeV(105), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id106.root" , getBkgXsec(106)*lumi8/evt_8TeV(106), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id107.root" , getBkgXsec(107)*lumi8/evt_8TeV(107), 0.0, false);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id101.root" , getBkgXsec(101)*lumi8/evt_8TeV(101), 0.0, false);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"MC/8TeV/yesRegrYesCalibYesMu/hzzTree_id100.root" , getBkgXsec(100)*lumi8/evt_8TeV(100), 0.0, false);

    ofstream file8;
    file8.open("cat0/info8TeV.txt");

    hmpi8.analyze(0, file8);
    hmpi8.analyze(1, file8);
    hmpi8.analyze(2, file8);
    
    file8.close();
}

