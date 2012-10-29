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

using namespace RooFit;

struct QQZZCardInfo {

    float a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13;

    QQZZCardInfo() {
        a0  = 0.0;
        a1  = 0.0;
        a2  = 0.0;
        a3  = 0.0;
        a4  = 0.0;
        a5  = 0.0;
        a6  = 0.0;
        a7  = 0.0;
        a8  = 0.0;
        a9  = 0.0;
        a10 = 0.0;
        a11 = 0.0;
        a12 = 0.0;
        a13 = 0.0;
    }


    QQZZCardInfo(QQZZFitMaker fm) {
        a0  = fm.getVarA0();
        a1  = fm.getVarA1();
        a2  = fm.getVarA2();
        a3  = fm.getVarA3();
        a4  = fm.getVarA4();
        a5  = fm.getVarA5();
        a6  = fm.getVarA6();
        a7  = fm.getVarA7();
        a8  = fm.getVarA8();
        a9  = fm.getVarA9();
        a10 = fm.getVarA10();
        a11 = fm.getVarA11();
        a12 = fm.getVarA12();
        a13 = fm.getVarA13();
    }


    std::string printInfo() {
        
        stringstream info;
        info << "\t\t\t\tqqa0  = " << a0  << ";" << std::endl;
        info << "\t\t\t\tqqa1  = " << a1  << ";" << std::endl;
        info << "\t\t\t\tqqa2  = " << a2  << ";" << std::endl;
        info << "\t\t\t\tqqa3  = " << a3  << ";" << std::endl;
        info << "\t\t\t\tqqa4  = " << a4  << ";" << std::endl;
        info << "\t\t\t\tqqa5  = " << a5  << ";" << std::endl;
        info << "\t\t\t\tqqa6  = " << a6  << ";" << std::endl;
        info << "\t\t\t\tqqa7  = " << a7  << ";" << std::endl;
        info << "\t\t\t\tqqa8  = " << a8  << ";" << std::endl;
        info << "\t\t\t\tqqa9  = " << a9  << ";" << std::endl;
        info << "\t\t\t\tqqa10 = " << a10 << ";" << std::endl;
        info << "\t\t\t\tqqa11 = " << a11 << ";" << std::endl;
        info << "\t\t\t\tqqa12 = " << a12 << ";" << std::endl;
        info << "\t\t\t\tqqa13 = " << a13 << ";" << std::endl;
        info << "\n\n";

        return info.str();
    }
};

struct GGZZCardInfo {

    float a0, a1, a2, a3, a4, a5, a6, a7, a8, a9;

    GGZZCardInfo() {
        a0  = 0.0;
        a1  = 0.0;
        a2  = 0.0;
        a3  = 0.0;
        a4  = 0.0;
        a5  = 0.0;
        a6  = 0.0;
        a7  = 0.0;
        a8  = 0.0;
        a9  = 0.0;
    }


    GGZZCardInfo(GGZZFitMaker fm) {
        a0  = fm.getVarA0();
        a1  = fm.getVarA1();
        a2  = fm.getVarA2();
        a3  = fm.getVarA3();
        a4  = fm.getVarA4();
        a5  = fm.getVarA5();
        a6  = fm.getVarA6();
        a7  = fm.getVarA7();
        a8  = fm.getVarA8();
        a9  = fm.getVarA9();
    }

    std::string printInfo() {
  
        stringstream info;
        info << "\t\t\t\tgga0  = " << a0  << ";" << std::endl;
        info << "\t\t\t\tgga1  = " << a1  << ";" << std::endl;
        info << "\t\t\t\tgga2  = " << a2  << ";" << std::endl;
        info << "\t\t\t\tgga3  = " << a3  << ";" << std::endl;
        info << "\t\t\t\tgga4  = " << a4  << ";" << std::endl;
        info << "\t\t\t\tgga5  = " << a5  << ";" << std::endl;
        info << "\t\t\t\tgga6  = " << a6  << ";" << std::endl;
        info << "\t\t\t\tgga7  = " << a7  << ";" << std::endl;
        info << "\t\t\t\tgga8  = " << a8  << ";" << std::endl;
        info << "\t\t\t\tgga9  = " << a9  << ";" << std::endl;
        info << "\n\n";

        return info.str();
    }

};

struct ZXCardInfo {

    float mean, sigma;

    ZXCardInfo() {
        mean   = 0.0;
        sigma  = 0.0;
    }

    ZXCardInfo(ZXFitMaker fm) {
        mean   = fm.getVarMean()  ;
        sigma  = fm.getVarSigma() ;
    }

    std::string printInfo() {
  
        stringstream info;
        info << "\t\t\t\tzxme  = " << mean   << ";" << std::endl;
        info << "\t\t\t\tzxsi = " << sigma  << ";" << std::endl;
        info << "\n\n";

        return info.str();
    }

};

struct HiggsMassPointInfo {

    float z1min;
    float z2min;
    float massLow;
    float massHigh;
    float melacut;
    bool  do7TeV;
    int   njetscut;
    bool  njetveto;
    float ptcut;
    float ptmax;

    ZXYieldMaker   ymaker_zxss;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;

    void analyze(int ch, ofstream& file) {

        float yield_qq  = ymaker_qqzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        float yield_gg  = ymaker_ggzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        float yield_zj  = ymaker_zxss.getYield(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);

        if (ch == 2) {
            yield_qq += ymaker_qqzz.getYield(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
            yield_gg += ymaker_ggzz.getYield(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
            yield_zj += ymaker_zxss.getYield(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax);
        }



        QQZZFitMaker   fitmaker_qz ("bkg_qqzz" , massLow, massHigh);
        GGZZFitMaker   fitmaker_gz ("bkg_ggzz" , massLow, massHigh);
        ZXFitMaker     fitmaker_zx ("bkg_zjets", massLow, massHigh);

        fitmaker_gz.add(ymaker_ggzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
        fitmaker_qz.add(ymaker_qqzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
        fitmaker_zx.add(ymaker_zxss.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));

        if (ch == 2) {
            fitmaker_gz.add(ymaker_ggzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
            fitmaker_qz.add(ymaker_qqzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
            fitmaker_zx.add(ymaker_zxss.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut, njetscut, njetveto, ptcut, ptmax));
        }        

        std::string chstr;
        if (ch == 0) chstr = "4mu";
        if (ch == 1) chstr = "4e";
        if (ch == 2) chstr = "2e2mu";

        std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
        fitmaker_zx.fit();
        fitmaker_zx.print("cat0/zjets_"+chstr+"_"+tevstr+".pdf", 15);
        ZXCardInfo zci(fitmaker_zx);

        fitmaker_qz.fit();
        fitmaker_qz.print("cat0/qqZZ_"+chstr+"_"+tevstr+".pdf", 70);
        QQZZCardInfo qzci(fitmaker_qz);

        fitmaker_gz.fit();
        fitmaker_gz.print("cat0/ggZZ_"+chstr+"_"+tevstr+".pdf", 70);
        GGZZCardInfo gzci(fitmaker_gz);

      
        file << "\n\n\n***********************************************************";
        file << "Channel = " << chstr << "********************************************************************\n\n";
        file << zci.printInfo() ;
        file << qzci.printInfo();
        file << gzci.printInfo();


        file << "\n--------------- Yields ------------------ \n";
        file << "\t\t\t\tqq->ZZ : " << yield_qq << "\n";
        file << "\t\t\t\tgg->ZZ : " << yield_gg << "\n";
        file << "\t\t\t\tZ+jets : " << yield_zj << "\n";

 
    } 


};

void getParams() {

    float lumi7     = 5.05;
    float lumi8     = 12.2;

    std::string treeFolder7 = "/home/avartak/CMS/Higgs/HCP/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";
    std::string treeFolder8 = "/home/avartak/CMS/Higgs/HCP/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/";

    bool doSS = true;
    float Z1min = doSS ? 40  : 81.186;
    float Z1max = doSS ? 120 : 101.186;
    
    HiggsMassPointInfo hmpi7;
    hmpi7.z1min    = 40.;
    hmpi7.z2min    = 12.;
    hmpi7.massLow  = 100.;
    hmpi7.massHigh = 1600;
    hmpi7.melacut  = -1.;
    hmpi7.do7TeV   = true;
    hmpi7.njetscut = -1;
    hmpi7.ptcut    = 0;
    hmpi7.ptmax    = 1000000;
    hmpi7.njetveto = true;

    init(hmpi7.do7TeV);

    FakeRateCalculator FR_7TeV(treeFolder7+"data.root", hmpi7.do7TeV, Z1min, Z1max, 0.0, 0.0, true);
    
    hmpi7.ymaker_zxss.fill(treeFolder7+"data.root"          , 1.0, FR_7TeV, doSS);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id102.root" , getBkgXsec(102)*lumi7/evt_7TeV(102), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id103.root" , getBkgXsec(103)*lumi7/evt_7TeV(103), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id104.root" , getBkgXsec(104)*lumi7/evt_7TeV(104), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id105.root" , getBkgXsec(105)*lumi7/evt_7TeV(105), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id106.root" , getBkgXsec(106)*lumi7/evt_7TeV(106), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id107.root" , getBkgXsec(107)*lumi7/evt_7TeV(107), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"hzzTree_id101.root" , getBkgXsec(101)*lumi7/evt_7TeV(101), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"hzzTree_id100.root" , getBkgXsec(100)*lumi7/evt_7TeV(100), 0.0, false);

    ofstream file7;
    file7.open("cat0/info7TeV.txt");

    hmpi7.analyze(0, file7);
    hmpi7.analyze(1, file7);
    hmpi7.analyze(2, file7);
    
    file7.close();

    HiggsMassPointInfo hmpi8;
    hmpi8.z1min    = 40.;
    hmpi8.z2min    = 12.;
    hmpi8.massLow  = 100;
    hmpi8.massHigh = 1600;
    hmpi8.do7TeV   = false;
    hmpi8.melacut  = -1.;
    hmpi8.njetscut = -1;
    hmpi8.ptcut    = 0;
    hmpi8.ptmax    = 1000000.;
    hmpi8.njetveto = true;

    init(hmpi8.do7TeV);

    FakeRateCalculator FR_8TeV(treeFolder8+"data.root", hmpi8.do7TeV, Z1min, Z1max, 0.0, 0.0, true);
    
    hmpi8.ymaker_zxss.fill(treeFolder8+"data.root"          , 1.0, FR_8TeV, doSS);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id102.root" , getBkgXsec(102)*lumi8/evt_8TeV(102), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id103.root" , getBkgXsec(103)*lumi8/evt_8TeV(103), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id104.root" , getBkgXsec(104)*lumi8/evt_8TeV(104), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id105.root" , getBkgXsec(105)*lumi8/evt_8TeV(105), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id106.root" , getBkgXsec(106)*lumi8/evt_8TeV(106), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id107.root" , getBkgXsec(107)*lumi8/evt_8TeV(107), 0.0, false);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"hzzTree_id101.root" , getBkgXsec(101)*lumi8/evt_8TeV(101), 0.0, false);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"hzzTree_id100.root" , getBkgXsec(100)*lumi8/evt_8TeV(100), 0.0, false);

    ofstream file8;
    file8.open("cat0/info8TeV.txt");

    hmpi8.analyze(0, file8);
    hmpi8.analyze(1, file8);
    hmpi8.analyze(2, file8);
    
    file8.close();

}


