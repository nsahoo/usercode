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
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
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
        info << "\n---------------- qq->ZZ shape parameters --------------------------\n";
        info << "\t\t\t\ta0  = " << a0  << ";" << std::endl;
        info << "\t\t\t\ta1  = " << a1  << ";" << std::endl;
        info << "\t\t\t\ta2  = " << a2  << ";" << std::endl;
        info << "\t\t\t\ta3  = " << a3  << ";" << std::endl;
        info << "\t\t\t\ta4  = " << a4  << ";" << std::endl;
        info << "\t\t\t\ta5  = " << a5  << ";" << std::endl;
        info << "\t\t\t\ta6  = " << a6  << ";" << std::endl;
        info << "\t\t\t\ta7  = " << a7  << ";" << std::endl;
        info << "\t\t\t\ta8  = " << a8  << ";" << std::endl;
        info << "\t\t\t\ta9  = " << a9  << ";" << std::endl;
        info << "\t\t\t\ta10 = " << a10 << ";" << std::endl;
        info << "\t\t\t\ta11 = " << a11 << ";" << std::endl;
        info << "\t\t\t\ta12 = " << a12 << ";" << std::endl;
        info << "\t\t\t\ta13 = " << a13 << ";" << std::endl;
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
        info << "\n---------------- gg->ZZ shape parameters --------------------------\n";
        info << "\t\t\t\ta0  = " << a0  << ";" << std::endl;
        info << "\t\t\t\ta1  = " << a1  << ";" << std::endl;
        info << "\t\t\t\ta2  = " << a2  << ";" << std::endl;
        info << "\t\t\t\ta3  = " << a3  << ";" << std::endl;
        info << "\t\t\t\ta4  = " << a4  << ";" << std::endl;
        info << "\t\t\t\ta5  = " << a5  << ";" << std::endl;
        info << "\t\t\t\ta6  = " << a6  << ";" << std::endl;
        info << "\t\t\t\ta7  = " << a7  << ";" << std::endl;
        info << "\t\t\t\ta8  = " << a8  << ";" << std::endl;
        info << "\t\t\t\ta9  = " << a9  << ";" << std::endl;
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
        info << "\n---------------- Z+jets shape parameters --------------------------\n";
        info << "\t\t\t\tmean  = " << mean   << ";" << std::endl;
        info << "\t\t\t\tsigma = " << sigma  << ";" << std::endl;
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
    bool do7TeV;

    ZXYieldMaker   ymaker_zxss;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;

    void analyze(int ch, ofstream& file) {

        float yield_qq  = ymaker_qqzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
        float yield_gg  = ymaker_ggzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
        float yield_zj  = ymaker_zxss.getYield(ch, z1min, z2min, massLow, massHigh, melacut);

        if (ch == 2) {
            yield_qq += ymaker_qqzz.getYield(3, z1min, z2min, massLow, massHigh, melacut);
            yield_gg += ymaker_ggzz.getYield(3, z1min, z2min, massLow, massHigh, melacut);
            yield_zj += ymaker_zxss.getYield(3, z1min, z2min, massLow, massHigh, melacut);
        }



        QQZZFitMaker   fitmaker_qz ("bkg_qqzz" , massLow, massHigh);
        GGZZFitMaker   fitmaker_gz ("bkg_ggzz" , massLow, massHigh);
        ZXFitMaker     fitmaker_zx ("bkg_zjets", massLow, massHigh);

        fitmaker_gz.add(ymaker_ggzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
        fitmaker_qz.add(ymaker_qqzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
        fitmaker_zx.add(ymaker_zxss.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));

        if (ch == 2) {
            fitmaker_gz.add(ymaker_ggzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut));
            fitmaker_qz.add(ymaker_qqzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut));
            fitmaker_zx.add(ymaker_zxss.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut));
        }        

        std::string chstr;
        if (ch == 0) chstr = "4mu";
        if (ch == 1) chstr = "4e";
        if (ch == 2) chstr = "2e2mu";

        std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
        fitmaker_zx.fit();
        fitmaker_zx.print("zjets_"+chstr+"_"+tevstr+".pdf", 15);
        ZXCardInfo zci(fitmaker_zx);

        fitmaker_qz.fit();
        fitmaker_qz.print("qqZZ_"+chstr+"_"+tevstr+".pdf");
        QQZZCardInfo qzci(fitmaker_qz);

        fitmaker_gz.fit();
        fitmaker_gz.print("ggZZ_"+chstr+"_"+tevstr+".pdf");
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
    float lumi8     = 5.26;

    std::string treeFolder7 = "/home/avartak/CMS/Higgs/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/";
    std::string treeFolder8 = "/home/avartak/CMS/Higgs/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/";

    bool doSS = true;

    /*
    HiggsMassPointInfo hmpi7;
    hmpi7.z1min    = 40.;
    hmpi7.z2min    = 12.;
    hmpi7.massLow  = 100.;
    hmpi7.massHigh = 800.;
    hmpi7.melacut  = -1.;
    hmpi7.do7TeV   = true;

    init(hmpi7.do7TeV);

    FakeRateCalculator FR_7TeV(treeFolder7+"hzzTree.root", hmpi7.do7TeV, 40, 120, 0.0, 0.0, true);
    
    hmpi7.ymaker_zxss.fill(treeFolder7+"hzzTree.root"       , 1.0, FR_7TeV, doSS);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id121.root" , getBkgXsec(121)*lumi7/evt_7TeV(121), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id122.root" , getBkgXsec(122)*lumi7/evt_7TeV(122), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id123.root" , getBkgXsec(123)*lumi7/evt_7TeV(123), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id124.root" , getBkgXsec(124)*lumi7/evt_7TeV(124), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id125.root" , getBkgXsec(125)*lumi7/evt_7TeV(125), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"hzzTree_id126.root" , getBkgXsec(126)*lumi7/evt_7TeV(126), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"hzzTree_id101.root" , getBkgXsec(101)*lumi7/evt_7TeV(101), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"hzzTree_id100.root" , getBkgXsec(100)*lumi7/evt_7TeV(100), 0.0, false);

    ofstream file7;
    file7.open("info7TeV.txt");

    hmpi7.analyze(0, file7);
    hmpi7.analyze(1, file7);
    hmpi7.analyze(2, file7);
    
    file7.close();
    */        

    HiggsMassPointInfo hmpi8;
    hmpi8.z1min    = 40.;
    hmpi8.z2min    = 12.;
    hmpi8.massLow  = 100.;
    hmpi8.massHigh = 800.;
    hmpi8.do7TeV   = false;

    init(hmpi8.do7TeV);

    FakeRateCalculator FR_8TeV(treeFolder8+"hcp.root", hmpi8.do7TeV, 40, 120, 0.0, 0.0, true);
    
    hmpi8.ymaker_zxss.fill(treeFolder8+"hcp.root"           , 1.0, FR_8TeV, doSS);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id102.root" , getBkgXsec(102)*lumi8/evt_8TeV(102), 0.0, false, 1);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id103.root" , getBkgXsec(103)*lumi8/evt_8TeV(103), 0.0, false, 1);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id104.root" , getBkgXsec(104)*lumi8/evt_8TeV(104), 0.0, false, 1);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id105.root" , getBkgXsec(105)*lumi8/evt_8TeV(105), 0.0, false, 1);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id106.root" , getBkgXsec(106)*lumi8/evt_8TeV(106), 0.0, false, 1);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"hzzTree_id107.root" , getBkgXsec(107)*lumi8/evt_8TeV(107), 0.0, false, 1);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"hzzTree_id101.root" , getBkgXsec(101)*lumi8/evt_8TeV(101), 0.0, false, 2);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"hzzTree_id100.root" , getBkgXsec(100)*lumi8/evt_8TeV(100), 0.0, false, 2);

    ofstream file8;
    file8.open("info8TeV.txt");

    hmpi8.analyze(0, file8);
    hmpi8.analyze(1, file8);
    hmpi8.analyze(2, file8);
    
    file8.close();


    std::cout << getBkgXsec(102)*lumi8/evt_8TeV(102) << std::endl;
    std::cout << getBkgXsec(103)*lumi8/evt_8TeV(103) << std::endl;
    std::cout << getBkgXsec(104)*lumi8/evt_8TeV(104) << std::endl;
    std::cout << getBkgXsec(105)*lumi8/evt_8TeV(105) << std::endl;
    std::cout << getBkgXsec(106)*lumi8/evt_8TeV(106) << std::endl;
    std::cout << getBkgXsec(107)*lumi8/evt_8TeV(107) << std::endl;
    std::cout << getBkgXsec(101)*lumi8/evt_8TeV(101) << std::endl;
    std::cout << getBkgXsec(100)*lumi8/evt_8TeV(100) << std::endl;
       
}


