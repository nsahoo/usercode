#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

//#if !defined(__CINT__) || !defined(__MAKECINT__)
//Headers for the data items
//#include <WWAnalysis/DataFormats/interface/SkimEvent.h>
//#endif

#include <WWAnalysis/DataFormats/interface/SkimEvent.h>


#include <TSystem.h>
#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/VectorUtil.h"


#include <sstream>
#include <vector>
#include <map>
//using namespace reco;

using namespace std;

void setInput(std::string inputFolder,std::vector<std::string> samples,
        std::vector<std::string>& fileNames);

//int histoProducer(int sample=0) {
int main(int argc,char* argv[]) {

    gROOT->Reset();  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);


    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    if ( argc < 2 ) {
        std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
        return 0;
    }

    // Get the python configuration
    PythonProcessDesc builder(argv[1]);
    edm::ParameterSet const& pars = 
        builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("HistoProducerParams");

    // Now get each parameter
    string eventHypo(  pars.getParameter<string>("eventHypo") );

    string inputFolder(  );
    vector<string> sampleStrings(  );

    using namespace ROOT::Math::VectorUtil;

    vector<string> fileNames;
    setInput(pars.getParameter<string>("inputFolder"),pars.getParameter<vector<string> >("input"),fileNames);
    fwlite::ChainEvent ev(fileNames);

    map<int,int> counterSkim;
    map<int,int> counterEvents;
    map<int,bool> passer;
    int tempCount = 0;
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;
    int count5 = 0;

    vector<unsigned short> bits;
    if(ev.size()){
        for( ev.toBegin(); ! ev.atEnd(); ++ev) {

            fwlite::Handle<std::vector<reco::SkimEvent> > skimEventH;
            skimEventH.getByLabel(ev,eventHypo.c_str());

            passer.clear();
            bool tempBool=false;
            bool bool1=false;
            bool bool2=false;
            bool bool3=false;
            bool bool4=false;
            bool bool5=false;
            unsigned short tempBit;
            for(vector<reco::SkimEvent>::const_iterator mySkimEvent = skimEventH.ptr()->begin();
                    mySkimEvent != skimEventH.ptr()->end(); mySkimEvent++){

                if( !(mySkimEvent->leptEtaCut(pars.getParameter<double>("etaMu"),pars.getParameter<double>("etaEl"))) ) continue;
                if( mySkimEvent->ptMin() <= pars.getParameter<double>("ptMin") ) continue;
                if( mySkimEvent->q(0)*mySkimEvent->q(1) >= 0 ) continue;
                if( mySkimEvent->isSTA(0) ) continue;
                if( mySkimEvent->isSTA(1) ) continue;
                counterSkim[1]++; passer[1]=true;
                

                if( ! mySkimEvent->hasGoodVertex() ) continue;
                counterSkim[2]++; passer[2]=true;

//                 if( fabs(mySkimEvent->d0SnT(0) ) >= pars.getParameter<double>("d0") || fabs(mySkimEvent->d0SnT(1) ) >= pars.getParameter<double>("d0") ) continue;
                if( fabs(mySkimEvent->d0Reco(0)) >= pars.getParameter<double>("d0") || fabs(mySkimEvent->d0Reco(1)) >= pars.getParameter<double>("d0") ) continue;
                counterSkim[3]++; passer[3]=true;

//                 if( fabs(mySkimEvent->dZSnT(0)) >= pars.getParameter<double>("dZ") || fabs(mySkimEvent->dZSnT(1)) >= pars.getParameter<double>("dZ") ) continue;
                if( fabs(mySkimEvent->dZReco(0)) >= pars.getParameter<double>("dZ") || fabs(mySkimEvent->dZReco(1)) >= pars.getParameter<double>("dZ") ) continue;
                counterSkim[4]++; passer[4]=true;

                if( abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= pars.getParameter<double>("isoEl") ) continue;
                if( abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= pars.getParameter<double>("isoEl") ) continue;
                if( abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= pars.getParameter<double>("isoMu") ) continue;
                if( abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= pars.getParameter<double>("isoMu") ) continue;
                counterSkim[5]++; passer[5]=true;

                if( !(mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1)) ) continue;
                counterSkim[6]++; passer[6]=true;

                if( !(mySkimEvent->passesConversion(0) && mySkimEvent->passesConversion(1)) ) continue;
                counterSkim[7]++; passer[7]=true;
                tempBool=true;

                if(!(mySkimEvent->tcMet()>pars.getParameter<double>("met"))) continue; 
                counterSkim[8]++; passer[8]=true;

                if(!(mySkimEvent->mll()>pars.getParameter<double>("mll")) ) continue;
                counterSkim[9]++; passer[9]=true;

                if(!( fabs(mySkimEvent->mll()-91.1876)>pars.getParameter<double>("mZ")) ) continue;
                counterSkim[10]++; passer[10]=true;

                if( mySkimEvent->projTcMet() <= pars.getParameter<double>("pMet") ) continue;
                counterSkim[11]++; passer[11]=true;

                if(!(mySkimEvent->nCentralJets(pars.getParameter<double>("jetPt"),pars.getParameter<double>("jetEta"))<=pars.getParameter<int>("nCentralJet")) ) continue;
                counterSkim[12]++; passer[12]=true;

                if(mySkimEvent->nSoftMu(3.)>pars.getParameter<int>("nSoftMu") ) continue;
                counterSkim[13]++; passer[13]=true;

                if(mySkimEvent->nExtraLep()>pars.getParameter<int>("nExtraLep") ) continue;
                counterSkim[14]++; passer[14]=true;

                if(mySkimEvent->bTaggedJetsUnder(pars.getParameter<double>("jetPt"),pars.getParameter<double>("bValue"))>pars.getParameter<int>("nBtagJets") ) continue;
                counterSkim[15]++; passer[15]=true;

            } //end of loop over SkimEvent

            if(tempBool) tempCount++;
            if(bool1) count1++;
            if(bool2) count2++;
            if(bool3) count3++;
            if(bool4) count4++;
            if(bool5) count5++;

            std::map<int,bool>::const_iterator itp = passer.begin();
            for(;itp!=passer.end();++itp) {
                if(itp->second) {
                    counterEvents[itp->first]++;
                } 
            }

        } // end loop over edm::events 
    }


    cout << "=========== Skim Counts ===========" << endl;
    std::map<int,int>::const_iterator its = counterSkim.begin();
    for(;its!=counterSkim.end();++its) {
        cout << setw(10) << its->first << setw(25) << its->second << endl;
    }

    cout << "========== Events Counts ==========" << endl;
    std::map<int,int>::const_iterator ite = counterEvents.begin();
    for(;ite!=counterEvents.end();++ite) {
        cout << setw(10) << ite->first << setw(25) << ite->second << endl;
    }

    cout << count2 << endl;
    cout << count3 << endl;
    cout << count4 << endl;
    cout << count5 << endl;
    cout << count1 << endl;
    cout << tempCount << endl;

    return 0;
}


void setInput(std::string inputFolder,std::vector<std::string> samples,
        std::vector<std::string>& fileNames){  
    for(unsigned int i=0;i<samples.size(); ++i){ 
        std::string tmpString=inputFolder+samples[i];
        tmpString+=".root";
        cout << "input: " << tmpString << endl;    
        fileNames.push_back(tmpString.c_str());
    }
}


