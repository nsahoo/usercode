#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

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



void setInput(std::string inputFolder,std::vector<std::string> samples,
	      std::vector<std::string>& fileNames);


int main(int argc,char* argv[]) {
  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 


  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  // Get the python configuration
  PythonProcessDesc builder(argv[1]);
  edm::ParameterSet const& fwliteParameters = 
    builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("HistoProducerParams");
  
  // Now get each parameter
  using namespace std;
  string eventHypo0(  fwliteParameters.getParameter<string>("eventHypo0") );
  string eventHypoID(  fwliteParameters.getParameter<string>("eventHypoID") );
  string eventHypoISO(  fwliteParameters.getParameter<string>("eventHypoISO") );
  string eventHypoNOCONV(  fwliteParameters.getParameter<string>("eventHypoNOCONV") );
  string eventHypoIP(  fwliteParameters.getParameter<string>("eventHypoIP") );

  int nMu(  fwliteParameters.getParameter<int>("nMu") );
  int nEl(  fwliteParameters.getParameter<int>("nEl") );
  bool useNewExpectHits(  fwliteParameters.getParameter<bool>("useNewExpectHits") );
  double ptMinHighPre(  fwliteParameters.getParameter<double>("ptMinHighPre") );
  double ptMinLowPre(  fwliteParameters.getParameter<double>("ptMinLowPre") );
  double minProjMetPre(  fwliteParameters.getParameter<double>("minProjMetPre") );
  double minMet(  fwliteParameters.getParameter<double>("minMet") );
  double minMll(  fwliteParameters.getParameter<double>("minMll") );
  double minDiffMz(  fwliteParameters.getParameter<double>("minDiffMz") );
  bool vetoSoftMuons(  fwliteParameters.getParameter<bool>("vetoSoftMuons") );
  int nCentralJet(  fwliteParameters.getParameter<int>("nCentralJet") );
  double jetVetoEt(  fwliteParameters.getParameter<double>("jetVetoEt") );
  double jetVetoEta(  fwliteParameters.getParameter<double>("jetVetoEta") );

  double maxMll(  fwliteParameters.getParameter<double>("maxMll") );
  double ptMinHigh(  fwliteParameters.getParameter<double>("ptMinHigh") );
  double ptMinLow(  fwliteParameters.getParameter<double>("ptMinLow") );
  double minProjMet(  fwliteParameters.getParameter<double>("minProjMet") );
  double minRatioMetPtLL(  fwliteParameters.getParameter<double>("minRatioMetPtLL") );
  double maxDPhiLL(  fwliteParameters.getParameter<double>("maxDPhiLL") );

  string inputFolder( fwliteParameters.getParameter<string>("inputFolder") );
  vector<string> sampleStrings( fwliteParameters.getParameter<vector<string> >("input") );
  string fileOutName(  fwliteParameters.getParameter<string>("fileOutName") );
  int maxEvents(  fwliteParameters.getParameter<int>("maxEvents") );

  using namespace ROOT::Math::VectorUtil;


  stringstream name;
  stringstream title;


  vector<string> fileNames;
  setInput(inputFolder,sampleStrings,fileNames);
  fwlite::ChainEvent ev(fileNames);


  int counter(0);
  if(ev.size()){
    for( ev.toBegin(); ! ev.atEnd(); ++ev) {
      counter++; 
      if(counter == maxEvents) break;
      //if(counter == 20) break;
      if(counter % 1000 == 1) {   
	cout << "counter: " << counter << endl;
	//cout << "file name: " << ev.getTFile()->GetName() << endl;
      }
      
      fwlite::Handle<std::vector<reco::SkimEvent> > hSkimEvent0;
      fwlite::Handle<std::vector<reco::SkimEvent> > hSkimEventID;
      fwlite::Handle<std::vector<reco::SkimEvent> > hSkimEventISO;
      fwlite::Handle<std::vector<reco::SkimEvent> > hSkimEventNOCONV;
      fwlite::Handle<std::vector<reco::SkimEvent> > hSkimEventIP;

      hSkimEvent0.getByLabel(ev,eventHypo0.c_str());
      hSkimEventID.getByLabel(ev,eventHypoID.c_str());
      hSkimEventISO.getByLabel(ev,eventHypoISO.c_str());
      hSkimEventNOCONV.getByLabel(ev,eventHypoNOCONV.c_str());
      hSkimEventIP.getByLabel(ev,eventHypoIP.c_str());

      long int event = ev.eventAuxiliary().event();

      for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEvent0.ptr()->begin();
	  mySkimEvent != hSkimEvent0.ptr()->end(); mySkimEvent++){
	cout << "cut0, " << eventHypo0.c_str() << " ,event " << event << endl;

	double pmax = mySkimEvent->ptMax();
	double pmin = mySkimEvent->ptMin();
	if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
	if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	if(!(pmin>10.) ) continue;
	if(!(pmax>20.) ) continue;
	cout << "cut1, " << eventHypo0.c_str() << " ,event " << event << endl;		
      }

      
      for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventID.ptr()->begin();
	  mySkimEvent != hSkimEventID.ptr()->end(); mySkimEvent++){
	double pmax = mySkimEvent->ptMax();
	double pmin = mySkimEvent->ptMin();
	if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
	if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	if(!(pmin>10.) ) continue;
	if(!(pmax>20.) ) continue;
	cout << "cut2, " << eventHypo0.c_str() << " ,event " << event << endl;	
      }
      

      for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventISO.ptr()->begin();
	  mySkimEvent != hSkimEventISO.ptr()->end(); mySkimEvent++){
	double pmax = mySkimEvent->ptMax();
	double pmin = mySkimEvent->ptMin();
	if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
	if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	if(!(pmin>10.) ) continue;
	if(!(pmax>20.) ) continue;

	cout << "cut3, " << eventHypoISO.c_str() << " ,event " << event << endl;	
      }

      for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventNOCONV.ptr()->begin();
	  mySkimEvent != hSkimEventNOCONV.ptr()->end(); mySkimEvent++){

	double pmax = mySkimEvent->ptMax();
	double pmin = mySkimEvent->ptMin();
	if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
	if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	if(!(pmin>10.) ) continue;
	if(!(pmax>20.) ) continue;

	cout << "cut4, " << eventHypoNOCONV.c_str() << " ,event " << event << endl;	
      }

      for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventIP.ptr()->begin();
	  mySkimEvent != hSkimEventIP.ptr()->end(); mySkimEvent++){

	//buffering
	double met  = mySkimEvent->tcMet();
	double pmet = mySkimEvent->projTcMet();
	double pmax = mySkimEvent->ptMax();
	double pmin = mySkimEvent->ptMin();
	double mll  = mySkimEvent->mll();
	double dmZ  = fabs(mySkimEvent->mll()-91.1876);
	int nSoftMu = mySkimEvent->nSoftMu();
	int nJet = mySkimEvent->nCentralJets(jetVetoEt,jetVetoEta);
	//d8 = mySkimEvent->tcMet()/mySkimEvent->pTll();
	double dPhi = mySkimEvent->dPhill();
	int nExtraLep = mySkimEvent->nExtraLep();

	if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
	if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	if(!(pmin>10.) ) continue;
	if(!(pmax>20.) ) continue;
	cout << "cut5, " << eventHypoIP.c_str() << " ,event " << event << endl;	


	if(!(nExtraLep==0)) continue;
	cout << "cut6, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(met>minMet)) continue;
	cout << "cut7, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(mll>minMll)) continue;
	cout << "cut8, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(dmZ>minDiffMz)) continue;
	cout << "cut9, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(pmet>minProjMet)) continue;
	cout << "cut10, " << eventHypoIP.c_str() << " ,event " << event << endl;
	
	if(!(nJet==0)) continue;
	cout << "cut11, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(nSoftMu==0)) continue;
	cout << "cut12, " << eventHypoIP.c_str() << " ,event " << event << endl;

	
	if(!(mySkimEvent->bTaggedJetsUnder( jetVetoEt, 2.1) == 0)) continue;
	cout << "cut13, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(mll<50)) continue;
	cout << "cut14, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(pmax>30)) continue;
	cout << "cut15, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(pmin>25)) continue;
	cout << "cut16, " << eventHypoIP.c_str() << " ,event " << event << endl;

	if(!(dPhi<(60./180*M_PI))) continue;
	cout << "cut17, " << eventHypoIP.c_str() << " ,event " << event << endl;
	


	//--- last cut of the cut-based SELECTION
      }//end of loop over SkimEvent
      
    }// end loop over edm::events 
  }
  
  return 0;
}


 

void setInput(std::string inputFolder,std::vector<std::string> samples,
	      std::vector<std::string>& fileNames){  
  using namespace std;
  for(unsigned int i=0;i<samples.size(); ++i){ 
    std::string tmpString=inputFolder+samples[i]; 
    cout << "input: " << tmpString << endl;    
    fileNames.push_back(tmpString.c_str());
  }
  
}


