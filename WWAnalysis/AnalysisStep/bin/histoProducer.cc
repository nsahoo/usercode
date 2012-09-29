#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
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
//using namespace reco;


void fillHistos(TH1F* met, TH1F* projMet, TH1F* ptMax, TH1F* ptMin, TH1F* mll,		
		TH1F* dMllMz,  TH1F* nSoftMu,  TH1F* nJet,  TH1F* metOverPTll,
		TH1F* dPhill,
		double d0, double d1,double d2, double d3, double d4,double d5,
		int i6,int i7, double d8, double d9);
//std::vector<reco::SkimEvent>::const_iterator event);


void setInput(std::string inputFolder,std::vector<std::string> samples,
	      std::vector<std::string>& fileNames);


//int histoProducer(int sample=0) {
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
  string eventHypo(  fwliteParameters.getParameter<string>("eventHypo") );

  int nMu(  fwliteParameters.getParameter<int>("nMu") );
  int nEl(  fwliteParameters.getParameter<int>("nEl") );
  bool useNewExpectHits(  fwliteParameters.getParameter<bool>("useNewExpectHits") );
  double ptMinHighPre(  fwliteParameters.getParameter<double>("ptMinHighPre") );
  double ptMinLowPre(  fwliteParameters.getParameter<double>("ptMinLowPre") );
  double minProjMetPre(  fwliteParameters.getParameter<double>("minProjMetPre") );
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

  //TCanvas* myC = new TCanvas("myCanvas","myCanvas",600,600); myC->cd(1)
  const unsigned int nSelections(14);
  TH1F* met[nSelections];
  TH1F* projMet[nSelections];
  TH1F* ptmax[nSelections];
  TH1F* ptmin[nSelections];
  TH1F* mll[nSelections];
  TH1F* dMllMz[nSelections];
  TH1F* nSoftMu[nSelections];
  TH1F* nJet[nSelections];
  TH1F* metOverPTll[nSelections];
  TH1F* dPhill[nSelections];


  stringstream name;
  stringstream title;

  for(unsigned int index=0; index<nSelections; ++index){
    name.str(string()); title.str(string());
    name << "met" << index;
    title << "tcMET (GeV)";
    met[index] = new TH1F(name.str().c_str(),title.str().c_str(),50, 0., 200.);

    name.str(string()); title.str(string());
    name << "projMet" << index;
    title << "projected tcMET (GeV)";
    projMet[index] = new TH1F(name.str().c_str(),title.str().c_str(),50, 0., 200.);

    name.str(string()); title.str(string());
    name << "ptMax" << index;
    title << "lepton pt_{max} (GeV/c)";
    ptmax[index] = new TH1F(name.str().c_str(),title.str().c_str(),50, 0., 150.);

    name.str(string()); title.str(string());
    name << "ptMin" << index;
    title << "lepton pt_{min} (GeV/c)";
    ptmin[index] = new TH1F(name.str().c_str(),title.str().c_str(),50, 0., 150.);

    name.str(string()); title.str(string());
    name << "mll" << index;
    title << "mass(l,l) (GeV)";
    mll[index] = new TH1F(name.str().c_str(),title.str().c_str(),50, 0., 150.);

    name.str(string()); title.str(string());
    name << "dMllMz" << index;
    title << "| mass(l,l) - m_Z | (GeV)";
    dMllMz[index] = new TH1F(name.str().c_str(),title.str().c_str(),50, 0., 150.);

    name.str(string()); title.str(string());
    name << "nSoftMu" << index;
    title << "# extra soft muons";
    nSoftMu[index] = new TH1F(name.str().c_str(),title.str().c_str(),5, -1., 4.);

    name.str(string()); title.str(string());
    name << "nJet" << index;
    title << "# central Jets";
    nJet[index] = new TH1F(name.str().c_str(),title.str().c_str(),5, 0., 5.);

    name.str(string()); title.str(string());
    name << "metOverPTll" << index;
    title << "met / pt_{ll}";
    metOverPTll[index] = new TH1F(name.str().c_str(),title.str().c_str(),30, 0., 5.);

    name.str(string()); title.str(string());
    name << "dPhill" << index;
    title << "DeltaPhi_{ll}";
    dPhill[index] = new TH1F(name.str().c_str(),title.str().c_str(),30, 0., M_PI);    
  }
  

  vector<string> fileNames;
  setInput(inputFolder,sampleStrings,fileNames);
  fwlite::ChainEvent ev(fileNames);

  double d0,d1,d2,d3,d4,d5;
  int i6,i7,i10;
  double d8,d9;

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
      
      fwlite::Handle<std::vector<reco::SkimEvent> > skimEventH;
      skimEventH.getByLabel(ev,eventHypo.c_str());

      if(skimEventH->size()==0) continue;
      
      for(vector<reco::SkimEvent>::const_iterator mySkimEvent = skimEventH.ptr()->begin();
	  mySkimEvent != skimEventH.ptr()->end(); mySkimEvent++){
	//cout << "Beginning, event: " << ev.eventAuxiliary().event() << endl;


      
	//if(!(mySkimEvent->nMu()==nMu && mySkimEvent->nEl()==nEl)) continue;
	
	//buffering
	d0 = mySkimEvent->tcMet();
	d1 = mySkimEvent->projTcMet();
	d2 = mySkimEvent->ptMax();
	d3 = mySkimEvent->ptMin();
	d4 = mySkimEvent->mll();
	d5 = fabs(mySkimEvent->mll()-91.1876);
	i6 = mySkimEvent->nSoftMu();
	i7 = mySkimEvent->nCentralJets(jetVetoEt,jetVetoEta);
	d8 = mySkimEvent->tcMet()/mySkimEvent->pTll();
	d9 = mySkimEvent->dPhill();
	i10 = mySkimEvent->nExtraLep();
	
	// ----//as it is after SkimEvent producer
	
	// requires mode ELEL,ELMU or MUMU
	if(!(mySkimEvent->q(0)*mySkimEvent->q(1) == -1) ) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	//if( !(mySkimEvent->eleExpHitCut(useNewExpectHits)) ) continue;
	fillHistos(met[1],projMet[1],ptmax[1],ptmin[1],mll[1],
		   dMllMz[1],nSoftMu[1],nJet[1],metOverPTll[1],dPhill[1],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);
	//mySkimEvent);
	
      
	//if( !mySkimEvent->leptEtaCut()) {
	//	cout << "here we are!" << endl;
	//cout << "eta(0),eta(1): " << mySkimEvent->eta(0) 
	//   << " , " << mySkimEvent->eta(1) << endl;
	//}
	
	
	// pre-sel on pt
	//if(!(d2>=ptMinHighPre && d3>ptMinLowPre) ) continue;
	if(!(d2>=20.0 && d3>20.0) ) continue;
	fillHistos(met[2],projMet[2],ptmax[2],ptmin[2],mll[2],
		   dMllMz[2],nSoftMu[2],nJet[2],metOverPTll[2],dPhill[2],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

	
	
	// pre-sel on projTcMet
	if(!(d0>20)) continue; //cut on standard met
	//if(!(d1>minProjMetPre)) continue; //cut on projected met
	fillHistos(met[3],projMet[3],ptmax[3],ptmin[3],mll[3],
		 dMllMz[3],nSoftMu[3],nJet[3],metOverPTll[3],dPhill[3],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

	

	// pre-sel on mll
	if(!(d4>minMll) ) continue;
	fillHistos(met[4],projMet[4],ptmax[4],ptmin[4],mll[4],
		   dMllMz[4],nSoftMu[4],nJet[4],metOverPTll[4],dPhill[4],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

	
	// pre-sel on mll-mz
	if(!( d5>minDiffMz) ) continue;
	fillHistos(met[5],projMet[5],ptmax[5],ptmin[5],mll[5],
		   dMllMz[5],nSoftMu[5],nJet[5],metOverPTll[5],dPhill[5],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

	
	
	// selection cut on projTcMet+met
	if(!(d1 >20 && d0> minProjMet ) ) continue;
	fillHistos(met[6],projMet[6],ptmax[6],ptmin[6],mll[6],
		   dMllMz[6],nSoftMu[6],nJet[6],metOverPTll[6],dPhill[6],
		 d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

	
	
	// pre-sel on central jet
	if(!(i7<=nCentralJet) ) continue;
	fillHistos(met[7],projMet[7],ptmax[7],ptmin[7],mll[7],
		   dMllMz[7],nSoftMu[7],nJet[7],metOverPTll[7],dPhill[7],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

	
	
	// pre-sel softMuon veto
	if(vetoSoftMuons && i6>0 ) continue;
	fillHistos(met[8],projMet[8],ptmax[8],ptmin[8],mll[8],
		   dMllMz[8],nSoftMu[8],nJet[8],metOverPTll[8],dPhill[8],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

	// extra lepton veto
	if(i10>0 ) continue;
	fillHistos(met[9],projMet[9],ptmax[9],ptmin[9],mll[9],
		   dMllMz[9],nSoftMu[9],nJet[9],metOverPTll[9],dPhill[9],
		   d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);

//     if( abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) > 0.15 ) continue;
//     if( abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) > 0.15 ) continue;
//     if( abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) > 0.10 ) continue;
//     if( abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) > 0.10 ) continue;
	cout << "event passing final selection: " << ev.eventAuxiliary().event() << endl;


	//--- last cut of the cut-based SELECTION
      }//end of loop over SkimEvent
      
    }// end loop over edm::events 
  }


  // write histograms on disk
  TFile output(fileOutName.c_str(),"RECREATE");  
  //projPfMet0->Draw(); gPad->Update();
  
  for(unsigned int index=0; index<14; ++index){
    met[index]->Write();
    projMet[index]->Write();
    ptmax[index]->Write();
    ptmin[index]->Write();
    mll[index]->Write();
    dMllMz[index]->Write();
    nSoftMu[index]->Write();
    nJet[index]->Write();
    metOverPTll[index]->Write();
    dPhill[index]->Write();    
  }

  return 0;
}


 

void setInput(std::string inputFolder,std::vector<std::string> samples,
	      std::vector<std::string>& fileNames){  
  using namespace std;
  for(unsigned int i=0;i<samples.size(); ++i){ 
    std::string tmpString=inputFolder+samples[i];
    tmpString+=".root";
    cout << "input: " << tmpString << endl;    
    fileNames.push_back(tmpString.c_str());
  }
}


void fillHistos(TH1F* met, TH1F* projMet, TH1F* ptMax, TH1F* ptMin, TH1F* mll,
		TH1F* dMllMz,  TH1F* nSoftMu,  TH1F* nJet,  TH1F* metOverPTll,TH1F* dPhill,
		double d0,double d1,double d2, double d3, double d4,double d5,
		int i6,int i7, double d8, double d9)
//std::vector<reco::SkimEvent>::const_iterator event){
{
  met->Fill(d0);
  projMet->Fill(d1);
  ptMax->Fill(d2);
  ptMin->Fill(d3);
  mll->Fill(d4);
  dMllMz->Fill(d5);
  nSoftMu->Fill(i6);
  nJet->Fill(i7);
  metOverPTll->Fill(d8);
  dPhill->Fill(d9);

  /*
  projMet->Fill(event->projTcMet(0,1));
  ptMax->Fill(event->pt(0));
  ptMin->Fill(event->pt(1));
  mll->Fill(event->mll());
  dMllMz->Fill(fabs(event->mll()-91.188));
  nSoftMu->Fill(event->nExtraMu());
  nJet->Fill(event->nCentralJets(20,3.0));
  metOverPTll->Fill(event->tcMet()/event->pTll());
  dPhill->Fill(event->dPhill());
  */
}

