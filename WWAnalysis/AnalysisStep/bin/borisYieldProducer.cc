#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include <WWAnalysis/DataFormats/interface/SkimEvent.h>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

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

void setTriggerBits(const edm::EventBase& event,
		    reco::SkimEvent::primaryDatasetType pdType,
		    bool& valueDoubleMu, bool& valueDoubleEl, 
		    bool& valueElMu, bool& valueSingleMu);

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

  //int nMu(  fwliteParameters.getParameter<int>("nMu") );
  //int nEl(  fwliteParameters.getParameter<int>("nEl") );
  //double ptMinHighPre(  fwliteParameters.getParameter<double>("ptMinHighPre") );
  //double ptMinLowPre(  fwliteParameters.getParameter<double>("ptMinLowPre") );
  //double minProjMetPre(  fwliteParameters.getParameter<double>("minProjMetPre") );
  double minMet(  fwliteParameters.getParameter<double>("minMet") );
  double minMll(  fwliteParameters.getParameter<double>("minMll") );
  double minDiffMz(  fwliteParameters.getParameter<double>("minDiffMz") );
  //bool vetoSoftMuons(  fwliteParameters.getParameter<bool>("vetoSoftMuons") );
  int nCentralJet(  fwliteParameters.getParameter<int>("nCentralJet") );
  double jetVetoEt(  fwliteParameters.getParameter<double>("jetVetoEt") );
  double jetVetoEta(  fwliteParameters.getParameter<double>("jetVetoEta") );

  double maxMll(  fwliteParameters.getParameter<double>("maxMll") );
  double ptMinHigh(  fwliteParameters.getParameter<double>("ptMinHigh") );
  double ptMinLow(  fwliteParameters.getParameter<double>("ptMinLow") );
  double minProjMet(  fwliteParameters.getParameter<double>("minProjMet") );
  //double minRatioMetPtLL(  fwliteParameters.getParameter<double>("minRatioMetPtLL") );
  double maxDPhiLLDegrees(  fwliteParameters.getParameter<double>("maxDPhiLLDegrees") );

  string inputFolder( fwliteParameters.getParameter<string>("inputFolder") );
  vector<string> sampleStrings( fwliteParameters.getParameter<vector<string> >("input") );
  //string fileOutName(  fwliteParameters.getParameter<string>("fileOutName") );
  int maxEvents(  fwliteParameters.getParameter<int>("maxEvents") );


  //bool skipDuplicates(  fwliteParameters.getParameter<bool>("skipDuplicates") );
  bool printEventNumber(  fwliteParameters.getParameter<bool>("printEventNumber") );

  reco::SkimEvent::primaryDatasetType pdType=reco::SkimEvent::MC; //FIXME: make this configurable

  using namespace ROOT::Math::VectorUtil;


  stringstream name;
  stringstream title;


  vector<string> fileNames;
  setInput(inputFolder,sampleStrings,fileNames);
  fwlite::ChainEvent ev(fileNames);
  
  const int numberCuts(20);
  int yieldCounter[numberCuts]; for(int index=0;index<numberCuts;index++) yieldCounter[index]=0;
  long int lastEventBookkeeper[numberCuts]; for(int index=0;index<numberCuts;index++) lastEventBookkeeper[index]=-1;

  int counter(0);
  if(ev.size()){
    for( ev.toBegin(); ! ev.atEnd(); ++ev) {
      counter++; 
      if(counter == maxEvents) break;
      //if(counter == 5674) break;
      
      if(counter % 5000 == 1) {   
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
	if(event !=lastEventBookkeeper[0]){ 
	  yieldCounter[0]++; 
	  lastEventBookkeeper[0]=event;
	  if(printEventNumber) cout << "cut0, " << eventHypo0.c_str() << " ,event " << event << endl;
	}

	double pmax = mySkimEvent->ptMax();
	double pmin = mySkimEvent->ptMin();
	if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
	if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	if(!(pmin>10.) ) continue;
	if(!(pmax>20.) ) continue;

	if(event !=lastEventBookkeeper[1]){ 
	  yieldCounter[1]++; 
	  lastEventBookkeeper[1]=event;
	  if(printEventNumber) cout << "cut1, " << eventHypo0.c_str() << " ,event " << event << endl;
	}
	//if(skipDuplicates) break;
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

	if(event !=lastEventBookkeeper[2]){
	  yieldCounter[2]++; 
	  lastEventBookkeeper[2]=event;
	  if(printEventNumber) cout << "cut2, " << eventHypo0.c_str() << " ,event " << event << endl;
	}
	//if(skipDuplicates) break;	
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


	if(event !=lastEventBookkeeper[3]){
	  yieldCounter[3]++; 
	  lastEventBookkeeper[3]=event;
	  if(printEventNumber) cout << "cut3, " << eventHypoISO.c_str() << " ,event " << event << endl;
	}
	//if(skipDuplicates) break;	
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


	if(event !=lastEventBookkeeper[4]){
	  yieldCounter[4]++; 
	  lastEventBookkeeper[4]=event;
	  if(printEventNumber) cout << "cut4, " << eventHypoNOCONV.c_str() << " ,event " << event << endl;
	}
	//if(skipDuplicates) break;	
      }

      for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventIP.ptr()->begin();
	  mySkimEvent != hSkimEventIP.ptr()->end(); mySkimEvent++){

	//buffering
	reco::SkimEvent::hypoType eventType = (reco::SkimEvent::hypoType) mySkimEvent->hypo();
	double met  = mySkimEvent->pfMet();
	//double pmet = mySkimEvent->projPfMet();
	double pmet = min(mySkimEvent->projPfMet(),mySkimEvent->projChargedMetSmurf());
	double pmax = mySkimEvent->ptMax();
	double pmin = mySkimEvent->ptMin();
	double mll  = mySkimEvent->mll();
	double dmZ  = fabs(mySkimEvent->mll()-91.1876);
	int nSoftMu = mySkimEvent->nSoftMu();
	int nJet = mySkimEvent->nCentralJets(jetVetoEt,jetVetoEta,true,0);
	//d8 = mySkimEvent->tcMet()/mySkimEvent->pTll();
	double dPhi = mySkimEvent->dPhillInDegrees();
	double mT = mySkimEvent->mTHiggs(reco::SkimEvent::PFMET);
	int nExtraLep = mySkimEvent->nExtraLep();
	bool passTriggerMatching = mySkimEvent->triggerMatchingCut(pdType);
	if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
	if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
	if( !(mySkimEvent->leptEtaCut()) ) continue;
	if(!(pmin>10.) ) continue;
	if(!(pmax>20.) ) continue;

	if(event !=lastEventBookkeeper[5]){
	  yieldCounter[5]++; 
	  lastEventBookkeeper[5]=event;
	  if(printEventNumber) cout << "cut5, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}


	
	/*
	  if(skimEventType_ ==2 || skimEventType_ ==3){
	bool vDoubleMu(false),vDoubleEl(false),vElMu(false),vSingleMu(false);
	edm::EventBase const & eventBase = ev;
	setTriggerBits(eventBase,pdType,vDoubleMu,vDoubleEl,vElMu,vSingleMu);
	bool triggerDecision(false);
	if(eventType == reco::SkimEvent::WWMUMU)
	  triggerDecision = vDoubleMu || vSingleMu;
	if(eventType == reco::SkimEvent::WWELEL)
	  triggerDecision = vDoubleEl;
	if(eventType == reco::SkimEvent::WWELMU || eventType == reco::SkimEvent::WWMUEL)
	  triggerDecision = vSingleMu || vElMu;
	*/
	//if(!triggerDecision) continue;
	//if(event !=lastEventBookkeeper[7]){
	//  yieldCounter[7]++; 
	//  lastEventBookkeeper[7]=event;
	  //if(printEventNumber) cout << "cut7, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	//}



	//if(!passTriggerMatching) continue;
	//if(event !=lastEventBookkeeper[8]){
	//  yieldCounter[8]++; 
	//  lastEventBookkeeper[8]=event;
	  //if(printEventNumber) cout << "cut8, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	//}


	if(!(met>minMet)) continue;
	if(event !=lastEventBookkeeper[6]){
	  yieldCounter[6]++; 
	  lastEventBookkeeper[6]=event;
	  if(printEventNumber) cout << "cut6, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(mll>minMll)) continue;
	if(event !=lastEventBookkeeper[7]){
	  yieldCounter[7]++; 
	  lastEventBookkeeper[7]=event;
	  if(printEventNumber) cout << "cut7, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(dmZ>minDiffMz)) continue;
	if(event !=lastEventBookkeeper[8]){
	  yieldCounter[8]++; 
	  lastEventBookkeeper[8]=event;
	  if(printEventNumber) cout << "cut8, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(pmet>minProjMet)) continue;
	if(event !=lastEventBookkeeper[9]){
	  yieldCounter[9]++; 
	  lastEventBookkeeper[9]=event;
	  if(printEventNumber) cout << "cut9, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}
	
	if(!(nJet==nCentralJet)) continue;
	if(event !=lastEventBookkeeper[10]){
	  yieldCounter[10]++; 
	  lastEventBookkeeper[10]=event;
	  if(printEventNumber) cout << "cut10, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(nSoftMu==0)) continue;
	if(event !=lastEventBookkeeper[11]){
	  yieldCounter[11]++; 
	  lastEventBookkeeper[11]=event;
	  if(printEventNumber) cout << "cut11, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(nExtraLep==0)) continue;
	if(event !=lastEventBookkeeper[12]){
	  yieldCounter[12]++; 
	  lastEventBookkeeper[12]=event;
	  if(printEventNumber) cout << "cut12, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	
	if(!(mySkimEvent->bTaggedJetsUnder( jetVetoEt, 2.1) == 0)) continue;
	if(event !=lastEventBookkeeper[13]){
	  yieldCounter[13]++; 
	  lastEventBookkeeper[13]=event;
	  if(printEventNumber) cout << "cut13, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(mll<maxMll)) continue;
	if(event !=lastEventBookkeeper[14]){
	  yieldCounter[14]++; 
	  lastEventBookkeeper[14]=event;
	  if(printEventNumber) cout << "cut14, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(pmax>ptMinHigh)) continue;
	if(event !=lastEventBookkeeper[15]){
	  yieldCounter[15]++; 
	  lastEventBookkeeper[15]=event;
	  if(printEventNumber) cout << "cut15, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(pmin>ptMinLow)) continue;
	if(event !=lastEventBookkeeper[16]){
	  yieldCounter[16]++; 
	  lastEventBookkeeper[16]=event;
	  if(printEventNumber) cout << "cut16, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(dPhi<maxDPhiLLDegrees) ) continue;
	if(event !=lastEventBookkeeper[17]){
	  yieldCounter[17]++; 
	  lastEventBookkeeper[17]=event;
	  if(printEventNumber) cout << "cut17, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}

	if(!(mT>90 && mT<160) ) continue;
	if(event !=lastEventBookkeeper[18]){
	  yieldCounter[18]++; 
	  lastEventBookkeeper[18]=event;
	  if(printEventNumber) cout << "cut18, " << eventHypoIP.c_str() << " ,event " << event << endl;	
	}
	


	//--- last cut of the cut-based SELECTION
      }//end of loop over SkimEvent
      
    }// end loop over edm::events 
  }
  
  cout << endl;
  cout << "======= SUMMARY ======" << endl;
  cout << "# analyzed events: " << counter << endl << endl;
  cout << "yields: " << endl;
  for(int index=0;index<numberCuts;index++){
    cout << "cut" << index << ": " <<  yieldCounter[index] << endl;
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


void setTriggerBits(const edm::EventBase& event,
		    reco::SkimEvent::primaryDatasetType pdType,
		    bool& valueDoubleMu, bool& valueDoubleEl, 
		    bool& valueElMu, bool& valueSingleMu)
{
  valueDoubleMu=false;
  valueDoubleEl=false;
  valueElMu=false;
  valueSingleMu=false;


  edm::Handle< edm::TriggerResults > hltTriggerResults;
  if(pdType==reco::SkimEvent::MC)
    event.getByLabel( edm::InputTag("TriggerResults","","REDIGI311X"), hltTriggerResults);
  else
    event.getByLabel( edm::InputTag("TriggerResults","","HLT"), hltTriggerResults);

  edm::TriggerNames triggerNames = event.triggerNames(*hltTriggerResults );

  if(pdType==reco::SkimEvent::MC){
    unsigned int bitSingleMu      = triggerNames.triggerIndex("HLT_Mu21_v1");
    unsigned int bitDoubleMu      = triggerNames.triggerIndex("HLT_DoubleMu5_v1");
    unsigned int bitDoubleEl      = triggerNames.triggerIndex("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2");
    unsigned int bitCrossElMu     = triggerNames.triggerIndex("HLT_Mu5_Ele17_v2");
    unsigned int bitCrossMuEl     = triggerNames.triggerIndex("HLT_Mu11_Ele8_v1");

    unsigned int sizeNames=triggerNames.size();
    if(bitSingleMu != sizeNames) valueSingleMu = hltTriggerResults->accept(bitSingleMu);
    if(bitDoubleMu != sizeNames) valueDoubleMu = hltTriggerResults->accept(bitDoubleMu);
    if(bitDoubleEl != sizeNames) valueDoubleEl = hltTriggerResults->accept(bitDoubleEl);

    if(bitCrossElMu != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitCrossElMu);
    if(bitCrossMuEl != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitCrossMuEl);



  }

  /*
  //unsigned int bitDoubleMu  = triggerNames.triggerIndex("HLT_DoubleMu7_v*");
  //unsigned int bitDoubleEl  = 
  triggerNames.triggerIndex("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*");
  unsigned int bitMuEl1      = triggerNames.triggerIndex("HLT_Mu17_Ele8_CaloIdL_v1");
  unsigned int bitMuEl2      = triggerNames.triggerIndex("HLT_Mu17_Ele8_CaloIdL_v2");
  unsigned int bitMuEl3      = triggerNames.triggerIndex("HLT_Mu17_Ele8_CaloIdL_v3");

  unsigned int bitElMu1      = triggerNames.triggerIndex("HLT_Mu8_Ele17_CaloIdL_v1");
  unsigned int bitElMu2      = triggerNames.triggerIndex("HLT_Mu8_Ele17_CaloIdL_v2");
  unsigned int bitElMu3      = triggerNames.triggerIndex("HLT_Mu8_Ele17_CaloIdL_v3");

  unsigned int sizeNames=triggerNames.size();
  if(bitMuEl1 != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitMuEl1);
  if(bitMuEl2 != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitMuEl2);
  if(bitMuEl3 != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitMuEl3);

  if(bitElMu1 != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitElMu1);
  if(bitElMu2 != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitElMu2);
  if(bitElMu3 != sizeNames) valueElMu = valueElMu || hltTriggerResults->accept(bitElMu3);
  //unsigned int bitSingleMu = triggerNames.triggerIndex("HLT_Mu24_v*");
	
  //valueDoubleMu = hltTriggerResults->accept(bitDoubleMu);
  //valueDoubleEl = hltTriggerResults->accept(bitDoubleEl);
  //valueElMu = (hltTriggerResults->accept(bitElMu1) || hltTriggerResults->accept(bitElMu2));
  //valueSingleMu = hltTriggerResults->accept(bitSingleMu);
  */
}
