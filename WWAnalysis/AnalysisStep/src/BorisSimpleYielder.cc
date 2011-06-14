#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "WWAnalysis/AnalysisStep/interface/BorisSimpleYielder.h"
#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"


/// default constructor
BorisSimpleYielder::BorisSimpleYielder(const edm::ParameterSet& cfg, TFileDirectory& fs): 
  edm::BasicAnalyzer::BasicAnalyzer(cfg, fs),
  eventHypo0_(cfg.getParameter<edm::InputTag>("eventHypo0")),
  eventHypoID_(cfg.getParameter<edm::InputTag>("eventHypoID")),
  eventHypoISO_(cfg.getParameter<edm::InputTag>("eventHypoISO")),
  eventHypoNOCONV_(cfg.getParameter<edm::InputTag>("eventHypoNOCONV")),
  eventHypoIP_(cfg.getParameter<edm::InputTag>("eventHypoIP")),
  minMet(  cfg.getParameter<double>("minMet") ),
  minMll(  cfg.getParameter<double>("minMll") ),
  minDiffMz(  cfg.getParameter<double>("minDiffMz") ),
  //bool vetoSoftMuons(  cfg.getParameter<bool>("vetoSoftMuons") ),
  nCentralJet(  cfg.getParameter<int>("nCentralJet") ),
  jetVetoEt(  cfg.getParameter<double>("jetVetoEt") ),
  jetVetoEta(  cfg.getParameter<double>("jetVetoEta") ),
  maxMll(  cfg.getParameter<double>("maxMll") ),
  ptMinHigh(  cfg.getParameter<double>("ptMinHigh") ),
  ptMinLow(  cfg.getParameter<double>("ptMinLow") ),
  minProjMet(  cfg.getParameter<double>("minProjMet") ),
  //minRatioMetPtLL(  cfg.getParameter<double>("minRatioMetPtLL") ),
  maxDPhiLLDegrees(  cfg.getParameter<double>("maxDPhiLLDegrees") ),
  printEventNumber(  cfg.getParameter<bool>("printEventNumber") )
{
  for(int index=0;index<numberCuts;index++) yieldCounter[index]=0;
  for(int index=0;index<numberCuts;index++) lastEventBookkeeper[index]=-1;
  //  hists_["dphi"  ] = fs.make<TH1F>("dphi"  , "dphi"  ,  100,  0., M_PI);
}

/// everything that needs to be done during the event loop
void 
BorisSimpleYielder::analyze(const edm::EventBase& ev)
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  //using pat::Muon;
  using namespace std;



  edm::Handle<std::vector<reco::SkimEvent> > hSkimEvent0;
  edm::Handle<std::vector<reco::SkimEvent> > hSkimEventID;
  edm::Handle<std::vector<reco::SkimEvent> > hSkimEventISO;
  edm::Handle<std::vector<reco::SkimEvent> > hSkimEventNOCONV;
  edm::Handle<std::vector<reco::SkimEvent> > hSkimEventIP;
  
  ev.getByLabel(eventHypo0_,hSkimEvent0);
  ev.getByLabel(eventHypoID_,hSkimEventID);
  ev.getByLabel(eventHypoISO_,hSkimEventISO);
  ev.getByLabel(eventHypoNOCONV_,hSkimEventNOCONV);
  ev.getByLabel(eventHypoIP_,hSkimEventIP);
  
  long int event = ev.eventAuxiliary().event();

  for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEvent0->begin();
      mySkimEvent != hSkimEvent0->end(); mySkimEvent++){
    if(event !=lastEventBookkeeper[0]){ 
      yieldCounter[0]++; 
      lastEventBookkeeper[0]=event;
      if(printEventNumber) cout << "cut0 ,event " << event << endl;
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
      if(printEventNumber) cout << "cut1 ,event " << event << endl;
    }
    //if(skipDuplicates) break;
  }

      
  for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventID->begin();
      mySkimEvent != hSkimEventID->end(); mySkimEvent++){
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
      if(printEventNumber) cout << "cut2 ,event " << event << endl;
    }
    //if(skipDuplicates) break;	
  }
      

  for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventISO->begin();
      mySkimEvent != hSkimEventISO->end(); mySkimEvent++){
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
      if(printEventNumber) cout << "cut3 ,event " << event << endl;
    }
    //if(skipDuplicates) break;	
  }

  for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventNOCONV->begin();
      mySkimEvent != hSkimEventNOCONV->end(); mySkimEvent++){
    
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
      if(printEventNumber) cout << "cut4 ,event " << event << endl;
    }
    //if(skipDuplicates) break;	
  }

  for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hSkimEventIP->begin();
      mySkimEvent != hSkimEventIP->end(); mySkimEvent++){
    
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
    //bool passTriggerMatching = mySkimEvent->triggerMatchingCut(pdType);
    if( mySkimEvent->q(0)*mySkimEvent->q(1) == 1 ) continue;
    if( mySkimEvent->isSTA(0) || mySkimEvent->isSTA(1)) continue;
    if( !(mySkimEvent->leptEtaCut()) ) continue;
    if(!(pmin>10.) ) continue;
    if(!(pmax>20.) ) continue;
    
    if(event !=lastEventBookkeeper[5]){
      yieldCounter[5]++; 
      lastEventBookkeeper[5]=event;
      if(printEventNumber) cout << "cut5 ,event " << event << endl;
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
      if(printEventNumber) cout << "cut6 ,event " << event << endl;
    }
    
    if(!(mll>minMll)) continue;
    if(event !=lastEventBookkeeper[7]){
      yieldCounter[7]++; 
      lastEventBookkeeper[7]=event;
      if(printEventNumber) cout << "cut7 ,event " << event << endl;
    }
    
    if(!(dmZ>minDiffMz)) continue;
    if(event !=lastEventBookkeeper[8]){
      yieldCounter[8]++; 
	  lastEventBookkeeper[8]=event;
      if(printEventNumber) cout << "cut8 ,event " << event << endl;
    }
    
    if(!(pmet>minProjMet)) continue;
    if(event !=lastEventBookkeeper[9]){
      yieldCounter[9]++; 
      lastEventBookkeeper[9]=event;
      if(printEventNumber) cout << "cut9 ,event " << event << endl;
    }
	
    if(!(nJet==nCentralJet)) continue;
    if(event !=lastEventBookkeeper[10]){
      yieldCounter[10]++; 
      lastEventBookkeeper[10]=event;
      if(printEventNumber) cout << "cut10 ,event " << event << endl;
    }

    if(!(nSoftMu==0)) continue;
    if(event !=lastEventBookkeeper[11]){
      yieldCounter[11]++; 
      lastEventBookkeeper[11]=event;
      if(printEventNumber) cout << "cut11 ,event " << event << endl;
    }

    if(!(nExtraLep==0)) continue;
    if(event !=lastEventBookkeeper[12]){
      yieldCounter[12]++; 
      lastEventBookkeeper[12]=event;
      if(printEventNumber) cout << "cut12 ,event " << event << endl;
    }

    
    if(!(mySkimEvent->bTaggedJetsUnder( jetVetoEt, 2.1) == 0)) continue;
    if(event !=lastEventBookkeeper[13]){
      yieldCounter[13]++; 
      lastEventBookkeeper[13]=event;
      if(printEventNumber) cout << "cut13 ,event " << event << endl;
    }

    if(!(mll<maxMll)) continue;
    if(event !=lastEventBookkeeper[14]){
      yieldCounter[14]++; 
      lastEventBookkeeper[14]=event;
      if(printEventNumber) cout << "cut14 ,event " << event << endl;
    }

    if(!(pmax>ptMinHigh)) continue;
    if(event !=lastEventBookkeeper[15]){
      yieldCounter[15]++; 
      lastEventBookkeeper[15]=event;
      if(printEventNumber) cout << "cut15 ,event " << event << endl;
    }

    if(!(pmin>ptMinLow)) continue;
    if(event !=lastEventBookkeeper[16]){
      yieldCounter[16]++; 
      lastEventBookkeeper[16]=event;
      if(printEventNumber) cout << "cut16 ,event " << event << endl;
    }

    if(!(dPhi<maxDPhiLLDegrees) ) continue;
    if(event !=lastEventBookkeeper[17]){
      yieldCounter[17]++; 
      lastEventBookkeeper[17]=event;
      if(printEventNumber) cout << "cut17 ,event " << event << endl;
    }

    if(!(mT>90 && mT<160) ) continue;
    if(event !=lastEventBookkeeper[18]){
      yieldCounter[18]++; 
      lastEventBookkeeper[18]=event;
      if(printEventNumber) cout << "cut18 ,event " << event << endl;
    }
    

    
    //--- last cut of the cut-based SELECTION
  }//end of loop over SkimEvent
  

  


}

BorisSimpleYielder::~BorisSimpleYielder(){
  using namespace std;
  cout << endl;
  cout << "======= SUMMARY YIELDS======" << endl;
  for(int index=0;index<numberCuts;index++){
    cout << "cut" << index << ": " <<  yieldCounter[index] << endl;
  }
}


void 
BorisSimpleYielder::printout(const edm::EventBase& event,std::vector<reco::SkimEvent>::const_iterator it){
  std::cout << "type,event: " << it->hypo() << " , " 
	    << event.id().event() << " , " 
	    << event.id().luminosityBlock () << " , " 
	    << event.id().run () << std::endl;
}
