#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "WWAnalysis/AnalysisStep/interface/BorisSimpleYielder.h"
#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"


/// default constructor
BorisSimpleYielder::BorisSimpleYielder(const edm::ParameterSet& cfg, TFileDirectory& fs): 
  edm::BasicAnalyzer::BasicAnalyzer(cfg, fs),
  skimEventsTag_(cfg.getParameter<edm::InputTag>("skimEventTag")),
  pdType_(cfg.getParameter<int>("pdType"))
{
  //  hists_["dphi"  ] = fs.make<TH1F>("dphi"  , "dphi"  ,  100,  0., M_PI);
}

/// everything that needs to be done during the event loop
void 
BorisSimpleYielder::analyze(const edm::EventBase& event)
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  using pat::Muon;
  
  using namespace std;


  reco::SkimEvent::primaryDatasetType pdType;
  if(pdType_==1){
    pdType = reco::SkimEvent::SingleMuon;
  }
  if(pdType_==2)
    pdType = reco::SkimEvent::DoubleMuon;
  if(pdType_==3)
    pdType = reco::SkimEvent::MuEG;
  if(pdType_==4)
    pdType = reco::SkimEvent::DoubleElectron;



  edm::Handle<vector<reco::SkimEvent> > skimEventsMuMu;
  event.getByLabel(edm::InputTag("wwmumuCONVLHT"),skimEventsMuMu);

  edm::Handle<vector<reco::SkimEvent> > skimEventsMuEl;
  event.getByLabel(edm::InputTag("wwmuelCONVLHT"),skimEventsMuEl);

  edm::Handle<vector<reco::SkimEvent> > skimEventsElMu;
  event.getByLabel(edm::InputTag("wwelmuCONVLHT"),skimEventsElMu);

  edm::Handle<vector<reco::SkimEvent> > skimEventsElEl;
  event.getByLabel(edm::InputTag("wwelelCONVLHT"),skimEventsElEl);

  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(edm::InputTag("TriggerResults","","HLT"),triggerResults);

  for(vector<reco::SkimEvent>::const_iterator it = skimEventsMuMu->begin();
      it != skimEventsMuMu->end(); ++it){
    if(!cut(it,pdType,*triggerResults,event)) continue;
    printout(event,it);
  }

  for(vector<reco::SkimEvent>::const_iterator it = skimEventsMuEl->begin();
      it != skimEventsMuEl->end(); ++it){
    if(!cut(it,pdType,*triggerResults,event)) continue;
    printout(event,it);
  }

  for(vector<reco::SkimEvent>::const_iterator it = skimEventsElMu->begin();
      it != skimEventsElMu->end(); ++it){
    if(!cut(it,pdType,*triggerResults,event)) continue;
    printout(event,it);
  }

  for(vector<reco::SkimEvent>::const_iterator it = skimEventsElEl->begin();
      it != skimEventsElEl->end(); ++it){
    if(!cut(it,pdType,*triggerResults,event)) continue;
    printout(event,it);
  }





}

bool BorisSimpleYielder::cut(std::vector<reco::SkimEvent>::const_iterator it, 
			     reco::SkimEvent::primaryDatasetType pdType,
			     edm::TriggerResults triggerResults,
			     const edm::EventBase& event){
  using namespace std;
  vector<string> pathsSingleMu; pathsSingleMu.push_back("HLT_Mu24_v");
  vector<string> pathsDoubleMu; pathsDoubleMu.push_back("HLT_DoubleMu7_v");
  vector<string> pathsDoubleEl; 
  pathsDoubleEl.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v");
  pathsDoubleEl.push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v");
  vector<string> pathsElMu; 
  pathsElMu.push_back("HLT_Mu8_Ele17_CaloIdL_v");
  pathsElMu.push_back("HLT_Mu17_Ele8_CaloIdL_v");

  TriggerBitChecker checkerSingleMu(pathsSingleMu);
  TriggerBitChecker checkerDoubleMu(pathsDoubleMu);
  TriggerBitChecker checkerDoubleEl(pathsDoubleEl);
  TriggerBitChecker checkerElMu(pathsElMu);




  if(!(it->q(0)*it->q(1)<0 && !it->isSTA(0) && !it->isSTA(1) && it->leptEtaCut(2.4,2.5)
       && it->ptMin() > 10 && it->ptMax() > 20)) return false;
  if(!(it->passesIP())) return false;
  if(!(it->nExtraLep(10) == 0)) return false;
  if(!(it->pfMet() > 20)) return false;
  if(!(it->mll()>12 ) ) return false;
  if(it->hypo()==6 || it->hypo()==3){ //same flavor
    if(!(abs(it->mll()-91.1876)>15.0 ) ) return false;
    if(!(min(it->projPfMet(),it->projChargedMet()) > 35 ) ) return false;
  }    
  if(it->hypo()==4 || it->hypo()==5){ //opposite flavor
    if(!(min(it->projPfMet(),it->projChargedMet()) > 20 ) ) return false;
  }

  if(it->hypo()==3 || it->hypo()==5){
    if(!(it->ptMin()>15    )) return false;
  }

  if(!(it->nCentralJets( 30.0, 5.0, 1 ) == 0 ) ) return false;
  if(!(it->nSoftMu( 3.0 ) == 0 ) ) return false;
  if(!(it->bTaggedJetsUnder( 30.0,  2.1 ) == 0 ) ) return false;
  if(!(it->mll()<45 ) ) return false;
  if(!(it->ptMax()>25 ) ) return false;
  if(!(it->ptMin()>10 ) ) return false;
  if(!(it->dPhillInDegrees()<90 ) ) return false;
  if(!(2.*(it->gammaMRStar())/130. > 0.5 && 2.*(it->gammaMRStar())/130 < 1.1)) return false;

  //---------
  if(!(it->triggerMatchingCut(pdType))) return false;
  bool boolSingleMu(checkerSingleMu.check(event,triggerResults));
  bool boolDoubleMu(checkerDoubleMu.check(event,triggerResults));  
  bool boolDoubleEl(checkerDoubleEl.check(event,triggerResults));  
  bool boolElMu(checkerElMu.check(event,triggerResults));  
  bool result(false);

  if(it->hypo()==6){//mumu
    if(pdType==reco::SkimEvent::DoubleMuon)
      result = boolDoubleMu;   
    if(pdType==reco::SkimEvent::SingleMuon)
      result = boolSingleMu && (!boolDoubleMu);
  }
  
  if(it->hypo()==4 || it->hypo()==5 ){//elmu or muel
    if(pdType==reco::SkimEvent::SingleMuon)
      result = boolSingleMu;
    if(pdType==reco::SkimEvent::MuEG)
      result = boolElMu && (!boolSingleMu);
  }
  
  if(it->hypo()==3){//elel
    if(pdType==reco::SkimEvent::DoubleElectron)
      result = boolDoubleEl;       
  }
  if(!result) return false;
  // ------------

  return true;
}

void 
BorisSimpleYielder::printout(const edm::EventBase& event,std::vector<reco::SkimEvent>::const_iterator it){
  std::cout << "type,event: " << it->hypo() << " , " 
	    << event.id().event() << " , " 
	    << event.id().luminosityBlock () << " , " 
	    << event.id().run () << std::endl;
}
