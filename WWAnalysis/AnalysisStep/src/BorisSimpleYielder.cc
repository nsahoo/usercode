#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "WWAnalysis/AnalysisStep/interface/BorisSimpleYielder.h"
#include <WWAnalysis/DataFormats/interface/SkimEvent.h>


/// default constructor
BorisSimpleYielder::BorisSimpleYielder(const edm::ParameterSet& cfg, TFileDirectory& fs): 
  edm::BasicAnalyzer::BasicAnalyzer(cfg, fs),
  skimEventsTag_(cfg.getParameter<edm::InputTag>("skimEventTag"))
{
  hists_["dphi"  ] = fs.make<TH1F>("dphi"  , "dphi"  ,  100,  0., M_PI);
}

/// everything that needs to be done during the event loop
void 
BorisSimpleYielder::analyze(const edm::EventBase& event)
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  using pat::Muon;
  
  using namespace std;

  edm::Handle<vector<reco::SkimEvent> > skimEvents;
  event.getByLabel(skimEventsTag_,skimEvents);

  for(vector<reco::SkimEvent>::const_iterator it = skimEvents->begin();
      it != skimEvents->end(); ++it){
    hists_["dphi" ]->Fill( it->dPhill() );
  }

}
