#include "WWAnalysis/AnalysisStep/interface/SkimEventProducer.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include<vector> 
#include "Math/VectorUtil.h"


SkimEventProducer::SkimEventProducer(const edm::ParameterSet& cfg) :
  branchAlias_(cfg.getParameter<std::string>("branchAlias")),
  hypoType_(SkimEvent::hypoTypeByName(cfg.getParameter<std::string>("hypoType")))
{
    if (cfg.exists("muTag"     )) muTag_      = cfg.getParameter<edm::InputTag>("muTag"     ); 
    else                          muTag_      = edm::InputTag("","","");
    if (cfg.exists("elTag"     )) elTag_      = cfg.getParameter<edm::InputTag>("elTag"     ); 
    else                          elTag_      = edm::InputTag("","","");
    if (cfg.exists("extraMuTag")) extraMuTag_ = cfg.getParameter<edm::InputTag>("extraMuTag"); 
    else                          extraMuTag_ = edm::InputTag("","","");
    if (cfg.exists("extraElTag")) extraElTag_ = cfg.getParameter<edm::InputTag>("extraElTag"); 
    else                          extraElTag_ = edm::InputTag("","","");
    if (cfg.exists("jetTag"    )) jetTag_     = cfg.getParameter<edm::InputTag>("jetTag"    ); 
    else                          jetTag_     = edm::InputTag("","","");
    if (cfg.exists("pfMetTag"  )) pfMetTag_   = cfg.getParameter<edm::InputTag>("pfMetTag"  ); 
    else                          pfMetTag_   = edm::InputTag("","","");
    if (cfg.exists("tcMetTag"  )) tcMetTag_   = cfg.getParameter<edm::InputTag>("tcMetTag"  ); 
    else                          tcMetTag_   = edm::InputTag("","","");
    if (cfg.exists("vtxTag"    )) vtxTag_     = cfg.getParameter<edm::InputTag>("vtxTag"    ); 
    else                          vtxTag_     = edm::InputTag("","","");

    produces<std::vector<reco::SkimEvent> >().setBranchAlias(branchAlias_);
}

void SkimEventProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<reco::SkimEvent> > skimEvent(new std::vector<reco::SkimEvent> );

    //SkimEvent::hypoType type = SkimEvent::hypoTypeByName(hypoType_);
    
    edm::Handle<pat::JetCollection> jetH;
    iEvent.getByLabel(jetTag_,jetH);

    edm::Handle<reco::PFMETCollection> pfMetH;
    iEvent.getByLabel(pfMetTag_,pfMetH);

    edm::Handle<reco::METCollection> tcMetH;
    iEvent.getByLabel(tcMetTag_,tcMetH);


    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxTag_,vtxH);


    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByLabel(muTag_,muons);
    edm::Handle<pat::MuonCollection> extraMuH;
    iEvent.getByLabel(extraMuTag_,extraMuH);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByLabel(elTag_,electrons);

    /*
    edm::Handle<pat::ElectronCollection> extraElH;
    iEvent.getByLabel(extraElTag_,extraElH);
    */


    if(hypoType_==3){//ELEL
      for(pat::ElectronCollection::const_iterator ele1=electrons->begin(); ele1!=electrons->end(); ++ele1){
	for(pat::ElectronCollection::const_iterator ele2=ele1+1; ele2!=electrons->end(); ++ele2){
	  skimEvent->push_back( *(new SkimEvent(hypoType_) ) );      
	  skimEvent->back().setJets(jetH);
	  skimEvent->back().setPFMet(pfMetH);
	  skimEvent->back().setTCMet(tcMetH);
	  skimEvent->back().setVertex(vtxH);
	  skimEvent->back().setLepton(*ele1);
	  skimEvent->back().setLepton(*ele2);
	  for(pat::ElectronCollection::const_iterator ele3=electrons->begin(); ele3!=electrons->end(); ++ele3){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(ele3->p4(),ele1->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(ele3->p4(),ele2->p4());
	    if(ele3!=ele1 && ele3!=ele2 && delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*ele3);
	  }
	  for(pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(mu->p4(),ele1->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(mu->p4(),ele2->p4());
	    if(delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*mu);
	  }
	  for(pat::MuonCollection::const_iterator smu=extraMuH->begin(); smu!=extraMuH->end(); ++smu){
	    skimEvent->back().setSoftMuon(*smu);
	  }
	}
      }//end loop on main lepton collection
    }else if(hypoType_==4){//ELMU
      for(pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
	for(pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
	  skimEvent->push_back( *(new SkimEvent(hypoType_) ) );      
	  skimEvent->back().setJets(jetH);
	  skimEvent->back().setPFMet(pfMetH);
	  skimEvent->back().setTCMet(tcMetH);
	  skimEvent->back().setVertex(vtxH);
	  skimEvent->back().setLepton(*ele);
	  skimEvent->back().setLepton(*mu);
	  for(pat::ElectronCollection::const_iterator ele2=electrons->begin(); ele2!=electrons->end(); ++ele2){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(ele2->p4(),ele->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(ele2->p4(),mu->p4());
	    if(ele2!=ele && delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*ele2);
	  }
	  for(pat::MuonCollection::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(mu2->p4(),ele->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(mu2->p4(),mu->p4());
	    if(mu2!=mu && delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*mu2);
	  }
	  for(pat::MuonCollection::const_iterator smu=extraMuH->begin(); smu!=extraMuH->end(); ++smu){
	    if(smu->pt() != mu->pt() && smu->eta() != mu->eta())
	      skimEvent->back().setSoftMuon(*smu);
	  }
	}
      }//end loop on main lepton collection
      
    }else if(hypoType_==5){//MUMU
      for(pat::MuonCollection::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1){
	for(pat::MuonCollection::const_iterator mu2=mu1+1; mu2!=muons->end(); ++mu2){
      	  skimEvent->push_back( *(new SkimEvent(hypoType_) ) );      
	  skimEvent->back().setJets(jetH);
	  skimEvent->back().setPFMet(pfMetH);
	  skimEvent->back().setTCMet(tcMetH);
	  skimEvent->back().setVertex(vtxH);
	  skimEvent->back().setLepton(*mu1);
	  skimEvent->back().setLepton(*mu2);
	  for(pat::MuonCollection::const_iterator mu3=muons->begin(); mu3!=muons->end(); ++mu3){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(mu3->p4(),mu1->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(mu3->p4(),mu2->p4());
	    if(mu3!=mu1 && mu3!=mu2 && delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*mu3);
	  }
	  for(pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(ele->p4(),mu1->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(ele->p4(),mu2->p4());
	    if(delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*ele);
	  }

	  for(pat::MuonCollection::const_iterator smu=extraMuH->begin(); smu!=extraMuH->end(); ++smu){
	    if(smu->pt() != mu1->pt() && smu->eta() != mu1->eta() &&
	       smu->pt() != mu2->pt() && smu->eta() != mu2->eta())
	      skimEvent->back().setSoftMuon(*smu);
	  }
	}
      }//end loop on main lepton collection
    }else{
      throw cms::Exception("BadInput") 
        << "ERROR: event type " << hypoType_ << " is not known" ;
    }





    iEvent.put(skimEvent);
}

SkimEventProducer::~SkimEventProducer() { } 
void SkimEventProducer::beginJob() { }
void SkimEventProducer::endJob() { }

DEFINE_FWK_MODULE(SkimEventProducer);
