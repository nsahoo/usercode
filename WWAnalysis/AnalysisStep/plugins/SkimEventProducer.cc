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
#include "DataFormats/Common/interface/ValueMap.h"

#include<vector> 
#include "Math/VectorUtil.h"


SkimEventProducer::SkimEventProducer(const edm::ParameterSet& cfg) :
//   branchAlias_(cfg.getParameter<std::string>("branchAlias")),
  hypoType_(reco::SkimEvent::hypoTypeByName(cfg.getParameter<std::string>("hypoType"))),
  triggerTag_(cfg.getParameter<edm::InputTag>("triggerTag")),
  singleMuData_  ( cfg.getParameter<std::vector<std::string> >("singleMuDataPaths") ),
  doubleMuData_  ( cfg.getParameter<std::vector<std::string> >("doubleMuDataPaths") ),
  doubleElData_  ( cfg.getParameter<std::vector<std::string> >("doubleElDataPaths") ),
  muEGData_      ( cfg.getParameter<std::vector<std::string> >("muEGDataPaths") ),
  singleMuMC_    ( cfg.getParameter<std::vector<std::string> >("singleMuMCPaths") ),
  doubleMuMC_    ( cfg.getParameter<std::vector<std::string> >("doubleMuMCPaths") ),
  doubleElMC_    ( cfg.getParameter<std::vector<std::string> >("doubleElMCPaths") ),
  muEGMC_        ( cfg.getParameter<std::vector<std::string> >("muEGMCPaths") )
{
  muTag_          = cfg.getParameter<edm::InputTag>("muTag"     ); 
  elTag_          = cfg.getParameter<edm::InputTag>("elTag"     ); 
  softMuTag_      = cfg.getParameter<edm::InputTag>("softMuTag" ); 
  jetTag_         = cfg.getParameter<edm::InputTag>("jetTag"    ); 
  tagJetTag_      = cfg.getParameter<edm::InputTag>("tagJetTag" ); 
  pfMetTag_       = cfg.getParameter<edm::InputTag>("pfMetTag"  ); 
  tcMetTag_       = cfg.getParameter<edm::InputTag>("tcMetTag"  ); 
  chargedMetTag_  = cfg.getParameter<edm::InputTag>("chargedMetTag" ); 
  vtxTag_         = cfg.getParameter<edm::InputTag>("vtxTag"        ); 
  chCandsTag_     = cfg.getParameter<edm::InputTag>("chCandsTag"    ); 

  if (cfg.exists("sptTag"    )) sptTag_     = cfg.getParameter<edm::InputTag>("sptTag"    ); 
  else                          sptTag_     = edm::InputTag("","","");
  if (cfg.exists("spt2Tag"   )) spt2Tag_    = cfg.getParameter<edm::InputTag>("spt2Tag"   ); 
  else                          spt2Tag_    = edm::InputTag("","","");

  produces<std::vector<reco::SkimEvent> >().setBranchAlias(cfg.getParameter<std::string>("@module_label"));
}

void SkimEventProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<reco::SkimEvent> > skimEvent(new std::vector<reco::SkimEvent> );

    //SkimEvent::hypoType type = SkimEvent::hypoTypeByName(hypoType_);
    
    edm::Handle<pat::JetCollection> jetH;
    iEvent.getByLabel(jetTag_,jetH);

    edm::Handle<pat::JetCollection> tagJetH;
    if(!(tagJetTag_==edm::InputTag(""))) iEvent.getByLabel(tagJetTag_,tagJetH);

    edm::Handle<reco::PFMETCollection> pfMetH;
    iEvent.getByLabel(pfMetTag_,pfMetH);

    edm::Handle<reco::METCollection> tcMetH;
    iEvent.getByLabel(tcMetTag_,tcMetH);

    edm::Handle<edm::ValueMap<reco::PFMET> > chargedMetH;
    iEvent.getByLabel(chargedMetTag_,chargedMetH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxTag_,vtxH);

    edm::Handle<reco::CandidateCollection> candsH;
    iEvent.getByLabel(chCandsTag_,candsH);

    edm::Handle<edm::ValueMap<float> > sptH;
    if(!(sptTag_==edm::InputTag(""))) iEvent.getByLabel(sptTag_,sptH);

    edm::Handle<edm::ValueMap<float> > spt2H;
    if(!(spt2Tag_==edm::InputTag(""))) iEvent.getByLabel(spt2Tag_,spt2H);

    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerTag_,triggerResults);

    // May God have mercy on my soul ...
    std::vector<bool> passBits;
    passBits.push_back( singleMuData_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleMuData_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleElData_.check(iEvent,*triggerResults) );
    passBits.push_back( muEGData_.check(    iEvent,*triggerResults) );
    passBits.push_back( singleMuMC_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleMuMC_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleElMC_.check(iEvent,*triggerResults) );
    passBits.push_back( muEGMC_.check(    iEvent,*triggerResults) );


    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByLabel(muTag_,muons);
    edm::Handle<pat::MuonCollection> softMuH;
    iEvent.getByLabel(softMuTag_,softMuH);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByLabel(elTag_,electrons);

    /*
    edm::Handle<pat::ElectronCollection> extraElH;
    iEvent.getByLabel(extraElTag_,extraElH);
    */


    if(hypoType_==reco::SkimEvent::WWELEL){//ELEL
      for(pat::ElectronCollection::const_iterator ele1=electrons->begin(); ele1!=electrons->end(); ++ele1){
	for(pat::ElectronCollection::const_iterator ele2=ele1+1; ele2!=electrons->end(); ++ele2){
	  skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
	  skimEvent->back().setTriggerBits(passBits);
	  skimEvent->back().setJets(jetH);
//       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
	  if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
	  else                  skimEvent->back().setTagJets(jetH);
	  skimEvent->back().setPFMet(pfMetH);
	  skimEvent->back().setTCMet(tcMetH);
	  skimEvent->back().setChargedMet(chargedMetH->get(0));
	  skimEvent->back().setVertex(vtxH);
	  if(sptH.isValid()) skimEvent->back().setVtxSumPts(sptH);
	  if(spt2H.isValid()) skimEvent->back().setVtxSumPt2s(spt2H);
	  skimEvent->back().setLepton(*ele1);
	  skimEvent->back().setLepton(*ele2);
	  skimEvent->back().setChargedMetSmurf(doChMET(candsH,&*ele1,&*ele2));
	  for(pat::ElectronCollection::const_iterator ele3=electrons->begin(); ele3!=electrons->end(); ++ele3){
	    //float delta1 = ROOT::Math::VectorUtil::DeltaR(ele3->p4(),ele1->p4());
	    //float delta2 = ROOT::Math::VectorUtil::DeltaR(ele3->p4(),ele2->p4());
	    //if(ele3!=ele1 && ele3!=ele2 && delta1 > 0.1 && delta2 > 0.1)
	    if(ele3!=ele1 && ele3!=ele2) //I would prefere the line above. Now synch with ST selection
	      skimEvent->back().setExtraLepton(*ele3);
	  }
	  for(pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(mu->p4(),ele1->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(mu->p4(),ele2->p4());
	    if(delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*mu);
	  }
	  for(pat::MuonCollection::const_iterator smu=softMuH->begin(); smu!=softMuH->end(); ++smu){
	    skimEvent->back().setSoftMuon(*smu);
	  }
	}
      }//end loop on main lepton collection
    }else if(hypoType_ == reco::SkimEvent::WWELMU){
      for(pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
	for(pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
      if( mu->pt() >= ele->pt() ) continue;
	  skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
      skimEvent->back().setTriggerBits(passBits);
	  skimEvent->back().setJets(jetH);
//       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
	  if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
	  else                  skimEvent->back().setTagJets(jetH);
	  skimEvent->back().setPFMet(pfMetH);
	  skimEvent->back().setTCMet(tcMetH);
	  skimEvent->back().setChargedMet(chargedMetH->get(0));
	  skimEvent->back().setVertex(vtxH);
	  if(sptH.isValid()) skimEvent->back().setVtxSumPts(sptH);
	  if(spt2H.isValid()) skimEvent->back().setVtxSumPt2s(spt2H);
	  skimEvent->back().setLepton(*ele);
	  skimEvent->back().setLepton(*mu);
	  skimEvent->back().setChargedMetSmurf(doChMET(candsH,&*ele,&*mu));
	  for(pat::ElectronCollection::const_iterator ele2=electrons->begin(); ele2!=electrons->end(); ++ele2){
	    //float delta1 = ROOT::Math::VectorUtil::DeltaR(ele2->p4(),ele->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(ele2->p4(),mu->p4());
	    //if(ele2!=ele && delta1 > 0.1 && delta2 > 0.1)
	    if(ele2!=ele && delta2 > 0.1)//I would prefere the line above. Now synch with ST selection
	      skimEvent->back().setExtraLepton(*ele2);
	  }
	  for(pat::MuonCollection::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(mu2->p4(),ele->p4());
	    //float delta2 = ROOT::Math::VectorUtil::DeltaR(mu2->p4(),mu->p4());
	    //if(mu2!=mu && delta1 > 0.1 && delta2 > 0.1)
	    if(mu2!=mu && delta1 > 0.1)//I would prefere the line above. Now synch with ST selection
	      skimEvent->back().setExtraLepton(*mu2);
	  }
	  for(pat::MuonCollection::const_iterator smu=softMuH->begin(); smu!=softMuH->end(); ++smu){
	    if(smu->pt() != mu->pt() && smu->eta() != mu->eta())
	      skimEvent->back().setSoftMuon(*smu);
	  }
	}
      }//end loop on main lepton collection
      
    }else if(hypoType_ == reco::SkimEvent::WWMUEL){
      for(pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
	for(pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
      if( mu->pt() < ele->pt() ) continue;
	  skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
      skimEvent->back().setTriggerBits(passBits);
	  skimEvent->back().setJets(jetH);
//       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
	  if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
	  else                  skimEvent->back().setTagJets(jetH);
	  skimEvent->back().setPFMet(pfMetH);
	  skimEvent->back().setTCMet(tcMetH);
	  skimEvent->back().setChargedMet(chargedMetH->get(0));
	  skimEvent->back().setVertex(vtxH);
	  if(sptH.isValid()) skimEvent->back().setVtxSumPts(sptH);
	  if(spt2H.isValid()) skimEvent->back().setVtxSumPt2s(spt2H);
	  skimEvent->back().setLepton(*mu);
	  skimEvent->back().setLepton(*ele);
	  skimEvent->back().setChargedMetSmurf(doChMET(candsH,&*mu,&*ele));
	  for(pat::ElectronCollection::const_iterator ele2=electrons->begin(); ele2!=electrons->end(); ++ele2){
	    //float delta1 = ROOT::Math::VectorUtil::DeltaR(ele2->p4(),ele->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(ele2->p4(),mu->p4());
	    //if(ele2!=ele && delta1 > 0.1 && delta2 > 0.1)
	    if(ele2!=ele && delta2 > 0.1)//I would prefere the line above. Now synch with ST selection
	      skimEvent->back().setExtraLepton(*ele2);
	  }
	  for(pat::MuonCollection::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(mu2->p4(),ele->p4());
	    //float delta2 = ROOT::Math::VectorUtil::DeltaR(mu2->p4(),mu->p4());
	    //if(mu2!=mu && delta1 > 0.1 && delta2 > 0.1)
	    if(mu2!=mu && delta1 > 0.1)//I would prefere the line above. Now synch with ST selection
	      skimEvent->back().setExtraLepton(*mu2);
	  }
	  for(pat::MuonCollection::const_iterator smu=softMuH->begin(); smu!=softMuH->end(); ++smu){
	    if(smu->pt() != mu->pt() && smu->eta() != mu->eta())
	      skimEvent->back().setSoftMuon(*smu);
	  }
	}
      }//end loop on main lepton collection
      
    }else if(hypoType_==reco::SkimEvent::WWMUMU){//MUMU
      for(pat::MuonCollection::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1){
	for(pat::MuonCollection::const_iterator mu2=mu1+1; mu2!=muons->end(); ++mu2){
      	  skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
          skimEvent->back().setTriggerBits(passBits);
	  skimEvent->back().setJets(jetH);
//       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
	  if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
	  else                  skimEvent->back().setTagJets(jetH);
	  skimEvent->back().setPFMet(pfMetH);
	  skimEvent->back().setTCMet(tcMetH);
	  skimEvent->back().setChargedMet(chargedMetH->get(0));
	  skimEvent->back().setVertex(vtxH);
	  if(sptH.isValid()) skimEvent->back().setVtxSumPts(sptH);
	  if(spt2H.isValid()) skimEvent->back().setVtxSumPt2s(spt2H);
	  skimEvent->back().setLepton(*mu1);
	  skimEvent->back().setLepton(*mu2);
	  skimEvent->back().setChargedMetSmurf(doChMET(candsH,&*mu1,&*mu2));
	  for(pat::MuonCollection::const_iterator mu3=muons->begin(); mu3!=muons->end(); ++mu3){
	    //float delta1 = ROOT::Math::VectorUtil::DeltaR(mu3->p4(),mu1->p4());
	    //float delta2 = ROOT::Math::VectorUtil::DeltaR(mu3->p4(),mu2->p4());
	    //if(mu3!=mu1 && mu3!=mu2 && delta1 > 0.1 && delta2 > 0.1)
	    if(mu3!=mu1 && mu3!=mu2)//I would prefere the line above. Now synch with ST selection
	      skimEvent->back().setExtraLepton(*mu3);
	  }
	  for(pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
	    float delta1 = ROOT::Math::VectorUtil::DeltaR(ele->p4(),mu1->p4());
	    float delta2 = ROOT::Math::VectorUtil::DeltaR(ele->p4(),mu2->p4());
	    if(delta1 > 0.1 && delta2 > 0.1)
	      skimEvent->back().setExtraLepton(*ele);
	  }

	  for(pat::MuonCollection::const_iterator smu=softMuH->begin(); smu!=softMuH->end(); ++smu){
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



reco::MET SkimEventProducer::doChMET(edm::Handle<reco::CandidateCollection> candsH,
				  const reco::Candidate* cand1,const reco::Candidate* cand2){
  using namespace std;
  reco::Candidate::LorentzVector totalP4;
  for(reco::CandidateCollection::const_iterator it= candsH->begin(); it!=candsH->end(); ++it){
    if(fabs(ROOT::Math::VectorUtil::DeltaR(it->p4(),cand1->p4())) <=0.1) continue;
    if(fabs(ROOT::Math::VectorUtil::DeltaR(it->p4(),cand2->p4())) <=0.1) continue;
    totalP4 += it->p4();
  }
  totalP4 +=cand1->p4();
  totalP4 +=cand2->p4();
  reco::Candidate::LorentzVector invertedP4(-totalP4);
  reco::MET met(invertedP4,reco::Candidate::Point(0,0,0));
  return met;
}

DEFINE_FWK_MODULE(SkimEventProducer);
