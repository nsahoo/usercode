#include "OctoberXTracking/JPsiAnalysis/interface/JPsiAnalysisPAT.h"

//Headers for the data items
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/Candidate/interface/LeafCandidate.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TMath.h"
#include "Math/VectorUtil.h"

//#include "OctoberXTracking/DataFormats/interface/QuarkoniaHypothesis.h"
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

JPsiAnalysisPAT::JPsiAnalysisPAT(const edm::ParameterSet& iConfig)
{
  produces<pat::CompositeCandidateCollection>();  
}


JPsiAnalysisPAT::~JPsiAnalysisPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JPsiAnalysisPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;


  std::auto_ptr<pat::CompositeCandidateCollection> hypoOut(new pat::CompositeCandidateCollection);
  
  //  Handle<TriggerResults> trigger;
  //iEvent.getByLabel("TriggerResults",trigger);

  Handle< vector<Track> > tracks;
  iEvent.getByLabel("generalTracks",tracks);

  Handle< vector<Muon> > muons;
  iEvent.getByLabel("muons",muons);

  /*
  Handle< vector<GenParticle> > particles;
  iEvent.getByLabel("genParticles",particles);
  */
  

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter;

  
  /*
    vector<LeafCandidate> extraCands;
    //cleaning track collection
    for(vector<Track>::const_iterator it = tracks->begin();
    it!=tracks->end();++it){	  
    bool matched(false);
    for(vector<Muon>::const_iterator it2 = muons->begin();
    it2!=muons->end();++it2){
    if(&(*it) == it2->track().get()) matched = true; 
    }
    if(matched) continue;
    if(it->found()<12) continue;
    double tE = sqrt(it->p()*it->p() +0.1057*0.1057 );
    LorentzVector tLV(it->px(),it->py(),it->pz(),tE);	     
    extraCands.push_back(LeafCandidate(it->charge(),tLV));
    }
  */
  
  /*
  // muons
  vector<LeafCandidate> muonCands;
  for(vector<Muon>::const_iterator it = muons->begin();
  it!=muons->end();++it){	 
  muonCands.push_back(*it);
  cout << "muon type: " << it->type() << endl;
  cout << "is global: " << it->isGlobalMuon() << endl;
  cout << "is tracker: " << it-> isTrackerMuon() << endl;
  cout << "is sta: " << it->isStandAloneMuon() << endl;
  cout << "is calo: " << it->isCaloMuon() << endl;
  }
  */

  // JPsi candidates only from muons
  for(vector<Muon>::const_iterator it = muons->begin();
      it!=muons->end();++it){
    for(vector<Muon>::const_iterator it2 = it+1;
	it2!=muons->end();++it2){
      //if(it->charge() == it2->charge()) continue;
      LorentzVector jpsi = it->p4() + it2->p4();
      //cout << "jpsi: " << jpsi.mass() << endl;
      
      // what is the purpose of this line?? I don't remember
      if(!(it->track().get() && it2->track().get())) continue; 
      //QuarkoniaHypothesis qh; 
      pat::CompositeCandidate myCand;


      vector<TransientTrack> t_tks;
      t_tks.push_back(theTTBuilder->build(it->track().get()));
      t_tks.push_back(theTTBuilder->build(it2->track().get()));
      TransientVertex myVertex = vtxFitter.vertex(t_tks);
      float vChi2 = myVertex.totalChiSquared();
      float vNDF  = myVertex.degreesOfFreedom();
      float vProb(TMath::Prob(vChi2,(int)vNDF));

      myCand.addUserFloat("vNChi2",vChi2/vNDF);
      myCand.addUserFloat("vProb",vProb);


      pat::Muon myMuonA(*it);  myMuonA.embedTrack();
      pat::Muon myMuonB(*it2); myMuonB.embedTrack();

      //qh.setRecoVertex(myVertex);      

      //qh.setRecoJPsi(jpsi);
      if(it->pt()<it2->pt()){
	myCand.addDaughter(myMuonA,"muon1");
	myCand.addDaughter(myMuonB,"muon2");
      }else{
	myCand.addDaughter(myMuonA,"muon2");
	myCand.addDaughter(myMuonB,"muon1");	
      }

      myCand.setP4(jpsi);
      
      hypoOut->push_back(myCand);
    }
  }
  

  iEvent.put(hypoOut);

  
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
JPsiAnalysisPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiAnalysisPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiAnalysisPAT);
