#include "OctoberXTracking/JPsiAnalysis/interface/JPsiAnalysis.h"

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

#include "Math/VectorUtil.h"

#include "OctoberXTracking/DataFormats/interface/QuarkoniaHypothesis.h"

JPsiAnalysis::JPsiAnalysis(const edm::ParameterSet& iConfig)
{
   produces<Hypotheses>();

   //if do put with a label
   //produces<ExampleData2>("label");

   //now do what ever other initialization is needed
  
}


JPsiAnalysis::~JPsiAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JPsiAnalysis::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;


  std::auto_ptr<Hypotheses> hypoOut(new Hypotheses);
  
  //  Handle<TriggerResults> trigger;
  //iEvent.getByLabel("TriggerResults",trigger);

  Handle< vector<Track> > tracks;
  iEvent.getByLabel("generalTracks",tracks);

  Handle< vector<Muon> > muons;
  iEvent.getByLabel("muons",muons);

  
  Handle< vector<GenParticle> > particles;
  iEvent.getByLabel("genParticles",particles);

  

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter;

  
  
  //int pdgId = 0;  //ppMuX
  int pdgId = 443;  //JPsi
  //int pdgId = 553;  //Upsilon
  //int pdgId = 100553;  //Upsilon2S
  //int pdgId = 100553;  //Upsilon3S

  LorentzVector simJPsi, simLeg1, simLeg2; 
  
  
  if(pdgId) {
    vector<LorentzVector>  tmpSimLegs;
    for(vector<GenParticle>::const_iterator it = particles->begin();
	it!=particles->end(); ++it){
      
      if( it->pdgId() == pdgId && it->status() == 2){      
	//cout << "sim mass, status: " << it->p4().mass() << " , "
	//	   << it->status() << endl;
	simJPsi = it->p4();
      
      }
      
      if(it->mother() && it->mother()->pdgId() == pdgId && abs(it->pdgId()) == 13){
	//cout << "doughter id: " << it->pdgId() << endl;
	tmpSimLegs.push_back(it->p4());
      }
      
    }
    if(tmpSimLegs.size() ==2){ 
      if(tmpSimLegs[0].pt() <= tmpSimLegs[1].pt()){
	simLeg1 = tmpSimLegs[0];
	simLeg2 = tmpSimLegs[1];
      }else{
	simLeg1 = tmpSimLegs[1];
	simLeg2 = tmpSimLegs[0];
      }
    }else{
      //doing nothing
    }
  }
  

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
      
      if(!(it->track().get() && it2->track().get())) continue;
      QuarkoniaHypothesis qh; 
      
      vector<TransientTrack> t_tks;
      t_tks.push_back(theTTBuilder->build(it->track().get()));
      t_tks.push_back(theTTBuilder->build(it2->track().get()));
      TransientVertex myVertex = vtxFitter.vertex(t_tks);
      qh.setRecoVertex(myVertex);      

     
      qh.setRecoJPsi(jpsi);
      if(it->pt()<it2->pt()){
	qh.setRecoLeg(&*it,1);
	qh.setRecoLeg(&*it2,2);
	qh.setRecoLegTk(it->track().get(),1);
	qh.setRecoLegTk(it2->track().get(),2);
      }else{
	qh.setRecoLeg(&*it,2);
	qh.setRecoLeg(&*it2,1);
	qh.setRecoLegTk(it->track().get(),2);
	qh.setRecoLegTk(it2->track().get(),1);
      }
      qh.setSimJPsi(simJPsi);
      qh.setSimLeg(simLeg1,1); qh.setSimLeg(simLeg2,2); 
      hypoOut->push_back(qh);
    }
  }
  

  iEvent.put(hypoOut);

  
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
JPsiAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiAnalysis::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiAnalysis);
