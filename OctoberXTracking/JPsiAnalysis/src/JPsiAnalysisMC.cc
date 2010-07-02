#include "OctoberXTracking/JPsiAnalysis/interface/JPsiAnalysisMC.h"

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

JPsiAnalysisMC::JPsiAnalysisMC(const edm::ParameterSet& iConfig)
{
   produces<Hypotheses>();

   //if do put with a label
   //produces<ExampleData2>("label");

   //now do what ever other initialization is needed
  
}


JPsiAnalysisMC::~JPsiAnalysisMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JPsiAnalysisMC::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace ROOT::Math::VectorUtil;
  typedef Candidate::LorentzVector LorentzVector;


  std::auto_ptr<Hypotheses> hypoOut(new Hypotheses);
  
  Handle< vector<Track> > tracks;
  iEvent.getByLabel("generalTracks",tracks);

  
  Handle< vector<GenParticle> > particles;
  iEvent.getByLabel("genParticles",particles);
  
  
  //int pdgId = 0;  //ppMuX
  int pdgId = 443;  //JPsi
  //int pdgId = 553;  //Upsilon
  //int pdgId = 100553;  //Upsilon2S
  //int pdgId = 100553;  //Upsilon3S

  LorentzVector simJPsi, simLeg1, simLeg2; 
  vector<int> charges;
  int charge1,charge2;

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
	charges.push_back(it->charge());
      }
      
    }
    if(tmpSimLegs.size() ==2){ 
      if(tmpSimLegs[0].pt() <= tmpSimLegs[1].pt()){
	simLeg1 = tmpSimLegs[0]; charge1 = charges[0];
	simLeg2 = tmpSimLegs[1]; charge2 = charges[1];
      }else{
	simLeg1 = tmpSimLegs[1]; charge1 = charges[1];
	simLeg2 = tmpSimLegs[0]; charge2 = charges[0];
      }
    }else{
      //doing nothing
    }
  }
  
  LorentzVector LV1,LV2;
  Track tk1,tk2;
  double deltaR1_tmp(1000.),deltaR2_tmp(1000.);
  double deltaPtMax(.04); //cut at 4%

  for(vector<Track>::const_iterator it = tracks->begin();
      it!=tracks->end();++it){
    double tkE = sqrt(it->p()*it->p() + 0.1057*0.1057 );
    LorentzVector tk4V(it->px(),it->py(),it->pz(),tkE);	     
    
    double deltaR;
    double deltaPt;

    deltaR =  DeltaR(tk4V.Vect(),simLeg1.Vect() );
    deltaPt  = fabs( (tk4V.pt() - simLeg1.pt())/simLeg1.pt());
    if(deltaR < deltaR1_tmp && deltaPt < deltaPtMax && it->charge() == charge1){
      deltaR1_tmp = deltaR;
      //cout << "deltaR1_tmp: " << deltaR1_tmp << endl;
      LV1=tk4V;
      tk1=*it;
    }

    deltaR =  DeltaR(tk4V.Vect(),simLeg2.Vect() );
    deltaPt  = fabs( (tk4V.pt() - simLeg2.pt())/simLeg2.pt());
    if(deltaR < deltaR2_tmp && deltaPt < deltaPtMax && it->charge() == charge2){
      deltaR2_tmp = deltaR;
      //cout << "deltaR2_tmp: " << deltaR2_tmp << endl;
      LV2=tk4V;
      tk2=*it;
    }
  }
  
  if(deltaR1_tmp < 0.1 && deltaR2_tmp < 0.1){
    QuarkoniaHypothesis qh; 
         
    LorentzVector jpsi = LV1+LV2;
    qh.setRecoJPsi(jpsi);

    qh.setRecoLegTk(&tk1,1);
    qh.setRecoLegTk(&tk2,2);
    qh.setSimJPsi(simJPsi);
    qh.setSimLeg(simLeg1,1); qh.setSimLeg(simLeg2,2); 
    hypoOut->push_back(qh);
  }
   
  iEvent.put(hypoOut);

}

// ------------ method called once each job just before starting event loop  ------------
void 
JPsiAnalysisMC::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiAnalysisMC::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiAnalysisMC);
