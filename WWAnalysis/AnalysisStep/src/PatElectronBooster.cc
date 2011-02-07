// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>


#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/GaussianSumVertexFit/interface/GsfVertexTrackCompatibilityEstimator.h"
#include "WWAnalysis/Tools/interface/VertexReProducer.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"

//for sim-reco IP bias
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//
// class declaration
//

class PatElectronBooster : public edm::EDProducer {
    public:
        explicit PatElectronBooster(const edm::ParameterSet&);
        ~PatElectronBooster();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag electronTag_;
        edm::InputTag vertexTag_;
        // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PatElectronBooster::PatElectronBooster(const edm::ParameterSet& iConfig) :
    electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
    vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag"))
{
    produces<pat::ElectronCollection>();  
}


PatElectronBooster::~PatElectronBooster()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
    void
PatElectronBooster::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    ParameterSet psetMatcher;
    psetMatcher.addParameter<double>("maxDPtRel",999.0);
    psetMatcher.addParameter<double>("maxDeltaR",0.015);

    reco::MatchByDRDPt<reco::Track,reco::Track> ttMatcher(psetMatcher);
    reco::MatchByDRDPt<reco::Track,reco::GsfTrack> tgMatcher(psetMatcher);

    ESHandle<TransientTrackBuilder> theTTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
    iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP",parametersDefinerTP); 

    Handle<pat::ElectronCollection> electrons;
    iEvent.getByLabel(electronTag_,electrons);

    Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(vertexTag_,vertices);

    Handle<reco::BeamSpot> bs;
    iEvent.getByLabel(InputTag("offlineBeamSpot"),bs);


    VertexReProducer revertex(vertices, iEvent);
    Handle<reco::TrackCollection> pvtracks;   iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
    Handle<reco::BeamSpot>        pvbeamspot; iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);


    //   if(vertices->size()==0) return;
    //   reco::Vertex vertex = vertices->front();

    reco::Vertex vertex;
    reco::Vertex vertex2;


    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);
    for(pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
      pat::Electron clone = *ele;

      // -------- here is the comparison with the TP
      /*
      TrackingParticle tp;
      tp.setP4(ele->genParticleRef()->p4());
      tp.setVertex(ele->genParticleRef()->vertex(),0.);
      tp.setThreeCharge(ele->genParticleRef()->charge()*3.0);

      ParticleBase::Vector momentumTP = parametersDefinerTP->momentum(iEvent,iSetup,tp);
      ParticleBase::Point vertexTP = parametersDefinerTP->vertex(iEvent,iSetup,tp);
      double chargeTP = tp.charge();

      double qoverpSim = chargeTP/sqrt(momentumTP.x()*momentumTP.x()+momentumTP.y()*momentumTP.y()+momentumTP.z()*momentumTP.z());
      double lambdaSim = M_PI/2-momentumTP.theta();
      double phiSim    = momentumTP.phi();
      double dxySim    = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
      double dzSim     = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y())/sqrt(momentumTP.perp2()) 
	* momentumTP.z()/sqrt(momentumTP.perp2());

      //double dxy = ele->gsfTrack()->dxy(vertices->front().position());
      double dxy = ele->gsfTrack()->dxy(bs->position());
      double dz = ele->gsfTrack()->dz(bs->position());

      double dxyPull = (dxy-dxySim)/ele->gsfTrack()->dxyError(); 
      double dzPull = (dz-dzSim)/ele->gsfTrack()->dzError(); 

      clone.addUserFloat(std::string("dxy"),dxy-dxySim);
      clone.addUserFloat(std::string("dz"),dz-dzSim);
      clone.addUserFloat(std::string("dxyPull"),dxyPull);
      clone.addUserFloat(std::string("dzPull"),dzPull);
      */

      // ------- here I add the information about the IP significance wrt the PV
      reco::TrackCollection newTkCollection;
      bool foundMatch(false);
      for(reco::TrackCollection::const_iterator itk = pvtracks->begin(); itk != pvtracks->end(); ++itk){    
	bool refMatching;
	if(ele->closestCtfTrack().ctfTrack.isNonnull())
	  refMatching = (&*itk == &*(ele->closestCtfTrack().ctfTrack) );
	else
	  refMatching = false;	
	if(refMatching && ele->closestCtfTrack().shFracInnerHits > 0.5){
	  float shFraction = ele->closestCtfTrack().shFracInnerHits;
	  //cout << "sharedHit fraction: " << shFraction << endl;
	  newTkCollection.assign(pvtracks->begin(), itk);
	  newTkCollection.insert(newTkCollection.end(),itk+1,pvtracks->end());
	  foundMatch = true;
	  break;
	}
      }

      if(!foundMatch) {
	cout << "WARNING: matching to electron was not found" << endl;
	if(vertices->empty()) {
	  vertex = reco::Vertex(reco::Vertex::Point(bs->position().x(),bs->position().y(),bs->position().z()),
				reco::Vertex::Error());
	} else {
	  vertex = vertices->front();
	}
      }else{      
	vector<TransientVertex> pvs = revertex.makeVertices(newTkCollection, *pvbeamspot, iSetup) ;
	if(pvs.empty()) {
	  vertex = reco::Vertex(reco::Vertex::Point(bs->position().x(),bs->position().y(),bs->position().z()),
				reco::Vertex::Error());
	} else {
	  vertex = reco::Vertex(pvs.front());
	}
      }

      if(vertices->empty()) {
	vertex2 = reco::Vertex(reco::Vertex::Point(bs->position().x(),bs->position().y(),bs->position().z()),
			      reco::Vertex::Error());
      } else {
	vertex2 = vertices->front();
      }

      reco::TransientTrack tt = theTTBuilder->build(ele->gsfTrack());
      Measurement1D ip = IPTools::absoluteTransverseImpactParameter(tt,vertex).second;
      Measurement1D ip3D = IPTools::absoluteImpactParameter3D(tt,vertex).second;

      Measurement1D ip_2 = IPTools::absoluteTransverseImpactParameter(tt,vertex2).second;
      Measurement1D ip3D_2 = IPTools::absoluteImpactParameter3D(tt,vertex2).second;


      clone.addUserFloat(std::string("tip"),ip.value());
      clone.addUserFloat(std::string("tipErr"),ip.error());
      clone.addUserFloat(std::string("ip"),ip3D.value());
      clone.addUserFloat(std::string("ipErr"),ip3D.error());

      clone.addUserFloat(std::string("tip2"),ip_2.value());
      clone.addUserFloat(std::string("tipErr2"),ip_2.error());
      clone.addUserFloat(std::string("ip2"),ip3D_2.value());
      clone.addUserFloat(std::string("ipErr2"),ip3D_2.error());


      
      // ------- here I add the information about the electron-PV compatibility
      try{
	GsfVertexTrackCompatibilityEstimator vertexTrackCompEst;       
	GsfVertexTrackCompatibilityEstimator::BDpair result = vertexTrackCompEst.estimate(vertex,tt);
	if(result.first)
	  clone.addUserFloat(std::string("vtxComp"),result.second);
	  else
	    clone.addUserFloat(std::string("vtxComp"),99999.);
	//cout << "here4" << endl;
      }catch(...){
	cout << "WARNING, failed to evaluate ele-PV compatibility in PatElectronBooster" << endl;
	clone.addUserFloat(std::string("vtxComp"),99999.);
      }

      try{
	GsfVertexTrackCompatibilityEstimator vertexTrackCompEst;       
	GsfVertexTrackCompatibilityEstimator::BDpair result = vertexTrackCompEst.estimate(vertex2,tt);
	if(result.first)
	  clone.addUserFloat(std::string("vtxComp2"),result.second);
	  else
	    clone.addUserFloat(std::string("vtxComp2"),99999.);
	//cout << "here4" << endl;
      }catch(...){
	cout << "WARNING, failed to evaluate ele-PV compatibility in PatElectronBooster" << endl;
	clone.addUserFloat(std::string("vtxComp2"),99999.);
      }




      // ------- here I add the information about the absolute IP  wrt the PV
      clone.addUserFloat( "dxyPV",clone.gsfTrack()->dxy(vertices->front().position()) );
      clone.addUserFloat( "dzPV",clone.gsfTrack()->dz(vertices->front().position()) );


      pOut->push_back(clone);

    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
    void 
PatElectronBooster::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PatElectronBooster::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronBooster);
