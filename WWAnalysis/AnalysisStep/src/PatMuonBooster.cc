// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/GaussianSumVertexFit/interface/GsfVertexTrackCompatibilityEstimator.h"

//
// class declaration
//

class PatMuonBooster : public edm::EDProducer {
    public:
        explicit PatMuonBooster(const edm::ParameterSet&);
        ~PatMuonBooster();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag muonTag_;
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
PatMuonBooster::PatMuonBooster(const edm::ParameterSet& iConfig) :
    muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
    vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag"))
{
    produces<pat::MuonCollection>();  
}


PatMuonBooster::~PatMuonBooster()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
    void
PatMuonBooster::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    ESHandle<TransientTrackBuilder> theTTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);


    Handle<pat::MuonCollection> muons;
    iEvent.getByLabel(muonTag_,muons);

    Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(vertexTag_,vertices);

    //   if(vertices->size()==0) return;
    //   reco::Vertex vertex = vertices->front();

    reco::Vertex vertex;
    if(vertices->size()==0) {
//         vertex.SetXYZ(0,0,0);
        vertex = reco::Vertex(reco::Vertex::Point(0,0,0), reco::Vertex::Error());
    } else {
        vertex = vertices->front();
    }



    //Use the ExampleData to create an ExampleData2 which 
    // is put into the Event
    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);
    for(pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
        pat::Muon clone = *mu;
        if(mu->type() != 8){
            //        reco::TransientTrack tt = theTTBuilder->build(mu->track());
            //        Measurement1D ip = IPTools::absoluteTransverseImpactParameter(tt,vertex).second;
            //        Measurement1D ip3D = IPTools::absoluteImpactParameter3D(tt,vertex).second;

            //        clone.addUserFloat(std::string("tip"),ip.value());
            //        clone.addUserFloat(std::string("tipErr"),ip.error());
            //        clone.addUserFloat(std::string("ip"),ip3D.value());
            //        clone.addUserFloat(std::string("ipErr"),ip3D.error());


            //        try{
            // 	 GsfVertexTrackCompatibilityEstimator vertexTrackCompEst;       
            // 	 GsfVertexTrackCompatibilityEstimator::BDpair result = vertexTrackCompEst.estimate(vertex,tt);
            // 	 if(result.first)
            // 	   clone.addUserFloat(std::string("vtxComp"),result.second);
            // 	 else
            // 	   clone.addUserFloat(std::string("vtxComp"),99999.);
            //        }catch(...){
            // 	 cout << "exception caught" << endl;
            // 	 clone.addUserFloat(std::string("vtxComp"),99999.);
            //        }

            clone.addUserFloat( "dxyPV",clone.track()->dxy(vertices->front().position()) );
            clone.addUserFloat( "dzPV",clone.track()->dz(vertices->front().position()) );

        } else {
            clone.addUserFloat( "dxyPV", 9999);
            clone.addUserFloat( "dzPV", 9999);
        }
        pOut->push_back(clone);
    }


    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
    void 
PatMuonBooster::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PatMuonBooster::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatMuonBooster);
