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
    ESHandle<TransientTrackBuilder> theTTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);


    Handle<pat::ElectronCollection> electrons;
    iEvent.getByLabel(electronTag_,electrons);

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
    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);
    for(pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){
        //      reco::TransientTrack tt = theTTBuilder->build(ele->gsfTrack());
        //cout << "here1" << endl;
        //      Measurement1D ip = IPTools::absoluteTransverseImpactParameter(tt,vertex).second;
        //      Measurement1D ip3D = IPTools::absoluteImpactParameter3D(tt,vertex).second;

        //cout << "here2" << endl;
        pat::Electron clone = *ele;
        //      clone.addUserFloat(std::string("tip"),ip.value());
        //      clone.addUserFloat(std::string("tipErr"),ip.error());
        //      clone.addUserFloat(std::string("ip"),ip3D.value());
        //      clone.addUserFloat(std::string("ipErr"),ip3D.error());


        //      try{
        //        GsfVertexTrackCompatibilityEstimator vertexTrackCompEst;       
        //        GsfVertexTrackCompatibilityEstimator::BDpair result = vertexTrackCompEst.estimate(vertex,tt);
        //        if(result.first)
        // 	 clone.addUserFloat(std::string("vtxComp"),result.second);
        //        else
        // 	 clone.addUserFloat(std::string("vtxComp"),99999.);
        //        //cout << "here4" << endl;
        //      }catch(...){
        //        cout << "exception caught" << endl;
        //        clone.addUserFloat(std::string("vtxComp"),99999.);
        //      }

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
