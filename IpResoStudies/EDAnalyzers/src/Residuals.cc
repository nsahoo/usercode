// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
 
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


#include <TrackingTools/TrajectoryState/interface/PerigeeConversions.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>

#include "IpResoStudies/EDAnalyzers/interface/VertexReProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "TH1F.h"
#include "TTree.h"

// structure to hold the tree raw data

struct treeRaw{
  double pt;
  double p;
  double eta;
  double phi;
  int nXLayers;
  int nMissedOut;
  int nMissedIn;
  int hasPXL;
  int    quality;
  double d0;
  double dz;
  double d0Err;
  double dzErr;
};
  
  
//
// class declaration
//
  
  
class Residuals : public edm::EDAnalyzer {
public:
  explicit Residuals(const edm::ParameterSet& pset);
  ~Residuals();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
    
  bool trackSelection(const reco::Track& track) const;
  bool vertexSelection(const reco::Vertex& vertex) const;
  
  // ----------member data ---------------------------
  edm::InputTag trackLabel;
  edm::InputTag vertexLabel;

  // --- track selection variables
  double tkMinPt;
  int tkMinXLayers,tkMaxMissedOuterLayers,tkMaxMissedInnerLayers;

  // --- vertex selection variables
  unsigned int vtxTracksSizeMin;  
  unsigned int vtxTracksSizeMax;  
  double vtxErrorXMin,vtxErrorXMax;
  double vtxErrorYMin,vtxErrorYMax;
  double vtxErrorZMin,vtxErrorZMax;

  TH1F *h_d0;
  TTree *tree;
  treeRaw raw;
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
Residuals::Residuals(const edm::ParameterSet& pset){
  //about configuration
  trackLabel  = pset.getParameter<edm::InputTag>("TrackLabel");    
  vertexLabel = pset.getParameter<edm::InputTag>("VertexLabel");    

  tkMinPt = pset.getParameter<double>("TkMinPt");    
  tkMinXLayers = pset.getParameter<int>("TkMinXLayers");
  tkMaxMissedOuterLayers = pset.getParameter<int>("TkMaxMissedOuterLayers");
  tkMaxMissedInnerLayers = pset.getParameter<int>("TkMaxMissedInnerLayers");

  vtxTracksSizeMin = pset.getParameter<int>("VtxTracksSizeMin");
  vtxTracksSizeMax = pset.getParameter<int>("VtxTracksSizeMax");
  vtxErrorXMin     = pset.getParameter<double>("VtxErrorXMin");
  vtxErrorXMax     = pset.getParameter<double>("VtxErrorXMax");
  vtxErrorYMin     = pset.getParameter<double>("VtxErrorYMin");
  vtxErrorYMax     = pset.getParameter<double>("VtxErrorYMax");
  vtxErrorZMin     = pset.getParameter<double>("VtxErrorZMin");
  vtxErrorZMax     = pset.getParameter<double>("VtxErrorZMax");
    

   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "tree"  , "recoTrack IP residuals");
  tree->Branch("raw",&raw.pt,"pt/D:p/D:eta/D:phi/D:nXLayers/I:nMissedOut/I:nMissedIn/I:hasPXL/I:quality/I:d0/D:dz:d0Err:dzErr");
}


Residuals::~Residuals()
{

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
Residuals::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

   Handle<VertexCollection> vtxH;
   iEvent.getByLabel(vertexLabel, vtxH);

   ESHandle<MagneticField> theMF;
   iSetup.get<IdealMagneticFieldRecord>().get(theMF);


   VertexReProducer revertex(vtxH, iEvent);
   Handle<TrackCollection> pvtracks;   iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
   Handle<BeamSpot>        pvbeamspot; iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);


   Handle<TrackCollection> tracks;  iEvent.getByLabel(trackLabel, tracks);
   if(tracks.id() != pvtracks.id())
     cout << "WARNING: the tracks originally used for PV are not the same used in this analyzer." 
	  << "Is this really what you want?" << endl;

   //if (pvbeamspot.id() != theBeamSpot.id()) 
   //  edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";


   // ------ check if the vertex is good enough -------
   if(vtxH->size()==0) return;
   if(! vertexSelection(vtxH->front()) ) return;
   // -------------------------------------------------

   /*
   cout << "original vtx x,y,z: " 
	<< vtxH->front().position().x() << " , "
	<< vtxH->front().position().y() << " , "
	<< vtxH->front().position().z() << endl;
   */

   for(TrackCollection::const_iterator itk = tracks->begin(); itk != tracks->end(); ++itk){

     // --- track selection ---
     if(! trackSelection(*itk)) continue;
     // ---
     

     TrackCollection newTkCollection;
     newTkCollection.assign(tracks->begin(), itk);
     newTkCollection.insert(newTkCollection.end(),itk+1,tracks->end());
     //newTkCollection.insert(newTkCollection.end(),itk,tracks->end()); // only for debugging purpose


     //cout << "before,after size: " << tkH->size() << " , " << newCollection.size() << endl;
  
    

     // --- from Giovanni to refit the prim.vertex     
     vector<TransientVertex> pvs = revertex.makeVertices(newTkCollection, *pvbeamspot, iSetup) ;
     if (pvs.empty()) continue;
     reco::Vertex newPV = reco::Vertex(pvs.front());
     Track::Point vtxPosition = Track::Point(newPV.position().x(),
					     newPV.position().y(),
					     newPV.position().z());
     // ---          
     if(! vertexSelection(newPV) ) continue;





     
     double d0 = itk->dxy(vtxPosition);
     double dz = itk->dz(vtxPosition);

     //Filling the tree
     raw.pt  = itk->pt();
     raw.p   = itk->p();
     raw.eta = itk->eta();
     raw.phi = itk->phi();
     raw.nXLayers   = itk->hitPattern().trackerLayersWithMeasurement();
     raw.nMissedOut = itk->trackerExpectedHitsOuter().numberOfLostHits();
     raw.nMissedIn  = itk->trackerExpectedHitsInner().numberOfLostHits();
     raw.hasPXL     = (itk->hitPattern().hasValidHitInFirstPixelBarrel() || 
		       itk->hitPattern().hasValidHitInFirstPixelEndcap());
     raw.quality = itk->qualityMask();
     raw.d0 = d0*10000.;
     raw.dz = dz*10000.;
     raw.d0Err = itk->d0Error()*10000;
     raw.dzErr = itk->dzError()*10000;

     tree->Fill();

   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
Residuals::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Residuals::endJob() {
}


bool
Residuals::trackSelection(const reco::Track& track) const {
  if( track.pt() < tkMinPt) return false;
  if( track.hitPattern().trackerLayersWithMeasurement() < tkMinXLayers) return false;
  if( track.trackerExpectedHitsOuter().numberOfLostHits() > tkMaxMissedOuterLayers) return false;
  if( track.trackerExpectedHitsInner().numberOfLostHits() > tkMaxMissedInnerLayers) return false;

  if( ! track.quality(reco::TrackBase::highPurity) ) return false;
  //if(!track.hitPattern().hasValidHitInFirstPixelBarrel()) return false;
  return true;
}

bool
Residuals::vertexSelection(const reco::Vertex& vertex) const{
  if(vertex.tracksSize()>vtxTracksSizeMax || vertex.tracksSize()<vtxTracksSizeMin) return false;
  if(vertex.xError() < vtxErrorXMin || vertex.xError() > vtxErrorXMax) return false;
  if(vertex.yError() < vtxErrorYMin || vertex.yError() > vtxErrorYMax) return false;
  if(vertex.zError() < vtxErrorZMin || vertex.zError() > vtxErrorZMax) return false;
  return true;
}



//define this as a plug-in
DEFINE_FWK_MODULE(Residuals);

