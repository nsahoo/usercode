#ifndef ElectronAlgoA_H
#define ElectronAlgoA_H


#include "DataFormats/EgammaCandidates/interface/GlobalCtfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SeedSuperClusterAssociation.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/PixelMatchTrackReco/interface/GsfTrackSeedAssociation.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

class MultiTrajectoryStateTransform;
class GsfPropagatorAdapter;
 
class ElectronAlgoA {

public:

  ElectronAlgoA(double maxEOverPBarrel, double maxEOverPBarrel, 
		double hOverEConeSize, double maxHOverE, 
		double maxDeltaEta, double maxDeltaPhi, double ptCut);
  
  ~ElectronAlgoA();

  void setupES(const edm::EventSetup& setup, const edm::ParameterSet& conf);
  void run(edm::Event&, reco::GlobalCtfElectronCollection&);

 private:

  // create electrons from tracks
  void process(edm::Handle<reco::TrackCollection> tracksH,
               edm::Handle<reco::SuperClusterCollection> superClustersBarrelH,
               edm::Handle<reco::SuperClusterCollection> superClustersEndcapH,
               HBHERecHitMetaCollection *mhbhe,
               reco::GlobalCtfElectronCollection & outEle);

  const reco::TrackRef superClusterMatching(reco::SuperClusterRef,
                                            edm::Handle<reco::TrackCollection>);


  // preselection method
  //  bool preSelection(const reco::SuperCluster& clus, const reco::GsfTrack& track,double HoE);
  bool preSelection(const reco::SuperCluster& clus, const GlobalVector&, const GlobalPoint&,double HoE);
  
  
  //Gsf mode calculations
  GlobalVector computeMode(const TrajectoryStateOnSurface &tsos);
  
  // preselection parameters
  // maximum E/p where E is the supercluster corrected energy and p the track momentum at innermost state  
  double maxEOverPBarrel_;   
  double maxEOverPEndcaps_;   
  // cone size for H/E
  double hOverEConeSize_; 
  // maximum H/E where H is the Hcal energy inside the cone centered on the seed cluster eta-phi position 
  double maxHOverE_; 
  // maximum eta difference between the supercluster position and the track position at the closest impact to the supercluster 
  double maxDeltaEta_;
  // maximum phi difference between the supercluster position and the track position at the closest impact to the supercluster
  // position to the supercluster
  double maxDeltaPhi_;

  double ptCut_;
 
  // input configuration
  std::string hbheLabel_;
  std::string hbheInstanceName_;
  std::string trackBarrelLabel_;
  std::string trackEndcapLabel_;
  std::string trackBarrelInstanceName_;
  std::string trackEndcapInstanceName_;
  std::string assBarrelLabel_;
  std::string assBarrelInstanceName_;
  std::string assEndcapLabel_;
  std::string assEndcapInstanceName_;
  std::string assBarrelTrTSLabel_;
  std::string assBarrelTrTSInstanceName_;
  std::string assEndcapTrTSLabel_;
  std::string assEndcapTrTSInstanceName_;

  edm::ESHandle<MagneticField>                theMagField;
  edm::ESHandle<GeometricSearchTracker>       theGeomSearchTracker;
  edm::ESHandle<CaloGeometry>                 theCaloGeom;
  edm::ESHandle<TrackerGeometry>              trackerHandle_;

  const MultiTrajectoryStateTransform *mtsTransform_;
  //const GsfPropagatorAdapter *geomPropBw_;
  //const GsfPropagatorAdapter *geomPropFw_;
  
  const Propagator* geomPropBw_;
  const Propagator* geomPropFw_;
};

#endif // ElectronAlgoA_H


