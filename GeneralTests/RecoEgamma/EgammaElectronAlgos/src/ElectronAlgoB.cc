
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronAlgoB.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronClassification.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronMomentumCorrector.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronEnergyCorrector.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"


#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/GSUtilities.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"

#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include <TMath.h>
#include <sstream>

#include <Math/VectorUtil.h>
#include <Math/Point3D.h>

#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"

using namespace edm;
using namespace std;
using namespace reco;
//using namespace math; // conflicts with DataFormat/Math/interface/Point3D.h!!!!

ElectronAlgoB::ElectronAlgoB(double maxEOverPBarrel, double maxEOverPEndcaps, 
			     double hOverEConeSize, double maxHOverE, 
			     double maxDeltaEta, double maxDeltaPhi, double ptcut):  
  maxEOverPBarrel_(maxEOverPBarrel), maxEOverPEndcaps_(maxEOverPEndcaps), 
  hOverEConeSize_(hOverEConeSize), maxHOverE_(maxHOverE), 
  maxDeltaEta_(maxDeltaEta), maxDeltaPhi_(maxDeltaPhi), ptCut_(ptcut)
{   
  geomPropBw_=0;	
  geomPropFw_=0;	
  mtsTransform_=0;
}

ElectronAlgoB::~ElectronAlgoB() {
  delete geomPropBw_;
  delete geomPropFw_;
  delete mtsTransform_;
}

void ElectronAlgoB::setupES(const edm::EventSetup& es, const edm::ParameterSet &conf) {

  //services
  es.get<TrackerRecoGeometryRecord>().get( theGeomSearchTracker );
  es.get<IdealMagneticFieldRecord>().get(theMagField);
  es.get<TrackerDigiGeometryRecord>().get(trackerHandle_);

  // get calo geometry
  es.get<IdealGeometryRecord>().get(theCaloGeom);
  
  mtsTransform_ = new MultiTrajectoryStateTransform;
  geomPropBw_ = new GsfPropagatorAdapter(AnalyticalPropagator(theMagField.product(), oppositeToMomentum));
  geomPropFw_ = new GsfPropagatorAdapter(AnalyticalPropagator(theMagField.product(), alongMomentum));

  // get nested parameter set for the TransientInitialStateEstimator
  ParameterSet tise_params = conf.getParameter<ParameterSet>("TransientInitialStateEstimatorParameters") ;
  hbheLabel_ = conf.getParameter<string>("hbheModule");
  hbheInstanceName_ = conf.getParameter<string>("hbheInstance");
  trackBarrelLabel_ = conf.getParameter<string>("TrackBarrelLabel");
  trackBarrelInstanceName_ = conf.getParameter<string>("TrackBarrelProducer");
  trackEndcapLabel_ = conf.getParameter<string>("TrackEndcapLabel");
  trackEndcapInstanceName_ = conf.getParameter<string>("TrackEndcapProducer");
  assBarrelLabel_ = conf.getParameter<string>("SCLBarrelLabel");
  assBarrelInstanceName_ = conf.getParameter<string>("SCLBarrelProducer");
  assEndcapLabel_ = conf.getParameter<string>("SCLEndcapLabel");
  assEndcapInstanceName_ = conf.getParameter<string>("SCLEndcapProducer");
  assBarrelTrTSLabel_ = conf.getParameter<string>("AssocTrTSBarrelLabel");
  assBarrelTrTSInstanceName_ = conf.getParameter<string>("AssocTrTBarrelProducer");
  assEndcapTrTSLabel_ = conf.getParameter<string>("AssocTrTEndcapLabel");
  assEndcapTrTSInstanceName_ = conf.getParameter<string>("AssocTrTEndcapProducer");
}

void  ElectronAlgoB::run(Event& e, PixelMatchGsfElectronCollection & outEle) {

  // get the input 
  edm::Handle<GsfTrackCollection> tracksBarrelH;
  edm::Handle<GsfTrackCollection> tracksEndcapH;
  // to check existance
  edm::Handle<HBHERecHitCollection> hbhe;
  HBHERecHitMetaCollection *mhbhe=0;
  if (hOverEConeSize_ > 0.) {
    e.getByLabel(hbheLabel_,hbheInstanceName_,hbhe);  
    mhbhe=  new HBHERecHitMetaCollection(*hbhe);
  }

  edm::Handle<reco::GsfTrackCollection> tracksH; 
  e.getByLabel("GsfElectrons",tracksH);

  edm::Handle<SuperClusterCollection> superClustersBarrelH; 
  e.getByLabel("correctedHybridSuperClusters",superClustersBarrelH);
  
  edm::Handle<SuperClusterCollection> superClustersEndcapH; 
  e.getByLabel("islandSuperClusters", "islandEndcapSuperClusters", superClustersEndcapH);

  process(tracksH, //trackcollection
          superClustersBarrelH, 
          superClustersEndcapH,      
          mhbhe,  //calo rechit collection. what is the purpose??
          outEle);

  return;
}

void ElectronAlgoB::process(edm::Handle<GsfTrackCollection> tracksH,
                            edm::Handle<reco::SuperClusterCollection> superClustersBarrelH,
                            edm::Handle<reco::SuperClusterCollection> superClustersEndcapH,
                            HBHERecHitMetaCollection *mhbhe,
                            PixelMatchGsfElectronCollection & outEle) {
  
  std::cout << "------- processing event" << std::endl;
  
  if (tracksH->size() == 0) {
    std::cout << "Electron lost: no track found. " << std::endl;
  } else {
    std::cout << "Number of tracks: " << tracksH->size() << std::endl;
  }
  
  std::cout << "SuperCluster: " << superClustersBarrelH->size() << "  " << 
    superClustersEndcapH->size() << std::endl;	

  reco::SuperClusterRefVector superClusters;

  for(int z=0; z<2; ++z) {

    superClusters.clear();
    if (z == 0) {
      for(reco::SuperClusterCollection::size_type i= 0; i<superClustersBarrelH->size(); ++i){
        reco::SuperClusterRef cluster(superClustersBarrelH, i);
        superClusters.push_back(cluster);
      }
      std::cout << superClustersBarrelH->size() << std::endl;
    }
    
    if (z == 1) {
      for(reco::SuperClusterCollection::size_type i= 0; i<superClustersEndcapH->size(); ++i){
        reco::SuperClusterRef cluster(superClustersEndcapH, i);
        superClusters.push_back(cluster);
      }
      std::cout << superClustersEndcapH->size() << std::endl;
    }
    
    //================= loop over SuperClusters ===============

    for(unsigned int i=0; i< superClusters.size(); ++i) {

      std::cout << "Start matching " << std::endl;	
      reco::SuperClusterRef theClus = superClusters[i];
      reco::GsfTrackRef track = superClusterMatching(theClus, tracksH);
      
      if(track.isNull()) {
        std::cout << "Electron lost: no supercluster match found: " << tracksH->size() << std::endl;
        continue;
      }
      std::cout << "End matching " << std::endl;
      
      
      // calculate HoE
      std::cout << "Start HoE " << std:: endl;	
      // calculate HoE
      double HoE;
      if (mhbhe) {
        CaloConeSelector sel(hOverEConeSize_, theCaloGeom.product(), DetId::Hcal);
        GlobalPoint pclu(theClus->x(),theClus->y(),theClus->z());
        double hcalEnergy = 0.;
        std::auto_ptr<CaloRecHitMetaCollectionV> chosen=sel.select(pclu,*mhbhe);
        for (CaloRecHitMetaCollectionV::const_iterator i=chosen->begin(); i!=chosen->end(); i++) {
          //std::cout << HcalDetId(i->detid()) << " : " << (*i) << std::endl;
          hcalEnergy += i->energy();
        }
        HoE = hcalEnergy/theClus->energy();
        LogDebug("") << "H/E : " << HoE;
      } else HoE=0;
      
      
      // calculate Trajectory StatesOnSurface....
      //at innermost point
      TrajectoryStateOnSurface innTSOS = mtsTransform_->innerStateOnSurface(*track, *(trackerHandle_.product()), theMagField.product());
      if (!innTSOS.isValid()) 
        continue;

      //at vertex
      // innermost state propagation to the nominal vertex
      TrajectoryStateOnSurface vtxTSOS =
        TransverseImpactPointExtrapolator(*geomPropBw_).extrapolate(innTSOS,GlobalPoint(0,0,0));
      if (!vtxTSOS.isValid()) 
        vtxTSOS=innTSOS;
      
      //at seed
      TrajectoryStateOnSurface outTSOS = mtsTransform_->outerStateOnSurface(*track, *(trackerHandle_.product()), theMagField.product());
      if (!outTSOS.isValid()) 
        continue;
      
      TrajectoryStateOnSurface seedTSOS = TransverseImpactPointExtrapolator(*geomPropFw_).extrapolate(outTSOS,GlobalPoint(theClus->seed()->position().x(),theClus->seed()->position().y(),theClus->seed()->position().z()));
      if (!seedTSOS.isValid()) 
        seedTSOS=outTSOS;
      
      //at scl
      TrajectoryStateOnSurface sclTSOS = TransverseImpactPointExtrapolator(*geomPropFw_).extrapolate(innTSOS,GlobalPoint(theClus->x(),theClus->y(),theClus->z()));
      if (!sclTSOS.isValid()) 
        sclTSOS=outTSOS;
      
      GlobalVector vtxMom=computeMode(vtxTSOS);
      GlobalPoint  sclPos=sclTSOS.globalPosition();

      if (preSelection(*theClus,vtxMom, sclPos, HoE)) {
        GlobalVector innMom=computeMode(innTSOS);
        GlobalPoint innPos=innTSOS.globalPosition();
        GlobalVector seedMom=computeMode(seedTSOS);
        GlobalPoint  seedPos=seedTSOS.globalPosition();
        GlobalVector sclMom=computeMode(sclTSOS);    
        GlobalPoint  vtxPos=vtxTSOS.globalPosition();
        GlobalVector outMom=computeMode(outTSOS);
        GlobalPoint  outPos=outTSOS.globalPosition();
        
        PixelMatchGsfElectron ele(theClus,
                                  track,
                                  sclPos,sclMom,
                                  seedPos,seedMom,
                                  innPos,innMom,
                                  vtxPos,vtxMom,
                                  outPos,outMom,
                                  HoE);
        /*
          GlobalGsfElectron ele(theClus,
          track,
          sclPos,sclMom,
          seedPos,seedMom,
          innPos,innMom,
          vtxPos,vtxMom,
          outPos,outMom,
          HoE);
        */
        
        // set corrections + classification
        ElectronClassification theClassifier;
        theClassifier.correct(ele);  
        ElectronEnergyCorrector theEnCorrector;
        theEnCorrector.correct(ele); 
        ElectronMomentumCorrector theMomCorrector;
        theMomCorrector.correct(ele,vtxTSOS); 
        //mCorr.getBestMomentum(),mCorr.getSCEnergyError(),mCorr.getTrackMomentumError());
        outEle.push_back(ele);
        //LogInfo("")<<"Constructed new electron with energy  "<< (*sclAss)[seed]->energy();
      }
    }  
  }
}

bool ElectronAlgoB::preSelection(const SuperCluster& clus, const GlobalVector& tsosVtxMom, const GlobalPoint& tsosSclPos, double HoE) 
{
  LogDebug("")<< "========== preSelection ==========";
 
  LogDebug("") << "E/p : " << clus.energy()/tsosVtxMom.mag();
  if (tsosVtxMom.perp()<ptCut_)   
    return false;
  
  //FIXME: how to get detId from a cluster??
  
  std::vector<DetId> vecId=clus.getHitsByDetId();
  int subdet =vecId[0].subdetId();  //FIXME: is the first one really the biggest??
  if ((subdet==EcalBarrel) && (clus.energy()/tsosVtxMom.mag() > maxEOverPBarrel_)) 
    return false;
  if ((subdet==EcalEndcap) && (clus.energy()/tsosVtxMom.mag() > maxEOverPEndcaps_)) 
    return false;
  LogDebug("") << "E/p criteria is satisfied ";
  // delta eta criteria
  double etaclu = clus.eta();
  double etatrk = tsosSclPos.eta();
  double deta = etaclu-etatrk;
  LogDebug("") << "delta eta : " << deta;
  if (fabs(deta) > maxDeltaEta_) 
    return false;
  LogDebug("") << "Delta eta criteria is satisfied ";
  // delta phi criteria
  double phiclu = clus.phi();
  double phitrk = tsosSclPos.phi();
  double dphi = phiclu-phitrk;
  LogDebug("") << "delta phi : " << dphi;
  if (fabs(dphi) > maxDeltaPhi_) 
    return false;
  LogDebug("") << "Delta phi criteria is satisfied ";

  if (HoE > maxHOverE_) 
    return false; //FIXME: passe dans tous les cas?
  LogDebug("") << "H/E criteria is satisfied ";

  LogDebug("") << "electron has passed preselection criteria ";
  LogDebug("") << "=================================================";

  return true;  
}  

GlobalVector ElectronAlgoB::computeMode(const TrajectoryStateOnSurface &tsos) {

  // mode computation	
  float mode_Px = 0.;
  float mode_Py = 0.;
  float mode_Pz = 0.;

  if ( tsos.isValid() ){
	  
    int Count = 0;
    unsigned int numb = tsos.components().size();
    float *Wgt   = new float[numb];
    float *Px    = new float[numb];
    float *Py    = new float[numb];
    float *Pz    = new float[numb];
    float *PxErr = new float[numb];
    float *PyErr = new float[numb];
    float *PzErr = new float[numb];
	  
    for (unsigned int ii = 0; ii < numb; ii ++){
      Wgt[ii]   = 0.;
      Px[ii]    = 0.;
      Py[ii]    = 0.;
      Pz[ii]    = 0.;
      PxErr[ii] = 0.;
      PyErr[ii] = 0.;
      PzErr[ii] = 0.;
    }
	  
    std::vector<TrajectoryStateOnSurface> comp = tsos.components();
    for (std::vector<TrajectoryStateOnSurface>::const_iterator it_comp = comp.begin(); it_comp!= comp.end(); it_comp++){
      Wgt[Count]    = it_comp->weight();
      Px[Count]     = it_comp->globalMomentum().x();
      Py[Count]     = it_comp->globalMomentum().y();
      Pz[Count]     = it_comp->globalMomentum().z();
      PxErr[Count]  = sqrt(it_comp->cartesianError().matrix()[3][3]);
      PyErr[Count]  = sqrt(it_comp->cartesianError().matrix()[4][4]);
      PzErr[Count]  = sqrt(it_comp->cartesianError().matrix()[5][5]);
      Count++;
    }
	  
    GSUtilities *myGSUtil_Px = new GSUtilities(numb, Wgt, Px, PxErr);
    GSUtilities *myGSUtil_Py = new GSUtilities(numb, Wgt, Py, PyErr);
    GSUtilities *myGSUtil_Pz = new GSUtilities(numb, Wgt, Pz, PzErr);
	  
    mode_Px = myGSUtil_Px->mode();
    mode_Py = myGSUtil_Py->mode();
    mode_Pz = myGSUtil_Pz->mode();
	 
    if ( myGSUtil_Px ) 
      { delete myGSUtil_Px; }
    if ( myGSUtil_Py ) 
      { delete myGSUtil_Py; }
    if ( myGSUtil_Pz ) 
      { delete myGSUtil_Pz; }
	  
    delete[] Wgt;
    delete[] Px;
    delete[] PxErr;
    delete[] Py;
    delete[] PyErr;
    delete[] Pz;
    delete[] PzErr;
  } else edm::LogInfo("") << "tsos not valid!!";

  return GlobalVector(mode_Px,mode_Py,mode_Pz);	
}

const reco::GsfTrackRef
ElectronAlgoB::superClusterMatching(reco::SuperClusterRef sc, edm::Handle<reco::GsfTrackCollection> tracks) {

  double minDr = 0.5;
  //reco::SuperClusterRef theClus = edm::Ref<SuperClusterCollection>();
  reco::GsfTrackRef theTrack = edm::Ref<reco::GsfTrackCollection>();


  for(reco::GsfTrackCollection::size_type i=0; i<tracks->size(); ++i){
    reco::GsfTrackRef track(tracks, i);
    math::XYZVector trackGlobalDir(track->momentum());   
    math::XYZVector clusterGlobalPos(sc->x() - track->vx(), sc->y() - track->vy(), sc->z() - track->vz());
 
    double tmpDr = ROOT::Math::VectorUtil::DeltaR(clusterGlobalPos, trackGlobalDir);
    if(tmpDr < minDr){
      minDr = tmpDr;
      theTrack = track;
    }
  }

  //std::cout << "returning null ref" << std::endl;
  return theTrack;
}
