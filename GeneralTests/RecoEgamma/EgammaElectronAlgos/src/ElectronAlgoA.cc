
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronAlgoA.h"
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


#include <Math/VectorUtil.h>
#include <Math/Point3D.h>

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include <TMath.h>
#include <sstream>

#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"

using namespace edm;
using namespace std;
using namespace reco;
//using namespace math; // conflicts with DataFormat/Math/interface/Point3D.h!!!!

ElectronAlgoA::ElectronAlgoA(double maxEOverPBarrel, double maxEOverPEndcaps, 
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

ElectronAlgoA::~ElectronAlgoA() {
  delete geomPropBw_;
  delete geomPropFw_;
  delete mtsTransform_;
}

void ElectronAlgoA::setupES(const edm::EventSetup& es, const edm::ParameterSet &conf) {

  //services
  es.get<TrackerRecoGeometryRecord>().get( theGeomSearchTracker );
  es.get<IdealMagneticFieldRecord>().get(theMagField);
  es.get<TrackerDigiGeometryRecord>().get(trackerHandle_);

  // get calo geometry
  es.get<IdealGeometryRecord>().get(theCaloGeom);
  
  mtsTransform_ = new MultiTrajectoryStateTransform;
  //geomPropBw_ = new GsfPropagatorAdapter(AnalyticalPropagator(theMagField.product(), oppositeToMomentum));
  //geomPropFw_ = new GsfPropagatorAdapter(AnalyticalPropagator(theMagField.product(), alongMomentum));

  geomPropBw_ = new AnalyticalPropagator(theMagField.product(), oppositeToMomentum);
  geomPropFw_ = new AnalyticalPropagator(theMagField.product(), alongMomentum);

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

void  ElectronAlgoA::run(Event& e, GlobalCtfElectronCollection & outEle) {

  // get the input 
  //BM edm::Handle<GsfTrackCollection> tracksBarrelH;
  //BM edm::Handle<GsfTrackCollection> tracksEndcapH;

  // to check existance
  edm::Handle<HBHERecHitCollection> hbhe;
  HBHERecHitMetaCollection *mhbhe=0;
  if (hOverEConeSize_ > 0.) {
    e.getByLabel(hbheLabel_,hbheInstanceName_,hbhe);  
    mhbhe=  new HBHERecHitMetaCollection(*hbhe);
  }
  //?????
  
  // get the tracks from the Event
  //e.getByLabel(trackBarrelLabel_,trackBarrelInstanceName_,tracksBarrelH);

  edm::Handle<reco::TrackCollection> tracksH; 
  e.getByLabel("ctfWithMaterialTracks",tracksH);

  edm::Handle<SuperClusterCollection> superClustersH; 
  e.getByLabel("correctedHybridSuperClusters",superClustersH);


  // get the Seed-SuperCluster association map
  /*
  edm::Handle<SeedSuperClusterAssociationCollection> barrelH;
  edm::Handle<SeedSuperClusterAssociationCollection> endcapH;
  e.getByLabel(assBarrelLabel_,assBarrelInstanceName_,barrelH);
  e.getByLabel(assEndcapLabel_,assEndcapInstanceName_,endcapH);
  */

  // get the GsfTrack-Seed association map
  /*
  edm::Handle<GsfTrackSeedAssociationCollection> barrelTSAssocH;
  edm::Handle<GsfTrackSeedAssociationCollection> endcapTSAssocH;
  e.getByLabel(assBarrelTrTSLabel_,assBarrelTrTSInstanceName_,barrelTSAssocH);
  e.getByLabel(assEndcapTrTSLabel_,assEndcapTrTSInstanceName_,endcapTSAssocH);
  edm::LogInfo("") 
    <<"\n\n Treating "<<e.id()<<", Number of seeds Barrel:"
    <<barrelH.product()->size()<<" Number of seeds Endcap:"<<endcapH.product()->size();
  */
  
  // create electrons from tracks in 2 steps: barrel + endcap
  //const SeedSuperClusterAssociationCollection  *sclAss=&(*barrelH);
  
  /*
  process(tracksBarrelH, //trackcollection
	  sclAss,   //seed-cluster assMap
	  barrelTSAssocH.product(), //track-seed assMap
	  mhbhe,  //calo rechit collection. what is the purpose??
	  outEle);

  sclAss=&(*endcapH);
  process(tracksEndcapH,sclAss,endcapTSAssocH.product(),mhbhe,outEle);
  delete mhbhe;
  */

  process(tracksH, //trackcollection
	  superClustersH,
	  mhbhe,  //calo rechit collection. what is the purpose??
	  outEle);


  /*
  std::ostringstream str;
  str << "========== ElectronAlgoA Info ==========";
  str << "Event " << e.id();
  str << "Number of final electron tracks: " << tracksBarrelH.product()->size()+ tracksEndcapH.product()->size();
  str << "Number of final electrons: " << outEle.size();
  for (vector<PixelMatchGsfElectron>::const_iterator it = outEle.begin(); it != outEle.end(); it++) {
    str << "New electron with charge, pt, eta, phi : "  << it->charge() << " , " 
        << it->pt() << " , " << it->eta() << " , " << it->phi();
  }
  str << "=================================================";
  LogDebug("ElectronAlgoA") << str.str();
  */
  return;
}

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

void ElectronAlgoA::process(edm::Handle<reco::TrackCollection> tracksH,
			   edm::Handle<reco::SuperClusterCollection> superClustersH,
			   HBHERecHitMetaCollection *mhbhe,
			   GlobalCtfElectronCollection & outEle) {
  
  //const GsfTrackCollection *tracks=tracksH.product();
  //const TrackCollection* tracks=tracksH.product();

  std::cout << "------- processing event" << std::endl;

  for (unsigned int i=0;i<tracksH->size();++i) {
    /*
    const GsfTrack & t=(*tracks)[i];
    const GsfTrackRef trackRef = edm::Ref<GsfTrackCollection>(tracksH,i);
    edm::Ref<TrajectorySeedCollection> seed = (*tsAss)[trackRef];
    const SuperCluster theClus=*((*sclAss)[seed]);
    */

    reco::TrackRef track(tracksH,i);
    reco::SuperClusterRef theCluster = superClusterMatching(track,superClustersH);
    if(theCluster.isNull()) continue;


    // calculate HoE
    double HoE;
    if (mhbhe) {
      CaloConeSelector sel(hOverEConeSize_, theCaloGeom.product(), DetId::Hcal);
      GlobalPoint pclu(theCluster->x(),theCluster->y(),theCluster->z());
      double hcalEnergy = 0.;
      std::auto_ptr<CaloRecHitMetaCollectionV> chosen=sel.select(pclu,*mhbhe);
      for (CaloRecHitMetaCollectionV::const_iterator i=chosen->begin(); i!=chosen->end(); i++) {
	//std::cout << HcalDetId(i->detid()) << " : " << (*i) << std::endl;
	hcalEnergy += i->energy();
      }
      HoE = hcalEnergy/theCluster->energy();
      LogDebug("") << "H/E : " << HoE;
    } else HoE=0;


    //reco::TransientTrack theTT(track,theMagField.product() );

    // calculate Trajectory StatesOnSurface....
    //at innermost point
    //GSF-specific TrajectoryStateOnSurface innTSOS = mtsTransform_->innerStateOnSurface(t, *(trackerHandle_.product()), theMagField.product());
    TrajectoryStateTransform transformer;
    TrajectoryStateOnSurface innTSOS  = transformer.innerStateOnSurface(*track,*trackerHandle_,theMagField.product());

    if (!innTSOS.isValid()) continue;

    //at vertex
    // innermost state propagation to the nominal vertex
    //GSF-specific TrajectoryStateOnSurface vtxTSOS =
    //TransverseImpactPointExtrapolator(*geomPropBw_).extrapolate(innTSOS,GlobalPoint(0,0,0));
    
    //TrajectoryStateOnSurface vtxTSOS = theTT.impactPointState(); TO FIX!!!!!!!!
    TrajectoryStateOnSurface vtxTSOS = innTSOS;
    if (!vtxTSOS.isValid()) vtxTSOS=innTSOS;

    //at seed
    //GSF-specific TrajectoryStateOnSurface outTSOS = mtsTransform_->outerStateOnSurface(t, *(trackerHandle_.product()), theMagField.product());
    //    TrajectoryStateOnSurface outTSOS = theTT.outermostMeasurementState();
    TrajectoryStateOnSurface outTSOS  = transformer.outerStateOnSurface(*track,*trackerHandle_,theMagField.product());
    if (!outTSOS.isValid()) continue;
    
    TrajectoryStateOnSurface seedTSOS = 
      TransverseImpactPointExtrapolator(*geomPropFw_).extrapolate(outTSOS,GlobalPoint(theCluster->seed()->position().x(),
										      theCluster->seed()->position().y(),
										      theCluster->seed()->position().z()));
    if (!seedTSOS.isValid()) seedTSOS=outTSOS;

    //at scl
    TrajectoryStateOnSurface sclTSOS = 
      TransverseImpactPointExtrapolator(*geomPropFw_).extrapolate(innTSOS,GlobalPoint(theCluster->x(),
										      theCluster->y(),
										      theCluster->z()));
    if (!sclTSOS.isValid()) sclTSOS=outTSOS;

    //GSF-specific GlobalVector vtxMom=computeMode(vtxTSOS);
    GlobalVector vtxMom = vtxTSOS.globalMomentum();
    const GlobalPoint  sclPos=sclTSOS.globalPosition();

    if (preSelection(*theCluster,vtxMom, sclPos, HoE)) {
      //GlobalVector innMom=computeMode(innTSOS);
      GlobalVector innMom= innTSOS.globalMomentum();
      GlobalPoint innPos=innTSOS.globalPosition();
      //GlobalVector seedMom=computeMode(seedTSOS);
      GlobalVector seedMom= seedTSOS.globalMomentum();
      GlobalPoint  seedPos=seedTSOS.globalPosition();
      //GlobalVector sclMom=computeMode(sclTSOS);    
      GlobalVector sclMom= sclTSOS.globalMomentum();    
      GlobalPoint  vtxPos=vtxTSOS.globalPosition();
      //GlobalVector outMom=computeMode(outTSOS);
      GlobalVector outMom= outTSOS.globalMomentum();
      GlobalPoint  outPos=outTSOS.globalPosition();

      GlobalCtfElectron ele(theCluster,
			    track,
			    sclPos,sclMom,
			    seedPos,seedMom,
			    innPos,innMom,
			    vtxPos,vtxMom,
			    outPos,outMom,
			    HoE);

      // set corrections + classification
      //ElectronClassification theClassifier;
      //theClassifier.correct(ele);
      //ElectronEnergyCorrector theEnCorrector;
      //theEnCorrector.correct(ele);
      //ElectronMomentumCorrector theMomCorrector;
      //theMomCorrector.correct(ele,vtxTSOS);
	//mCorr.getBestMomentum(),mCorr.getSCEnergyError(),mCorr.getTrackMomentumError());
      outEle.push_back(ele);
      //LogInfo("")<<"Constructed new electron with energy  "<< (*sclAss)[seed]->energy();
    }
  }  // loop over tracks
}

bool ElectronAlgoA::preSelection(const reco::SuperCluster& clus, const GlobalVector& tsosVtxMom, const GlobalPoint& tsosSclPos, double HoE) 
{
  LogDebug("")<< "========== preSelection ==========";
 
  LogDebug("") << "E/p : " << clus.energy()/tsosVtxMom.mag();
  if (tsosVtxMom.perp()<ptCut_)   return false;
  //FIXME: how to get detId from a cluster??
  std::vector<DetId> vecId=clus.getHitsByDetId();
  int subdet =vecId[0].subdetId();  //FIXME: is the first one really the biggest??
  if ((subdet==EcalBarrel) && (clus.energy()/tsosVtxMom.mag() > maxEOverPBarrel_)) return false;
  if ((subdet==EcalEndcap) && (clus.energy()/tsosVtxMom.mag() > maxEOverPEndcaps_)) return false;
  LogDebug("") << "E/p criteria is satisfied ";
  // delta eta criteria
  double etaclu = clus.eta();
  double etatrk = tsosSclPos.eta();
  double deta = etaclu-etatrk;
  LogDebug("") << "delta eta : " << deta;
  if (fabs(deta) > maxDeltaEta_) return false;
  LogDebug("") << "Delta eta criteria is satisfied ";
  // delta phi criteria
  double phiclu = clus.phi();
  double phitrk = tsosSclPos.phi();
  double dphi = phiclu-phitrk;
  LogDebug("") << "delta phi : " << dphi;
  if (fabs(dphi) > maxDeltaPhi_) return false;
  LogDebug("") << "Delta phi criteria is satisfied ";

  if (HoE > maxHOverE_) return false; //FIXME: passe dans tous les cas?
  LogDebug("") << "H/E criteria is satisfied ";

  LogDebug("") << "electron has passed preselection criteria ";
  LogDebug("") << "=================================================";
  return true;  
}  

GlobalVector ElectronAlgoA::computeMode(const TrajectoryStateOnSurface &tsos) {
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
	 
    if ( myGSUtil_Px ) { delete myGSUtil_Px; }
    if ( myGSUtil_Py ) { delete myGSUtil_Py; }
    if ( myGSUtil_Pz ) { delete myGSUtil_Pz; }
	  
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

const reco::SuperClusterRef 
ElectronAlgoA::superClusterMatching(reco::TrackRef track,edm::Handle<reco::SuperClusterCollection> collection) 

{
  double minDr = 10000.0;
  reco::SuperClusterRef theCluster = edm::Ref<SuperClusterCollection>();
  math::XYZVector trackGlobalDir(track->momentum());

  for(reco::SuperClusterCollection::size_type i=0; i<collection->size(); ++i){
    reco::SuperClusterRef cluster(collection, i);
    math::XYZVector clusterGlobalPos(cluster->x(),cluster->y(),cluster->z());

    double tmpDr = ROOT::Math::VectorUtil::DeltaR(clusterGlobalPos, trackGlobalDir);
    if(tmpDr < minDr){
      minDr = tmpDr;
      theCluster = cluster;
    }
  }
  //std::cout << "returning null ref" << std::endl;
  return theCluster;

}
