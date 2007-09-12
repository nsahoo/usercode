
#include "RecoEgammaBis/EgammaElectronAlgos/interface/ElectronAlgoB.h"
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
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"


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
			     double minEOverPBarrel, double minEOverPEndcaps,
			     double hOverEConeSize, double maxHOverE, 
			     double maxDeltaEta, double maxDeltaPhi, double ptcut,
			     bool highPtPresel, double highPtMin):  
  maxEOverPBarrel_(maxEOverPBarrel), maxEOverPEndcaps_(maxEOverPEndcaps), 
  minEOverPBarrel_(minEOverPBarrel), minEOverPEndcaps_(minEOverPEndcaps), 
  hOverEConeSize_(hOverEConeSize), maxHOverE_(maxHOverE), 
  maxDeltaEta_(maxDeltaEta), maxDeltaPhi_(maxDeltaPhi), ptCut_(ptcut),
  highPtPreselection_(highPtPresel), highPtMin_(highPtMin)
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
  e.getByLabel("correctedEndcapSuperClustersWithPreshower", superClustersEndcapH);

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
  
  //std::cout << "------- processing event" << std::endl;
  
  if (tracksH->size() == 0) {
    //std::cout << "Electron lost: no track found. " << std::endl;
  } else {
    //std::cout << "Number of tracks: " << tracksH->size() << std::endl;
  }
  
  //std::cout << "SuperCluster: " << superClustersBarrelH->size() << "  " << 
  //  superClustersEndcapH->size() << std::endl;	

  reco::SuperClusterRefVector superClusters;

  for(int z=0; z<2; ++z) {

    superClusters.clear();
    if (z == 0) {
      for(reco::SuperClusterCollection::size_type i= 0; i<superClustersBarrelH->size(); ++i){
        reco::SuperClusterRef cluster(superClustersBarrelH, i);
        superClusters.push_back(cluster);
      }
      //std::cout << superClustersBarrelH->size() << std::endl;
    }
    
    if (z == 1) {
      for(reco::SuperClusterCollection::size_type i= 0; i<superClustersEndcapH->size(); ++i){
        reco::SuperClusterRef cluster(superClustersEndcapH, i);
        superClusters.push_back(cluster);
      }
      //std::cout << superClustersEndcapH->size() << std::endl;
    }
    
    //================= loop over SuperClusters ===============

    for(unsigned int i=0; i< superClusters.size(); ++i) {

      //std::cout << "Start matching " << std::endl;	
      reco::SuperClusterRef theClus = superClusters[i];
      reco::GsfTrackRef track = superClusterMatching(theClus, tracksH);
      
      if(track.isNull()) {
        //std::cout << "Electron lost: no supercluster match found: " << tracksH->size() << std::endl;
        continue;
      }
      
      

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
        
	//create electron
	double scale = theClus->energy()/vtxMom.mag();    
	math::XYZTLorentzVectorD momentum= math::XYZTLorentzVector(vtxMom.x()*scale,
								   vtxMom.y()*scale,
								   vtxMom.z()*scale,
								   theClus->energy());
	PixelMatchGsfElectron ele(momentum,
				  theClus,
				  track,
				  sclPos,sclMom,
				  seedPos,seedMom,
				  innPos,innMom,
				  vtxPos,vtxMom,
				  outPos,outMom,
				  HoE);
        
	//and set various properties
	float trackEta = ecalEta(
				 track->innerMomentum().eta(),
				 track->innerPosition().z(),
				 track->innerPosition().Rho());
	
	float trackPhi = ecalPhi(
				 track->innerMomentum().Rho(),
				 track->innerMomentum().eta(),
				 track->innerMomentum().phi(),
				 track->charge(),
				 track->innerPosition().Rho());
	
	ele.setDeltaEtaSuperClusterAtVtx(theClus->position().eta() - trackEta);
	float dphi = theClus->position().phi() - trackPhi;
	if (fabs(dphi)>CLHEP::pi)
	  dphi = dphi < 0? CLHEP::pi2 + dphi : dphi - CLHEP::pi2;
	ele.setDeltaPhiSuperClusterAtVtx(dphi);
      
      
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

  // pt min
  LogDebug("") << "pT : " << tsosVtxMom.perp();
  if (tsosVtxMom.perp() < ptCut_)   return false;

  // E/p cut
  LogDebug("") << "E/p : " << clus.energy()/tsosVtxMom.mag();
  std::vector<DetId> vecId=clus.getHitsByDetId();
  int subdet =vecId[0].subdetId();  
  double rt2 = clus.x()*clus.x() + clus.y()*clus.y();
  double r2 = rt2 + clus.z()*clus.z();
  // no E/p preselection for high pT electrons
  if (!highPtPreselection_ || clus.energy()*sqrt(rt2/r2) <= highPtMin_) {
    if ((subdet==EcalBarrel) && (clus.energy()/tsosVtxMom.mag() > maxEOverPBarrel_)) return false;
    if ((subdet==EcalEndcap) && (clus.energy()/tsosVtxMom.mag() > maxEOverPEndcaps_)) return false;
    if ((subdet==EcalBarrel) && (clus.energy()/tsosVtxMom.mag() < minEOverPBarrel_)) return false;
    if ((subdet==EcalEndcap) && (clus.energy()/tsosVtxMom.mag() < minEOverPEndcaps_)) return false;
  }
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
  if (fabs(dphi)>CLHEP::pi) dphi = dphi < 0? CLHEP::pi2 + dphi : dphi - CLHEP::pi2;
  LogDebug("") << "delta phi : " << dphi;
  if (fabs(dphi) > maxDeltaPhi_) return false;
  LogDebug("") << "Delta phi criteria is satisfied ";

  //H/E cut
  if (HoE > maxHOverE_) return false; //FIXME: passe dans tous les cas?
  LogDebug("") << "H/E criteria is satisfied ";

  LogDebug("") << "electron has passed preselection criteria ";
  LogDebug("") << "=================================================";
  return true;  
}  

GlobalVector ElectronAlgoB::computeMode(const TrajectoryStateOnSurface &tsos) {

  // mode computation for momentum cartesian co-ordinates
  // change to 5D in local parameters??
  float mode_Px = 0.;
  float mode_Py = 0.;
  float mode_Pz = 0.;
  if ( tsos.isValid() ){
    std::vector<TrajectoryStateOnSurface> components(tsos.components());
    unsigned int numb = components.size();
    std::vector<SingleGaussianState1D> pxStates; pxStates.reserve(numb);
    std::vector<SingleGaussianState1D> pyStates; pyStates.reserve(numb);
    std::vector<SingleGaussianState1D> pzStates; pzStates.reserve(numb);
    for ( std::vector<TrajectoryStateOnSurface>::const_iterator ic=components.begin();
	  ic!=components.end(); ++ic ) {
      GlobalVector momentum(ic->globalMomentum());
      AlgebraicSymMatrix66 cov(ic->cartesianError().matrix());
      pxStates.push_back(SingleGaussianState1D(momentum.x(),cov(3,3),ic->weight()));
      pyStates.push_back(SingleGaussianState1D(momentum.y(),cov(4,4),ic->weight()));
      pzStates.push_back(SingleGaussianState1D(momentum.z(),cov(5,5),ic->weight()));
    }
    MultiGaussianState1D pxState(pxStates);
    MultiGaussianState1D pyState(pyStates);
    MultiGaussianState1D pzState(pzStates);
    GaussianSumUtilities1D pxUtils(pxState);
    GaussianSumUtilities1D pyUtils(pyState);
    GaussianSumUtilities1D pzUtils(pzState);
    mode_Px = pxUtils.mode().mean();
    mode_Py = pyUtils.mode().mean();
    mode_Pz = pzUtils.mode().mean();
  } else edm::LogInfo("") << "tsos not valid!!";
  return GlobalVector(mode_Px,mode_Py,mode_Pz);	

}

//FIXME!!
static const float R_ECAL           = 136.5;
static const float Z_Endcap         = 328.0;
static const float etaBarrelEndcap  = 1.479; 

float ElectronAlgoB::ecalEta(float EtaParticle , float Zvertex, float plane_Radius)
{
  if (EtaParticle!= 0.)
    {
      float Theta = 0.0  ;
      float ZEcal = (R_ECAL-plane_Radius)*sinh(EtaParticle)+Zvertex;
      
      if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
      if(Theta<0.0) Theta = Theta+Geom::pi() ;

      float ETA = - log(tan(0.5*Theta));
      
      if( fabs(ETA) > etaBarrelEndcap )
	{
	  float Zend = Z_Endcap ;
	  if(EtaParticle<0.0 )  Zend = -Zend ;
	  float Zlen = Zend - Zvertex ;
	  float RR = Zlen/sinh(EtaParticle);
	  Theta = atan((RR+plane_Radius)/Zend);
	  if(Theta<0.0) Theta = Theta+Geom::pi() ;
	  ETA = - log(tan(0.5*Theta));
	}
      return ETA;
    }
  else
    {
      edm::LogWarning("")  << "[EcalPositionFromTrack::etaTransformation] Warning: Eta equals to zero, not correcting" ;
      return EtaParticle;
    }
}

float ElectronAlgoB::ecalPhi(float PtParticle, float EtaParticle, float PhiParticle, int ChargeParticle, float Rstart)
{
  //Magnetic field
  const float RBARM = 1.357 ;  // was 1.31 , updated on 16122003
  const float ZENDM = 3.186 ;  // was 3.15 , updated on 16122003
  float Rbend = RBARM-(Rstart/100.); //Assumed Rstart in cm
  float Bend  = 0.3 * 4. * Rbend/ 2.0 ;

  //---PHI correction
  float PHI = 0.0 ;
  if( fabs(EtaParticle) <=  etaBarrelEndcap)
    {
      if (fabs(Bend/PtParticle)<=1.)
	{
	  PHI = PhiParticle - asin(Bend/PtParticle)*ChargeParticle;
	  if(PHI >  Geom::pi()) {PHI = PHI - Geom::twoPi();}
	  if(PHI < -Geom::pi()) {PHI = PHI + Geom::twoPi();}
	}
      else
	{
	  edm::LogWarning("") << "[EcalPositionFromTrack::phiTransformation] Warning:Too low Pt, giving up ";
	  return PhiParticle;
	}
    }
  
  if( fabs(EtaParticle) >  etaBarrelEndcap )
    {
      float Rhit = 0.0 ;
      Rhit = ZENDM / sinh(fabs(EtaParticle));
      if (fabs(((Rhit-(Rstart/100.))/Rbend)*Bend/PtParticle)<=1.)
	{
	  PHI = PhiParticle - asin(((Rhit-(Rstart/100.))/Rbend)*Bend/PtParticle)*ChargeParticle;
	  if(PHI >  Geom::pi()) {PHI = PHI - Geom::twoPi();}
	  if(PHI < -Geom::pi()) {PHI = PHI + Geom::twoPi();}
	}
      else
	{
	  edm::LogWarning("") <<"[EcalPositionFromTrack::phiTransformation] Warning:Too low Pt, giving up ";
	  return PhiParticle;
	}
      
    }
  
  //---Return the result
  return PHI;
}


const reco::GsfTrackRef
ElectronAlgoB::superClusterMatching(reco::SuperClusterRef sc, edm::Handle<reco::GsfTrackCollection> tracks) {

  double minDr = 0.5;
  double minDeop = 10.;
  //reco::SuperClusterRef theClus = edm::Ref<SuperClusterCollection>();
  reco::GsfTrackRef theTrack = edm::Ref<reco::GsfTrackCollection>();


  for(reco::GsfTrackCollection::size_type i=0; i<tracks->size(); ++i){
    reco::GsfTrackRef track(tracks, i);
    math::XYZVector trackGlobalDir(track->momentum());   
    math::XYZVector clusterGlobalDir(sc->x() - track->vx(), sc->y() - track->vy(), sc->z() - track->vz());
    //math::XYZVector clusterGlobalPos(sc->x(), sc->y(), sc->z());
    
    double tmpDr = ROOT::Math::VectorUtil::DeltaR(clusterGlobalDir, trackGlobalDir);
    if ( !(tmpDr < minDr) ) continue;

    TrajectoryStateOnSurface innTSOS = mtsTransform_->innerStateOnSurface(*track, *(trackerHandle_.product()), theMagField.product());
    GlobalVector innMom=computeMode(innTSOS);

    TrajectoryStateOnSurface outTSOS = mtsTransform_->outerStateOnSurface(*track, *(trackerHandle_.product()), theMagField.product());
    if (!outTSOS.isValid())   continue;

    TrajectoryStateOnSurface seedTSOS = TransverseImpactPointExtrapolator(*geomPropFw_).extrapolate(outTSOS,GlobalPoint(sc->seed()->position().x(),sc->seed()->position().y(),sc->seed()->position().z()));
    if (!seedTSOS.isValid()) seedTSOS=outTSOS;

    GlobalVector seedMom=computeMode(seedTSOS);

    double eOverPin  = sc->energy()/innMom.mag();
    double eOverPout = sc->seed()->energy()/seedMom.mag();
 
    double Deta = fabs(clusterGlobalDir.eta() - trackGlobalDir.eta());
    double Dphi = fabs(acos(cos(clusterGlobalDir.phi() - trackGlobalDir.phi())));

    //    if( !(eOverPout>0.5) ) continue;
    if( !(eOverPin<5) )  continue;
    if( !(Dphi < 0.2) )  continue;
    if( !(Deta < 0.02) ) continue;

    //    cout << " in matchbox, dphi, deta: " << Dphi << " , " << Deta << endl;
    //    cout << " in matchbox, E/Pin, out: " << eOverPin << " , " << eOverPout << endl;

    if( fabs(eOverPin-1.) < minDeop){
      minDeop = fabs(eOverPin-1.) ;
      theTrack = track;
    }
  }

  //cout << " in matchbox, minD(eop): " << minDeop << endl;
  //std::cout << "returning null ref" << std::endl;
  return theTrack;
}
