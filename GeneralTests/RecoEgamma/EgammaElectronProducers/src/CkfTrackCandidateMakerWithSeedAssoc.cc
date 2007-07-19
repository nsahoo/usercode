#include <memory>
#include <string>

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateSeedAssociation.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleanerBySharedHits.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoEgamma/EgammaElectronProducers/interface/CkfTrackCandidateMakerWithSeedAssoc.h"
#include "RecoTracker/CkfPattern/interface/TransientInitialStateEstimator.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"


using namespace edm;
using namespace std;

namespace cms{
  CkfTrackCandidateMakerWithSeedAssoc::CkfTrackCandidateMakerWithSeedAssoc(edm::ParameterSet const& conf) : 

    conf_(conf),theTrajectoryBuilder(0),theTrajectoryCleaner(0),
    theInitialState(0),theNavigationSchool(0)
  {  
    produces<TrackCandidateCollection>();  
    produces<reco::TrackCandidateSeedAssociationCollection>();
  }

  
  // Virtual destructor needed.
  CkfTrackCandidateMakerWithSeedAssoc::~CkfTrackCandidateMakerWithSeedAssoc() {
    delete theInitialState;  
    delete theNavigationSchool;
    delete theTrajectoryCleaner;    
  }  

  void CkfTrackCandidateMakerWithSeedAssoc::beginJob (EventSetup const & es)
  {
    //services
    es.get<TrackerRecoGeometryRecord>().get( theGeomSearchTracker );
    es.get<IdealMagneticFieldRecord>().get(theMagField);
      
    // get nested parameter set for the TransientInitialStateEstimator
    ParameterSet tise_params = conf_.getParameter<ParameterSet>("TransientInitialStateEstimatorParameters") ;
    theInitialState          = new TransientInitialStateEstimator( es,tise_params);
    theNavigationSchool      = new SimpleNavigationSchool(&(*theGeomSearchTracker),&(*theMagField));
    theTrajectoryCleaner     = new TrajectoryCleanerBySharedHits();          

    // set the correct navigation
    NavigationSetter setter( *theNavigationSchool);

    // set the TrajectoryBuilder
    std::string trajectoryBuilderName = conf_.getParameter<std::string>("TrajectoryBuilder");
    edm::ESHandle<TrackerTrajectoryBuilder> theTrajectoryBuilderHandle;
    es.get<CkfComponentsRecord>().get(trajectoryBuilderName,theTrajectoryBuilderHandle);
    theTrajectoryBuilder = theTrajectoryBuilderHandle.product();    
  }
  
  // Functions that gets called by framework every event
  void CkfTrackCandidateMakerWithSeedAssoc::produce(edm::Event& e, const edm::EventSetup& es)
  {        

    // Step A: set Event for the TrajectoryBuilder
    theTrajectoryBuilder->setEvent(e);        
    
    // Step B: Retrieve seeds
    
    std::string seedProducer = conf_.getParameter<std::string>("SeedProducer");
    std::string seedLabel = conf_.getParameter<std::string>("SeedLabel");
    edm::Handle<TrajectorySeedCollection> collseed;
    e.getByLabel(seedProducer, seedLabel, collseed);
    TrajectorySeedCollection theSeedColl = *collseed;
    
    // Step C: Create empty output collections
    std::auto_ptr<TrackCandidateCollection> output(new TrackCandidateCollection);    
    std::auto_ptr<reco::TrackCandidateSeedAssociationCollection> outAssoc(new reco::TrackCandidateSeedAssociationCollection);    
    // Step D: Invoke the building algorithm
      vector <int> seedLocations3;
    if ((*collseed).size()>0){
       vector<Trajectory> theFinalTrajectories;
       TrajectorySeedCollection::const_iterator iseed;
     
      vector<Trajectory> rawResult;
      vector <int> seedLocations;
      int seednr =0;
      for(iseed=theSeedColl.begin();iseed!=theSeedColl.end();iseed++){
	vector<Trajectory> theTmpTrajectories;
	theTmpTrajectories = theTrajectoryBuilder->trajectories(*iseed);
	
       
	LogDebug("CkfPattern") << "======== CkfTrajectoryBuilder returned " << theTmpTrajectories.size()
			       << " trajectories for this seed ========";

	theTrajectoryCleaner->clean(theTmpTrajectories);
      
	for(vector<Trajectory>::const_iterator it=theTmpTrajectories.begin();
	    it!=theTmpTrajectories.end(); it++){
	  if( it->isValid() ) {
	    rawResult.push_back(*it);
	    seedLocations.push_back(seednr);
	  }
	}
	LogDebug("CkfPattern") << "rawResult size after cleaning " << rawResult.size();
        seednr++;
      }
      
      // Step E: Clean the result
      seednr=0;
      vector<Trajectory> unsmoothedResult;
      vector <int> seedLocations2;
      theTrajectoryCleaner->clean(rawResult);
      
      for (vector<Trajectory>::const_iterator itraw = rawResult.begin();
	   itraw != rawResult.end(); itraw++) {
        if((*itraw).isValid()) {
	  unsmoothedResult.push_back( *itraw);
	  seedLocations2.push_back(seedLocations[seednr]);
	}
	  seednr++;
      }
      //analyseCleanedTrajectories(unsmoothedResult);
      

      // Step F: Convert to TrackCandidates
      seednr=0;
       for (vector<Trajectory>::const_iterator it = unsmoothedResult.begin();
	   it != unsmoothedResult.end(); it++) {
	
	OwnVector<TrackingRecHit> recHits;
	Trajectory::RecHitContainer thits = it->recHits();
	for (Trajectory::RecHitContainer::const_iterator hitIt = thits.begin();
	     hitIt != thits.end(); hitIt++) {
	  recHits.push_back( (**hitIt).hit()->clone());
	}
	
	//PTrajectoryStateOnDet state = *(it->seed().startingState().clone());
	std::pair<TrajectoryStateOnSurface, const GeomDet*> initState = 
	  theInitialState->innerState( *it);
      
	// temporary protection againt invalid initial states
	if (! initState.first.isValid() || initState.second == 0) {
          //cout << "invalid innerState, will not make TrackCandidate" << endl;
          continue;
        }

	PTrajectoryStateOnDet* state = TrajectoryStateTransform().persistentState( initState.first,
										   initState.second->geographicalId().rawId());
	//	FitTester fitTester(es);
	//	fitTester.fit( *it);
	
	output->push_back(TrackCandidate(recHits,it->seed(),*state));
	seedLocations3.push_back(seedLocations2[seednr]);
        delete state;
	seednr++;
       }
      
      
      edm::LogVerbatim("CkfPattern") << "========== CkfTrackCandidateMakerWithSeedAssoc Info ==========";
      edm::ESHandle<TrackerGeometry> tracker;
      es.get<TrackerDigiGeometryRecord>().get(tracker);
      edm::LogVerbatim("CkfPattern") << "number of Seed: " << theSeedColl.size();
      
      /*
      for(iseed=theSeedColl.begin();iseed!=theSeedColl.end();iseed++){
	DetId tmpId = DetId( iseed->startingState().detId());
	const GeomDet* tmpDet  = tracker->idToDet( tmpId );
	GlobalVector gv = tmpDet->surface().toGlobal( iseed->startingState().parameters().momentum() );
	
	edm::LogVerbatim("CkfPattern") << "seed perp,phi,eta : " 
				       << gv.perp() << " , " 
				       << gv.phi() << " , " 
				       << gv.eta() ;
      }
      */
      
      edm::LogVerbatim("CkfPattern") << "number of finalTrajectories: " << unsmoothedResult.size();
      for (vector<Trajectory>::const_iterator it = unsmoothedResult.begin();
	   it != unsmoothedResult.end(); it++) {
	edm::LogVerbatim("CkfPattern") << "n valid and invalid hit, chi2 : " 
	     << it->foundHits() << " , " << it->lostHits() <<" , " <<it->chiSquared();
      }
      edm::LogVerbatim("CkfPattern") << "=================================================";
          
    }

    // Step G: write output to file and create Associationmap
    const edm::OrphanHandle<TrackCandidateCollection> refprodTrackC = e.put(output);
 
    for (unsigned int i=0;i<seedLocations3.size();++i) {
      outAssoc->insert(edm::Ref<TrackCandidateCollection>(refprodTrackC,i),edm::Ref<TrajectorySeedCollection>(collseed,seedLocations3[i]));
     
    }
    
    e.put(outAssoc);
  }
}

