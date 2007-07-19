#ifndef CkfTrajectoryBuilderWithSeedAssoc_H
#define CkfTrajectoryBuilderWithSeedAssoc_H

#include <vector>

class Propagator;
class TrajectoryStateUpdator;
class MeasurementEstimator;
class TrajectorySeed;
class TrajectoryStateOnSurface;
class TrajectoryFilter;

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include "RecoTracker/CkfPattern/interface/TrackerTrajectoryBuilder.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"



class TransientTrackingRecHitBuilder;


class CkfTrajectoryBuilderWithSeedAssoc :public TrackerTrajectoryBuilder {
protected:
// short names
  typedef TrajectoryStateOnSurface TSOS;
  typedef TrajectoryMeasurement TM;

public:

  typedef std::vector<Trajectory>     TrajectoryContainer;

  //CkfTrajectoryBuilder( const edm::ParameterSet& conf,
  //			const edm::EventSetup& es,
  //		const MeasurementTracker* theInputMeasurementTracker);
  CkfTrajectoryBuilderWithSeedAssoc(const edm::ParameterSet&              conf,
		       const TrajectoryStateUpdator*         updator,
		       const Propagator*                     propagatorAlong,
		       const Propagator*                     propagatorOpposite,
		       const Chi2MeasurementEstimatorBase*   estimator,
		       const TransientTrackingRecHitBuilder* RecHitBuilder,
		       const MeasurementTracker*             measurementTracker);

  ~CkfTrajectoryBuilderWithSeedAssoc();
  
  /// trajectories building starting from a seed
  virtual TrajectoryContainer trajectories(const TrajectorySeed& seed) const;

  /// set Event for the internal MeasurementTracker data member
  virtual void setEvent(const edm::Event& event) const;

private:
  const TrajectoryStateUpdator*         theUpdator;
  const Propagator*                     thePropagatorAlong;
  const Propagator*                     thePropagatorOpposite;
  const Chi2MeasurementEstimatorBase*   theEstimator;
  const TransientTrackingRecHitBuilder* theTTRHBuilder;
  const MeasurementTracker*             theMeasurementTracker;
  const LayerMeasurements*              theLayerMeasurements;

  // these may change from seed to seed
  mutable const Propagator*             theForwardPropagator;
  mutable const Propagator*             theBackwardPropagator;
  
  TrajectoryFilter*              theMinPtCondition;
  TrajectoryFilter*              theMaxHitsCondition;

  int theMaxCand;               /**< Maximum number of trajectory candidates 
		                     to propagate to the next layer. */
  int theMaxLostHit;            /**< Maximum number of lost hits per trajectory candidate.*/
  int theMaxConsecLostHit;      /**< Maximum number of consecutive lost hits 
                                     per trajectory candidate. */
  float theLostHitPenalty;      /**< Chi**2 Penalty for each lost hit. */
  bool theIntermediateCleaning;	/**< Tells whether an intermediary cleaning stage 
                                     should take place during TB. */
  int theMinimumNumberOfHits;   /**< Minimum number of hits for a trajectory to be returned.*/
  bool theAlwaysUseInvalidHits;


  Trajectory createStartingTrajectory( const TrajectorySeed& seed) const;

  std::vector<TrajectoryMeasurement> seedMeasurements(const TrajectorySeed& seed) const;

  void limitedCandidates( Trajectory& startingTraj, TrajectoryContainer& result) const;

  std::vector<TrajectoryMeasurement> findCompatibleMeasurements( const Trajectory& traj) const;

  bool qualityFilter( const Trajectory& traj) const;

  void addToResult( Trajectory& traj, TrajectoryContainer& result) const; 
  
  void updateTrajectory( Trajectory& traj, const TM& tm) const;

  bool toBeContinued( const Trajectory& traj) const;

};

#endif
