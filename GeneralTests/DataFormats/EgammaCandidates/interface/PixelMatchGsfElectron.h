#ifndef PixelMatchGsfElectron_h
#define PixelMatchGsfElectron_h
/** \class reco::Electron 
 *
 * An Electron with GsfTrack seeded from an ElectronPixelSeed
 * adapted from the TRecElectron class in ORCA
 *
 * \author U.Berthon, ClaudeCharlot,LLR
 *
 * \version $Id: PixelMatchGsfElectron.h,v 1.15 2007/03/16 12:46:30 uberthon Exp $
 *
 */

//-------------------------------------------------------------------
//
// Package EgammaCandidates
//
/** \class PixelMatchGsfElectron
*/
//  adapted from the TRecElectron class in ORCA
//
// Author:
//
// Claude Charlot - CNRS & IN2P3, LLR Ecole polytechnique
// Ursula Berthon - LLR Ecole polytechnique
// 
// $Log: PixelMatchGsfElectron.h,v $
// Revision 1.15  2007/03/16 12:46:30  uberthon
// make PixelMatchElectrons inherit from RecoCandidate
//
// Revision 1.11  2007/01/31 15:11:48  uberthon
// adapt to changes in GsfTrack class
//
// Revision 1.10  2007/01/26 13:06:02  uberthon
// add doxygen header
//
// Revision 1.9  2007/01/23 12:46:58  uberthon
// add outer track momentum
//
// Revision 1.8  2007/01/17 10:23:29  llista
// added virtual member function pdgId()
//
// Revision 1.7  2007/01/04 06:35:09  wmtan
// Geometry/Vector moved to DataFormats/GeometryVector
//
// Revision 1.6  2006/12/19 18:37:28  uberthon
// add electron classification and corrections
//
// Revision 1.5  2006/12/18 16:59:43  uberthon
// change to gsf mode calculations
//
// Revision 1.4  2006/12/04 17:47:18  uberthon
// make PixelMatchElectron +PixelMatchGsfElectron separate classes
//
// Revision 1.3  2006/11/14 18:52:22  uberthon
// add some missing data (HoE etc)
//
// Revision 1.2  2006/10/27 15:02:49  uberthon
// add PixelMatchGsfElectron
//
// Revision 1.1  2006/10/18 15:29:56  uberthon
// add PixelMatchGsfElectron class interface
//
// initial version
//
//
//-------------------------------------------------------------------

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"


#include <vector>

namespace reco {

class PixelMatchGsfElectron : public RecoCandidate {

 public:
  
  PixelMatchGsfElectron() {;} 

  PixelMatchGsfElectron(const SuperClusterRef scl, const GsfTrackRef gsft,
			const GlobalPoint tssuperPos, const GlobalVector tssuperMom, 
                        const GlobalPoint tsseedPos, const GlobalVector tsseedMom, 
                        const GlobalPoint innPos, const GlobalVector innMom, 
                        const GlobalPoint vtxPos, const GlobalVector vtxMom, 
                        const GlobalPoint outPos, const GlobalVector outMom, 
                        double HoE);

  virtual ~PixelMatchGsfElectron(){};

  //Public methods

  // particle behaviour
   /** The electron classification.
      barrel  :   0: golden,  10: bigbrem,  20: narrow, 30-34: showering, 
                (30: showering nbrem=0, 31: showering nbrem=1, 32: showering nbrem=2 ,33: showering nbrem=3, 34: showering nbrem>=4)
                 40: crack
      endcaps : 100: golden, 110: bigbrem, 120: narrow, 130-134: showering
               (130: showering nbrem=0, 131: showering nbrem=1, 132: showering nbrem=2 ,133: showering nbrem=3, 134: showering nbrem>=4)
   */
  int classification() const {return electronClass_;}

  PixelMatchGsfElectron * clone() const;

  // supercluster and electron track related quantities
  //! the super cluster energy corrected by EnergyScaleFactor
  float caloEnergy() const {return superClusterEnergy_;}
  //! the super cluster position
  math::XYZPoint caloPosition() const {return superCluster()->position();}
  //! the track momentum at vertex
  math::XYZVector trackMomentumAtVtx() const {return trackMomentumAtVtx_;}
  //! the track impact point state position
  math::XYZVector TrackPositionAtVtx() const {return trackPositionAtVtx_;}
  //! the track momentum extrapolated at the supercluster position
  math::XYZVector trackMomentumAtCalo() const {return trackMomentumAtCalo_;}
  //! the track momentum extrapolated from outermost position at the seed cluster position
   math::XYZVector trackMomentumOut() const {return trackMomentumOut_;}
  //! the track extrapolated position at min distance to the supercluster position
  math::XYZVector TrackPositionAtCalo() const {return trackPositionAtCalo_;}
  //! the supercluster energy / track momentum at impact point
  float eSuperClusterOverP() const {return eSuperClusterOverP_;}
  //! the seed cluster energy / track momentum at calo from outermost state
  float eSeedClusterOverPout() const {return eSeedClusterOverPout_;}
  //! the supercluster eta - track eta from helix extrapolation from impact point
  float deltaEtaSuperClusterTrackAtVtx() const {return deltaEtaSuperClusterAtVtx_;}
  //! the seed cluster eta - track eta at calo from outermost state
  float deltaEtaSeedClusterTrackAtCalo() const {return deltaEtaSeedClusterAtCalo_;}
  //! the supercluster phi - track phi from helix extrapolation from impact point
  float deltaPhiSuperClusterTrackAtVtx() const {return deltaPhiSuperClusterAtVtx_;}
  //! the seed cluster phi - track phi at calo from outermost state
  float deltaPhiSeedClusterTrackAtCalo() const {return deltaPhiSeedClusterAtCalo_;}

  //! the hadronic over electromagnetic fraction
  float hadronicOverEm() const {return hadOverEm_;}

  // corrections
  //! tell if class dependant escale correction have been applied
  bool isEnergyScaleCorrected() const {return energyScaleCorrected_;}
  //! tell if class dependant E-p combination has been determined
  bool isMomentumCorrected() const {return momentumFromEpCombination_;}
  //! handle electron energy correction.  Rescales 4 momentum from corrected 
  //! energy value and sets momentumFromEpCombination_ to true
  void correctElectronFourMomentum(const math::XYZTLorentzVectorD & momentum,float & enErr, float  & tMerr);
  //! handle electron supercluster energy scale correction.  Propagates new 
  //! energy value to all electron attributes and sets energyScaleCorrected_ to true
  void correctElectronEnergyScale(const float newEnergy);
  //! determine the class of the electron
  void classifyElectron(const int myclass);

  //! the errors on the supercluster energy and track momentum
  float caloEnergyError() const {return energyError_;}
  float trackMomentumError() const {return trackMomentumError_;}

  //! get associated superCluster Pointer
  SuperClusterRef superCluster() const { return superCluster_; } 

  //! get associated GsfTrack pointer
  reco::GsfTrackRef gsfTrack() const { return track_; } 

  //! number of related brem clusters
  int numberOfClusters() const {return superCluster_->clustersSize();}

  //! array of pointers to the related brem clusters
  //  BasicClusterRefVector getBremClusters() const;
  basicCluster_iterator basicClustersBegin() const { return superCluster_->clustersBegin(); }
  basicCluster_iterator basicClustersEnd() const { return superCluster_->clustersEnd(); }

  /// PDG identifier
  virtual int pdgId() const { return - 11 * charge(); }
  
 private:

  // temporary
  float ecalEta(float EtaParticle , float Zvertex, float plane_Radius);
  float ecalPhi(float PtParticle, float EtaParticle, float PhiParticle, int ChargeParticle, float Rstart);

  math::XYZVector trackMomentumAtVtx_;
  math::XYZVector trackPositionAtVtx_;
  math::XYZVector trackMomentumAtCalo_;
  math::XYZVector trackPositionAtCalo_;
  math::XYZVector trackMomentumOut_; 

  float energyError_;
  float trackMomentumError_;

  int electronClass_;

  float superClusterEnergy_;
  float eSuperClusterOverP_;
  float eSeedClusterOverPout_;

  float deltaEtaSuperClusterAtVtx_;
  float deltaEtaSeedClusterAtCalo_;
  float deltaPhiSuperClusterAtVtx_;
  float deltaPhiSeedClusterAtCalo_;

  float hadOverEm_;

  reco::SuperClusterRef superCluster_;
  reco::GsfTrackRef track_;

  bool energyScaleCorrected_;
  bool momentumFromEpCombination_;

  /// check overlap with another candidate
  virtual bool overlap( const Candidate & ) const;

};

}
#endif
