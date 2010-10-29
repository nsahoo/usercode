#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include <iostream>
#include <algorithm>
#include "Math/VectorUtil.h"
#include <DataFormats/PatCandidates/interface/JetCorrFactors.h>

// To be kept in synch with the enumerator definitions in SkimEvent.h file
const std::string reco::SkimEvent::hypoTypeNames[] = { "undefined", "WELNU", "WMUNU", "WWELEL", "WWELMU", "WWMUMU"};


std::string reco::SkimEvent::hypoTypeName(reco::SkimEvent::hypoType a){
  if(int(a) < int(hypoTypeSize) && int(a)>0) return hypoTypeNames[int(a)];
  return "undefined";
}



reco::SkimEvent::hypoType reco::SkimEvent::hypoTypeByName(const std::string &name){
  hypoType size = hypoTypeSize;
  int index = std::find(hypoTypeNames, hypoTypeNames+size, name)-hypoTypeNames;
  if(index == size) return undefined; // better this or throw() ?

  // cast
  return hypoType(index);
}



reco::SkimEvent::SkimEvent() : 
        hypo_(-1), vtxPoint_(0,0,0) { }

reco::SkimEvent::SkimEvent(const reco::SkimEvent::hypoType &h) :
        hypo_(h), vtxPoint_(0,0,0) { 
}


/*
void reco::SkimEvent::setExtraMuons(const edm::Handle<pat::MuonCollection> & muExtraH)				    
{                 
  extraMus_.clear();
  
  for(size_t i=0;i<muExtraH->size();++i){ 
    extraMus_.push_back(edm::RefToBase<reco::Candidate>( pat::MuonRef(muExtraH,i) ));
    for(size_t j=0; j<mus_.size(); ++j){
      if(extraMus_.back()->pt()==mus_[j]->pt() && extraMus_.back()->charge()==mus_[j]->charge()){
	extraMus_.pop_back();
	break;
      }
    }
  }
  sortExtraMusByPt();
}
*/

void reco::SkimEvent::setLepton(const pat::Electron& ele){
  leps_.push_back(ele);
}

void reco::SkimEvent::setLepton(const pat::Muon& mu){
  leps_.push_back(mu);
}

void reco::SkimEvent::setExtraLepton(const pat::Electron& ele){
  extraLeps_.push_back(ele);
}

void reco::SkimEvent::setExtraLepton(const pat::Muon& mu){
  extraLeps_.push_back(mu);
}

void reco::SkimEvent::setSoftMuon(const pat::Muon& mu){
  softMuons_.push_back(mu);
}


void reco::SkimEvent::setJets(const edm::Handle<pat::JetCollection> & jH) {

    jets_.clear();

    for(size_t i=0;i<jH->size();++i) 
        jets_.push_back(edm::RefToBase<reco::Candidate>(pat::JetRef(jH,i)));

    sortJetsByPt();

}

void reco::SkimEvent::setTCMet(const edm::Handle<reco::METCollection> & mH) {
    tcMet_ = reco::METRef(mH,0);
}

void reco::SkimEvent::setPFMet(const edm::Handle<reco::PFMETCollection> & mH) {
    pfMet_ = reco::PFMETRef(mH,0);
}


void reco::SkimEvent::setVertex(const edm::Handle<reco::VertexCollection> & vtxH) {
    
    if(vtxH->size() != 0) {
        vtxPoint_.SetXYZ( vtxH->at(0).x(),vtxH->at(0).y(),vtxH->at(0).z() );
	/*
        for(size_t i=1;i<vtxH->size();++i) {  //BM, better to chose the vertex independently of the leptons
            if( fabs(dZ(0,vtxH->at(i).position())) < fabs(dZ(0)) ) 
	      vtxPoint_.SetXYZ( vtxH->at(i).x(),vtxH->at(i).y(),vtxH->at(i).z() );
        }
	*/
    }

}


//Lepton variables

const int reco::SkimEvent::nLep(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = leps_.size(); 
    else for(size_t i=0;i<leps_.size();++i) if(leps_[i].pt() > minPt) count++;
    return count;
} 

const int reco::SkimEvent::nExtraLep(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = extraLeps_.size(); 
    else for(size_t i=0;i<extraLeps_.size();++i) if(extraLeps_[i].pt() > minPt) count++;
    return count;
} 

const int reco::SkimEvent::nSoftMu(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = softMuons_.size(); 
    else for(size_t i=0;i<softMuons_.size();++i) if(softMuons_[i].pt() > minPt) count++;
    return count;
} 

/*
const int reco::SkimEvent::nMu(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = mus_.size(); 
    else for(size_t i=0;i<mus_.size();++i) if(mus_[i].pt() > minPt) count++;
    return count;
} 

const int reco::SkimEvent::nExtraMu(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = extraMus_.size(); 
    else for(size_t i=0;i<extraMus_.size();++i) if(extraMus_[i].pt() > minPt) count++;
    return count;
} 

const int reco::SkimEvent::nEl(float minPt) const {
    int count = 0;
    if(minPt < 0) count = els_.size(); 
    else for(size_t i=0;i<els_.size();++i) if(els_[i].pt() > minPt) count++;
    return count;
} 

const int reco::SkimEvent::nExtraEl(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = extraEls_.size(); 
    else for(size_t i=0;i<extraEls_.size();++i) if(extraEls_[i].pt() > minPt) count++;
    return count;
} 
*/

const int reco::SkimEvent::pdgId(size_t i) const { 
    if(i < leps_.size()) return leps_[i].pdgId();
    else return -9999;
}


/*
const float reco::SkimEvent::sigmaIEtaIEta(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->sigmaIetaIeta() : 
              -9999.0);
}

const float reco::SkimEvent::deltaEtaIn(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->deltaEtaSuperClusterTrackAtVtx() : 
              -9999.0);
}

const float reco::SkimEvent::deltaPhiIn(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->deltaPhiSuperClusterTrackAtVtx() : 
              -9999.0);
}

const float reco::SkimEvent::hOverE(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->hadronicOverEm() : 
              -9999.0);
}

const int reco::SkimEvent::classification(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->classification() : 
              -9999.0);
}

const float reco::SkimEvent::isoTracker(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->dr03TkSumPt();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->isolationR03().sumPt;

    return -9999.0;
              
}

const float reco::SkimEvent::isoEcal(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->dr03EcalRecHitSumEt();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->isolationR03().emEt;

    return -9999.0;
}

const float reco::SkimEvent::isoHcal(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->dr03HcalTowerSumEt();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->isolationR03().hadEt;

    return -9999.0;
}

const float reco::SkimEvent::isoSum(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;
    else return isoTracker(i) + isoEcal(i) + isoHcal(i);
}
*/

const float reco::SkimEvent::pt(size_t i) const {

    if(i < leps_.size()) return leps_[i].pt();
    else return -9999.0;

}

const float reco::SkimEvent::eta(size_t i) const {

    if(i < leps_.size()) return leps_[i].eta();
    else return -9999.0;

}

const float reco::SkimEvent::phi(size_t i) const {

    if(i < leps_.size()) return leps_[i].phi();
    else return -9999.0;

}

/*
const int reco::SkimEvent::missingHits(size_t i) const {
  
    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->track()->trackerExpectedHitsInner().numberOfLostHits();

    return -9999;
}
*/

const float reco::SkimEvent::dcot(size_t i) const {

    return 0;
}

const float reco::SkimEvent::dist(size_t i) const {

    return 0;
}

/*
const bool reco::SkimEvent::isEE(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->isEE() : 
              -9999.0);
}

const bool reco::SkimEvent::isEB(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->isEB() : 
              -9999.0);
}

const float reco::SkimEvent::eSeedOverPout(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::ElectronRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::ElectronRef>()->eSeedClusterOverPout() : 
              -9999.0);
}

const bool reco::SkimEvent::isTracker(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::MuonRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::MuonRef>()->isTrackerMuon() : 
              -9999.0);
}

const bool reco::SkimEvent::isGlobal(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::MuonRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::MuonRef>()->isGlobalMuon() : 
              -9999.0);
}

const int reco::SkimEvent::nTrackerHits(size_t i) const {

    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->hitPattern().numberOfValidHits();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->track()->hitPattern().numberOfValidHits();

    return -9999;
}

const int reco::SkimEvent::nPixelHits(size_t i) const {

    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->hitPattern().numberOfValidPixelHits();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->track()->hitPattern().numberOfValidPixelHits();

    return -9999;
}

const int reco::SkimEvent::nChambers(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::MuonRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::MuonRef>()->numberOfChambers() : 
              -9999.0);
}

const int reco::SkimEvent::nMatches(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::MuonRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::MuonRef>()->numberOfMatches() : 
              -9999.0);
}

const int reco::SkimEvent::ndof(size_t i) const {

    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->ndof();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->combinedMuon()->ndof();

    return -9999;
}

const float reco::SkimEvent::chi2ndof(size_t i) const {

    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->normalizedChi2();
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->combinedMuon()->normalizedChi2();

    return -9999;
}

const int reco::SkimEvent::muonHits(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::MuonRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::MuonRef>()->combinedMuon()->hitPattern().numberOfValidHits() : 
              -9999.0);
}

const float reco::SkimEvent::dXY(size_t i) const {

    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return -1.0 * leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->dxy(vtxPoint_);
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return -1.0 * leps_[i].castTo<pat::MuonRef>()->track()->dxy(vtxPoint_);

    return -9999;
}

const float reco::SkimEvent::dZ(size_t i) const {

    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->dz(vtxPoint_);
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->track()->dz(vtxPoint_);

    return -9999;
}

const float reco::SkimEvent::dZ(size_t i, const math::XYZPoint &p) const {

    if( i >= leps_.size() ) return -9999;
    else if ( leps_[i].castTo<pat::ElectronRef>().isNonnull() )
        return leps_[i].castTo<pat::ElectronRef>()->gsfTrack()->dz(p);
    else if ( leps_[i].castTo<pat::MuonRef>().isNonnull() )
        return leps_[i].castTo<pat::MuonRef>()->track()->dz(p);

    return -9999;
}

const bool reco::SkimEvent::isTMOneStation(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::MuonRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::MuonRef>()->muonID("TMOneStation") :
              -9999.0);
}

const bool reco::SkimEvent::isTM2DLoose(size_t i) const {

    return ( (i < leps_.size() && leps_[i].castTo<pat::MuonRef>().isNonnull() ) ? 
              leps_[i].castTo<pat::MuonRef>()->muonID("TM2DLoose") :
              -9999.0);
}
*/

const int reco::SkimEvent::q(size_t i) const {
  if(i < leps_.size()) return leps_[i].charge();
  else return -9999.0;
}


//Jet variables
const int reco::SkimEvent::nJets(float minPt) const {
    int count = 0;
    if(minPt < 0) count = jets_.size(); 
    else for(size_t i=0;i<jets_.size();++i) if(jets_[i]->pt() > minPt) count++;
    return count;

    return 0;
}

const int reco::SkimEvent::nCentralJets(float minPt,float eta,bool applyCorrection) const {
    int count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(applyCorrection){
	if(!( ((pat::Jet*)(&*jets_[i]))->correctedP4(pat::JetCorrFactors::L3).Et() > minPt && fabs(jets_[i]->eta()) < eta) ) continue;
      }else{
	using namespace std;
	if(!  (((pat::Jet*)(&*jets_[i]))->correctedJet(pat::JetCorrFactors::Raw).pt() > minPt && fabs(jets_[i]->eta()) < eta) ) continue;
      }

      bool thisJetIsLepton(false);
      for(size_t j=0; j<leps_.size();++j){
	double dR = fabs(ROOT::Math::VectorUtil::DeltaR(jets_[i]->p4(),leps_[j].p4()) );
	//if(dR < 0.3){ //as in V1 reference selection
	if(dR < 0.4){ //as in V1 reference selection
	  thisJetIsLepton = true;
	  break;
	}
      }
      
      if(!thisJetIsLepton)  count++;
    }
    return count;
}

const float reco::SkimEvent::jetPt(size_t i) const {

    if(i < jets_.size()) return jets_[i]->pt();
    else return -9999.0;
}


//Event variables
const float reco::SkimEvent::pfMet() const {

    if(pfMet_.isNonnull()) return pfMet_->pt();
    else return -9999.0;
}

const float reco::SkimEvent::tcMet() const {

    if(tcMet_.isNonnull()) return tcMet_->pt();
    else return -9999.0;
}

const float reco::SkimEvent::mll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).mass();
}

const float reco::SkimEvent::pTll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).pt();
}

const float reco::SkimEvent::dPhill() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0].p4(),leps_[1].p4()) );   
}

const float reco::SkimEvent::dRll() const {
  if(leps_.size()!=2) return -9999.0;
  return ROOT::Math::VectorUtil::DeltaR(leps_[0].p4(),leps_[1].p4());
}

const float reco::SkimEvent::dEtall() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(leps_[0].eta() - leps_[1].eta());
}

const float reco::SkimEvent::etall() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).eta();
}

const float reco::SkimEvent::yll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).y();
}

const float reco::SkimEvent::dPhillPfMet() const {
  if(leps_.size()!=2 || pfMet_.isNonnull()) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0].p4()+leps_[1].p4(),pfMet_->p4()) );
}

const float reco::SkimEvent::dPhillTcMet() const {
  if(leps_.size()!=2 || tcMet_.isNonnull()) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0].p4()+leps_[1].p4(),tcMet_->p4()) );
}

/*
const float reco::SkimEvent::dPhilPfMet(size_t i) const {
  if(i < leps_.size() && pfMet_.isNonnull()) {
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[i].p4(),pfMet_->p4()) );
  } else 
    return -9999.0;
}

const float reco::SkimEvent::dPhilTcMet(size_t i) const {
  if(i < leps_.size() && tcMet_.isNonnull()) {
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[i].p4(),tcMet_->p4()) );
  } else 
    return -9999.0;
}
*/

/*
const size_t reco::SkimEvent::dPhiPfMetMin(size_t i, size_t j) const {
    if(i == 0 && j == 0) {
        size_t min = 0;
        for(size_t k=1;k<leps_.size();++k) {
            if( fabs(dPhilPfMet(k)) < fabs(dPhilPfMet(min)) ) min = k;
        }
        return min;
    } else {
        return (( dPhilPfMet(i) < dPhilPfMet(j) ) ? i : j );
    }
}

const size_t reco::SkimEvent::dPhiTcMetMin(size_t i, size_t j) const {
    if(i == 0 && j == 0) {
        size_t min = 0;
        for(size_t k=1;k<leps_.size();++k) {
            if( fabs(dPhilTcMet(k)) < fabs(dPhilTcMet(min)) ) min = k;
        }
        return min;
    } else {
        return (( dPhilTcMet(i) < dPhilTcMet(j) ) ? i : j );
    }
}
*/

const float reco::SkimEvent::projPfMet() const {
  double smallestDphi(999.);
  double projectedMET(0);
  for(size_t l=0; l<leps_.size();++l){
    double dphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(pfMet_->p4(),leps_[l].p4()) );
    if( dphi < smallestDphi){
      if(dphi < M_PI/2.)
	projectedMET = pfMet_->pt()*sin(dphi);
       else
	 projectedMET = pfMet_->pt();       
      smallestDphi = dphi;
    }
  }
  //std::cout << "projectedMET: " << projectedMET << std::endl;
  return projectedMET;

  /*from Matt
    if(i < leps_.size() && pfMet_.isNonnull()) {
        if(j==0) {
            return pfMet() * ( (dPhilPfMet(i) < M_PI/2) ? sin(dPhilPfMet(i)) : 1.0 );
        } else {
            size_t k = dPhiPfMetMin(i,j);
            return pfMet() * ( (dPhilPfMet(k) < M_PI/2) ? sin(dPhilPfMet(k)) : 1.0 );
        }
    }
    else {
        return -9999.0;
    }
  */
}



const float reco::SkimEvent::projTcMet() const {
  double smallestDphi(999.);
  double projectedMET(0);
  for(size_t l=0; l<leps_.size();++l){
    double dphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(tcMet_->p4(),leps_[l].p4()) );
    if( dphi < smallestDphi){
      if(dphi < M_PI/2.)
	projectedMET = tcMet_->pt()*sin(dphi);
       else
	 projectedMET = tcMet_->pt();       
      smallestDphi = dphi;
    }
  }
  //std::cout << "projectedMET: " << projectedMET << std::endl;
  return projectedMET;
  

  // implementation from Matt
  /*
  if(i < leps_.size() && tcMet_.isNonnull()) {
    if(j==0) {
      return tcMet() * ( (dPhilTcMet(i) < M_PI/2) ? sin(dPhilTcMet(i)) : 1.0 );
    } else {
      size_t k = dPhiTcMetMin(i,j);
      double projectedMET=tcMet() * ( (dPhilTcMet(k) < M_PI/2) ? sin(dPhilTcMet(k)) : 1.0 );
      std::cout << "projectedMET: " << projectedMET << std::endl;
      return projectedMET;
    }
  }
  else {
    return -9999.0;
  }
  */
}

/*
const float reco::SkimEvent::pfMT(size_t i) const {

    if(i < leps_.size() && pfMet_.isNonnull()) 
        return sqrt(2*pt(i)*pfMet()*(1 - cos(dPhilPfMet(i))));
    else
        return -9999.0;
}

const float reco::SkimEvent::tcMT(size_t i) const {

    if(i < leps_.size() && tcMet_.isNonnull()) 
        return sqrt(2*pt(i)*tcMet()*(1 - cos(dPhilTcMet(i))));
    else
        return -9999.0;
}
*/

/*
const float reco::SkimEvent::pfMTll(size_t i, size_t j) const {

    if(i < leps_.size() && j < leps_.size() && i!=j && pfMet_.isNonnull()) 
        return sqrt(2*pTll(i,j)*pfMet()*(1 - cos(dPhillPfMet(i,j))));
    else
        return -9999.0;
}

const float reco::SkimEvent::tcMTll(size_t i, size_t j) const {

    if(i < leps_.size() && j < leps_.size() && i!=j && tcMet_.isNonnull()) 
        return sqrt(2*pTll(i,j)*tcMet()*(1 - cos(dPhillTcMet(i,j))));
    else
        return -9999.0;
}

const float reco::SkimEvent::llType(size_t i, size_t j) const {

    return 0;
}
*/
const float reco::SkimEvent::nTracks() const {

    return 0;
}

/*
const float reco::SkimEvent::dPhiJetPfMet() const {

    return 0;
}

const float reco::SkimEvent::dPhiJetTcMet() const {

    return 0;
}
*/

const float reco::SkimEvent::cosThetaStar(size_t i) const {

    return 0;
}

const float reco::SkimEvent::afb(size_t i) const {

    return 0;
}


const bool reco::SkimEvent::leptEtaCut(float maxAbsEtaMu,float maxAbsEtaEl) const{
  bool check0(true);
  bool check1(true);

  if(fabs(leps_[0].pdgId())==11 && fabs(leps_[0].eta())>=maxAbsEtaEl) check0=false;
  if(fabs(leps_[0].pdgId())==13 && fabs(leps_[0].eta())>=maxAbsEtaMu) check0=false;

  if(fabs(leps_[1].pdgId())==11 && fabs(leps_[1].eta())>=maxAbsEtaEl) check1=false;
  if(fabs(leps_[1].pdgId())==13 && fabs(leps_[1].eta())>=maxAbsEtaMu) check1=false;

  return (check0 && check1);
}



const bool reco::SkimEvent::eleExpHitCut(bool isNew) const{
  using namespace std;
  for(size_t j=0; j<leps_.size(); ++j){
    if(fabs(leps_[j].pdgId())!=11) continue;
    int number;
    if(isNew){
      number = static_cast<const pat::Electron&>(leps_[j]).userInt("expectedHitsEle");
      //number = leps_[j].castTo<pat::ElectronRef>()->userInt("expectedHitsEle");
    }else{
      number = leps_[j].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
      //number = leps_[j].castTo<pat::ElectronRef>()->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
      //cout << "expected hits no meas: " << number << endl;
    }
    if(number!=0)  return false; 
  }
  return true;
}

