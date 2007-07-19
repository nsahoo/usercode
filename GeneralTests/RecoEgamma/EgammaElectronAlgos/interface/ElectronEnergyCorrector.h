#ifndef ElectronEnergyCorrector_H
#define ElectronEnergyCorrector_H

//===================================================================
// Author: Federico Ferri - INFN Milano, Bicocca university 
//         Ivica Puljak - FESB, Split 
// 12/2005
//adapted to CMSSW by U.Berthon, dec 2006
//===================================================================

/*! \file ElectronEnergyCorrector.h
  Egamma class for Correction of electrons energy. 
*/
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GlobalCtfElectron.h"


class ElectronEnergyCorrector
{
 public:
  
  ElectronEnergyCorrector(){newEnergy_=0.;}

  float getCorrectedEnergy() const {return newEnergy_;}

  virtual void correct(reco::PixelMatchGsfElectron &);
  virtual void correct(reco::GlobalCtfElectron &);

 private:

  void setNewEnergy(const reco::PixelMatchGsfElectron &);
  void setNewEnergy(const reco::GlobalCtfElectron &);


  double fEtaBarrelBad(double scEta) const;
  double fEtaBarrelGood(double scEta) const;
  double fEtaEndcapBad(double scEta) const;
  double fEtaEndcapGood(double scEta) const;

  float newEnergy_;

};

#endif




