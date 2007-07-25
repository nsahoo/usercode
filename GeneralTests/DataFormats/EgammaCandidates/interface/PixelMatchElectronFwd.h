#ifndef EgammaReco_PixelMatchElectronFwd_h
#define EgammaReco_PixelMatchElectronFwd_h
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectron.h"
#include "DataFormats/EgammaCandidates/interface/GlobalCtfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GlobalGsfElectron.h"


namespace reco {
  class PixelMatchElectron;

  /// collection of PixelMatchElectron objects
  typedef std::vector<PixelMatchElectron> PixelMatchElectronCollection;

  /// reference to an object in a collection of PixelMatchElectron objects
  typedef edm::Ref<PixelMatchElectronCollection> PixelMatchElectronRef;

  /// reference to a collection of PixelMatchElectron objects
  typedef edm::RefProd<PixelMatchElectronCollection> PixelMatchElectronRefProd;

  /// vector of objects in the same collection of PixelMatchElectron objects
  typedef edm::RefVector<PixelMatchElectronCollection> PixelMatchElectronRefVector;

  /// iterator over a vector of reference to PixelMatchElectron objects
  typedef PixelMatchElectronRefVector::iterator pixelMatchElectron_iterator;


  //------------------------------------------------------
  /// collection of GlobalCtfElectron objects
  typedef std::vector<GlobalCtfElectron> GlobalCtfElectronCollection;

  /// reference to an object in a collection of GlobalCtfElectron objects
  typedef edm::Ref<GlobalCtfElectronCollection> GlobalCtfElectronRef;

  /// reference to a collection of GlobalCtfElectron objects
  typedef edm::RefProd<GlobalCtfElectronCollection> GlobalCtfElectronRefProd;

  /// vector of objects in the same collection of GlobalCtfElectron objects
  typedef edm::RefVector<GlobalCtfElectronCollection> GlobalCtfElectronRefVector;

  /// iterator over a vector of reference to GlobalCtfElectron objects
  typedef GlobalCtfElectronRefVector::iterator GlobalCtfElectron_iterator;


  //--------------- to be moved in a separete file -----------------
  /// collection of GlobalCtfElectron objects
  typedef std::vector<GlobalGsfElectron> GlobalGsfElectronCollection;

  /// reference to an object in a collection of GlobalGsfElectron objects
  typedef edm::Ref<GlobalGsfElectronCollection> GlobalGsfElectronRef;

  /// reference to a collection of GlobalGsfElectron objects
  typedef edm::RefProd<GlobalGsfElectronCollection> GlobalGsfElectronRefProd;

  /// vector of objects in the same collection of GlobalGsfElectron objects
  typedef edm::RefVector<GlobalGsfElectronCollection> GlobalGsfElectronRefVector;

  /// iterator over a vector of reference to GlobalGsfElectron objects
  typedef GlobalGsfElectronRefVector::iterator GlobalGsfElectron_iterator;


}

#endif
