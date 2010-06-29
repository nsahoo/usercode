#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "CommonTools/RecoAlgos/interface/MuonSelector.h"

struct MuonSelector1 {
  typedef reco::MuonCollection collection;

  typedef std::vector<const reco::Muon *> container;

  typedef container::const_iterator const_iterator;

  MuonSelector1 ( const edm::ParameterSet & cfg ) :
    ptMin_( cfg.getParameter<double>( "ptMin" ) ), 
    eta_( cfg.getParameter<double>( "absEtaMax" ) ), 
    charge_( cfg.getParameter<int>( "charge" ) ) 
  { }

  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }


  void select( const edm::Handle<collection> & c, const edm::Event & event, const edm::EventSetup& eventSetup) {
    selected_.clear();
    for( reco::MuonCollection::const_iterator mu = c->begin(); 
         mu != c->end(); ++ mu )
      if ( mu->pt() > ptMin_        &&
	   abs(mu->eta()) < eta_    &&
	   mu->charge() == charge_     ) selected_.push_back( & * mu );
  }
  size_t size() const { return selected_.size(); }

private:
  container selected_;
  double ptMin_,eta_;
  int charge_;
};






// define your producer name
typedef ObjectSelector<MuonSelector1> MuonSelection1;

// declare the module as plugin
DEFINE_FWK_MODULE( MuonSelection1 );
