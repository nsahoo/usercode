// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/Common/interface/ValueMap.h>
#include <DataFormats/Common/interface/View.h>

#include<vector>

template<typename T>
class PatLeptonUserFloatAdder : public edm::EDProducer {
    public:
        explicit PatLeptonUserFloatAdder(const edm::ParameterSet&);
        ~PatLeptonUserFloatAdder();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag leptonTag_;
        std::vector<std::pair<std::string,edm::InputTag> > mapTags_;

};

template<typename T>
PatLeptonUserFloatAdder<T>::PatLeptonUserFloatAdder(const edm::ParameterSet& iConfig) :
    leptonTag_(iConfig.getUntrackedParameter<edm::InputTag>("leptonTag")) 
{
    edm::ParameterSet pset = iConfig.getParameter<edm::ParameterSet>("valueMaps");
    std::vector<std::string> names = pset.getParameterNamesForType<edm::InputTag>();
    mapTags_.reserve(names.size());
    for (unsigned int i = 0, n = names.size(); i < n; ++i) {
        mapTags_.push_back(std::make_pair(names[i], pset.getParameter<edm::InputTag>(names[i])));
    }
    produces<std::vector<T> >();  
}


template<typename T>
PatLeptonUserFloatAdder<T>::~PatLeptonUserFloatAdder() {}


// ------------ method called to produce the data  ------------
template<typename T>
void PatLeptonUserFloatAdder<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<T> > leptons;
    iEvent.getByLabel(leptonTag_,leptons);

    std::vector<edm::Handle<edm::ValueMap<float> > > maps_(mapTags_.size());
    for (unsigned int i = 0, n = mapTags_.size(); i < n; ++i) {
        iEvent.getByLabel(mapTags_[i].second, maps_[i]);
    }

    std::auto_ptr<std::vector<T> > out(new std::vector<T>());
    out->reserve(leptons->size());
    for (unsigned int j = 0, m = leptons->size(); j < m; ++j) {
        edm::Ptr<T> lref(leptons->ptrAt(j));
        out->push_back((*leptons)[j]); 
        T &lep = out->back();
        for (unsigned int i = 0, n = mapTags_.size(); i < n; ++i) {
            lep.addUserFloat( mapTags_[i].first, (*maps_[i])[lref] );
        }    
    }
    
    iEvent.put(out);
}

template<typename T>
void PatLeptonUserFloatAdder<T>::beginJob() { }
template<typename T>
void PatLeptonUserFloatAdder<T>::endJob() { } 

typedef PatLeptonUserFloatAdder<pat::Muon>     PatMuonUserFloatAdder;
typedef PatLeptonUserFloatAdder<pat::Electron> PatElectronUserFloatAdder;
DEFINE_FWK_MODULE(PatMuonUserFloatAdder);
DEFINE_FWK_MODULE(PatElectronUserFloatAdder);
