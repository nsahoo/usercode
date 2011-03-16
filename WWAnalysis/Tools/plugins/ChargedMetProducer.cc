// system include files
#include <memory>
#include <cmath>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Utilities/interface/Exception.h"


class ChargedMetProducer : public edm::EDProducer {
    public:
        explicit ChargedMetProducer(const edm::ParameterSet&);
        ~ChargedMetProducer();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag collectionTag_;    
        edm::InputTag vertexTag_;
        double dzCut_;

};

ChargedMetProducer::ChargedMetProducer(const edm::ParameterSet& iConfig):
        collectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("collectionTag")),
        vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag")),
        dzCut_(iConfig.getUntrackedParameter<double>("dzCut")) {
        
    produces<edm::ValueMap<reco::PFMET> >();
}

void ChargedMetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::Vertex> > vtxH;
    iEvent.getByLabel(vertexTag_ ,vtxH);

    edm::Handle<edm::View<reco::Candidate> > leptonH;
    iEvent.getByLabel(collectionTag_,  leptonH);

    std::vector<math::XYZTLorentzVector> p4s(vtxH->size(),math::XYZTLorentzVector(0,0,0,0));
    std::vector<float> sumets(vtxH->size(),0);

    edm::View<reco::Candidate>::const_iterator itLep, endprobes = leptonH->end();
    float tempDz;
    for (itLep = leptonH->begin(); itLep != endprobes; ++itLep) {

        if( dynamic_cast<const reco::RecoCandidate*>(&*itLep) ) {
            tempDz =  dynamic_cast<const reco::RecoCandidate*>(&*itLep)->track().isNonnull() ?
                      dynamic_cast<const reco::RecoCandidate*>(&*itLep)->track()->vz() :
                      dynamic_cast<const reco::RecoCandidate*>(&*itLep)->gsfTrack()->vz() ;
        } else if( dynamic_cast<const reco::PFCandidate*>(&*itLep) ) {
            tempDz =  dynamic_cast<const reco::PFCandidate*>(&*itLep)->vz();
        } else throw cms::Exception("InvalidData") << "Data pointed to by collectionTag must inherit from reco::RecoCandidate or reco::PFCandidate\n";

        for(size_t j=0;j<vtxH->size();++j) {
            if(  fabs(tempDz - vtxH->at(j).z()) < dzCut_ ) {
                p4s[j] += itLep->p4();
                sumets[j] += itLep->pt();
            }
        }

    }

    std::vector<reco::PFMET> mets;
    CommonMETData output;
    for(size_t i=0;i<p4s.size();++i) {
        output.mex = -p4s[i].px();
        output.mey = -p4s[i].py();
        output.mez = -p4s[i].pz();
        output.met = p4s[i].pt();
        output.sumet = sumets[i];
        output.phi = atan2(-p4s[i].py(),-p4s[i].px());

        PFSpecificAlgo pf;
        mets.push_back(pf.addInfo(leptonH,output));
        
    }

    std::auto_ptr<edm::ValueMap<reco::PFMET> > metMap(new edm::ValueMap<reco::PFMET>);
    edm::ValueMap<reco::PFMET>::Filler filler(*metMap);
    filler.insert(vtxH, mets.begin(), mets.end());
    filler.fill();
    iEvent.put(metMap);

}

ChargedMetProducer::~ChargedMetProducer() { }
void ChargedMetProducer::beginJob() { }
void ChargedMetProducer::endJob() { } 
DEFINE_FWK_MODULE(ChargedMetProducer);
