#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <TMVA/Reader.h>
#include <cmath>

inline double crop(const double val, const double min, const double max) {
    return std::min(max, std::max(min, val));
}

class DanieleBenedettiBDTProducer : public edm::EDProducer {
    public:
        explicit DanieleBenedettiBDTProducer(const edm::ParameterSet&);
        ~DanieleBenedettiBDTProducer();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag vtxLabel_;
        edm::InputTag eleLabel_;

        struct EleVariables {
          float fbrem;
          float detain;
          float dphiin;
          float sieie;
          float hoe;
          float eop;
          float e1x5e5x5;
          float eleopout;
          float detaeleout;
          float kfchi2;
          float kfhits;
          float mishits;
          float dist;
          float dcot;
          float nvtx;
         // spectators
          float eta;
          float pt;
          int ecalseed;
        } vars_;
        std::auto_ptr<TMVA:: Reader> reader_;
        std::string                  bdtName_;
};



DanieleBenedettiBDTProducer::DanieleBenedettiBDTProducer(const edm::ParameterSet& iConfig) :
    vtxLabel_(iConfig.getParameter<edm::InputTag>("vtxLabel")),
    eleLabel_(iConfig.getParameter<edm::InputTag>("eleLabel")),
    reader_(new TMVA::Reader("!Color:Silent")),
    bdtName_(iConfig.getParameter<std::string>("bdtName"))
{
    std::string xml = edm::FileInPath(iConfig.getParameter<std::string>("bdtXML")).fullPath();
    reader_->AddVariable("fbrem",       & vars_.fbrem);
    reader_->AddVariable("detain",      & vars_.detain);
    reader_->AddVariable("dphiin",      & vars_.dphiin);
    reader_->AddVariable("sieie",       & vars_.sieie);
    reader_->AddVariable("hoe",         & vars_.hoe);
    reader_->AddVariable("eop",         & vars_.eop);
    reader_->AddVariable("e1x5e5x5",    & vars_.e1x5e5x5);
    reader_->AddVariable("eleopout",    & vars_.eleopout);
    reader_->AddVariable("detaeleout",  & vars_.detaeleout);
    reader_->AddVariable("kfchi2",      & vars_.kfchi2);
    reader_->AddVariable("kfhits",      & vars_.kfhits);
    reader_->AddVariable("mishits",     & vars_.mishits);
    reader_->AddVariable("dist",        & vars_.dist);
    reader_->AddVariable("dcot",        & vars_.dcot);
    reader_->AddVariable("nvtx",        & vars_.nvtx); 
    reader_->AddSpectator("eta",        & vars_.eta); 
    reader_->AddSpectator("pt",         & vars_.pt); 
    reader_->AddSpectator("ecalseed",   & vars_.ecalseed); 
    reader_->BookMVA(bdtName_.c_str(), xml.c_str());

    produces<edm::ValueMap<float> >();
}

void DanieleBenedettiBDTProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace std;

    edm::Handle<edm::View<reco::GsfElectron> > eleH;
    iEvent.getByLabel(eleLabel_,eleH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxLabel_,vtxH);

    std::vector<float> bdtV(eleH->size(), 0.);
    std::auto_ptr<edm::ValueMap<float> > bdtM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler bdtF(*bdtM);

    vars_.nvtx = vtxH->size();
    for(size_t i = 0, n = eleH->size(); i < n; ++i) {
        const reco::GsfElectron &ele = eleH->at(i);

        //if (!ele.ecalDrivenSeed()) continue;

        if (ele.gsfTrack().isNull())  { std::cerr << "ERROR: missing gsfTrack (shouldn't happen!)" << std::endl; continue; }

        vars_.fbrem =  max<float>(ele.fbrem(), -1.);
        vars_.detain = min<float>(abs( ele.deltaEtaSuperClusterTrackAtVtx() ), 0.06);
        vars_.dphiin = min<float>(abs( ele.deltaPhiSuperClusterTrackAtVtx() ), 0.06);
        vars_.sieie = ele.sigmaEtaEta();
        vars_.hoe = ele.hadronicOverEm();
        vars_.eop = min<float>( ele.eSuperClusterOverP(), 20.);
        vars_.eleopout = min<float>( ele.eEleClusterOverPout(), 20. );
        vars_.e1x5e5x5 = crop( (ele.e5x5() - ele.e1x5()) / ele.e5x5() , -1., 2.);
        vars_.detaeleout = min<float>(abs( ele.deltaEtaEleClusterTrackAtCalo() ), 0.2);
        if (ele.closestCtfTrackRef().isNull()) { 
            vars_.kfhits = -1; // this is what is done
            vars_.kfchi2 = 0.; // in the training
        } else {
            vars_.kfhits = ele.closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement();
            vars_.kfchi2 = crop( ele.closestCtfTrackRef()->normalizedChi2(), 0., 15.);
        }
        vars_.mishits = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
        vars_.dist = min<float>(abs( ele.convDist() ), 15.);;
        vars_.dcot = min<float>(abs( ele.convDcot() ), 3.);

        vars_.pt = ele.pt(); 
        vars_.eta = ele.eta(); 
        vars_.ecalseed = ele.ecalDrivenSeed(); 
        bdtV[i] = reader_->EvaluateMVA(bdtName_.c_str());
    }

    bdtF.insert(eleH,bdtV.begin(),bdtV.end());

    bdtF.fill();
    iEvent.put(bdtM);

}

DanieleBenedettiBDTProducer::~DanieleBenedettiBDTProducer() { }
void DanieleBenedettiBDTProducer::beginJob() { }
void DanieleBenedettiBDTProducer::endJob() { }
DEFINE_FWK_MODULE(DanieleBenedettiBDTProducer);
