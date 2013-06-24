#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>

#include<vector>
#include<TVector3.h>


class PatJetBooster : public edm::EDProducer {
    public:
        explicit PatJetBooster(const edm::ParameterSet&);
        ~PatJetBooster();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;


        edm::InputTag jetTag_;
        edm::InputTag vertexTag_;
        bool          storeJetId_;
        edm::InputTag jetIdTag_;
        edm::InputTag jetMvaTag_;
};

PatJetBooster::PatJetBooster(const edm::ParameterSet& iConfig) :
    jetTag_(iConfig.getParameter<edm::InputTag>("jetTag")),
    vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")) {

    storeJetId_ = iConfig.getUntrackedParameter<bool>("storeJetId",false);
    if ( storeJetId_ ) {
      jetIdTag_  = iConfig.getParameter<edm::InputTag>("jetIdTag") ;
      jetMvaTag_ = iConfig.getParameter<edm::InputTag>("jetMvaTag");
    }

    produces<pat::JetCollection>();  
}




void PatJetBooster::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::Candidate> > jetH;
    iEvent.getByLabel(jetTag_,jetH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vertexTag_,vtxH);
    reco::VertexCollection::const_iterator vtxLead = vtxH->begin();
    
    std::auto_ptr<pat::JetCollection> pOut(new pat::JetCollection);

    
    edm::Handle<edm::ValueMap<int> > JetIdH;
    edm::Handle<edm::ValueMap<float> > JetMvaH; 
    
    if ( storeJetId_ ) {
      iEvent.getByLabel(jetIdTag_,JetIdH);
      iEvent.getByLabel(jetMvaTag_,JetMvaH); 
    }

    reco::Vertex::Point pos(0,0,0);
    if(vtxH->size() > 0) pos = vtxH->at(0).position();

    float mom2, dz, sumDzMom2, sumMom2;
    for(edm::View<reco::Candidate>::const_iterator itJet=jetH->begin(); itJet!=jetH->end(); ++itJet) {    

        pat::Jet clone = *edm::RefToBase<reco::Candidate>(jetH,itJet-jetH->begin()).castTo<pat::JetRef>();
        pat::JetRef jetRef = edm::RefToBase<reco::Candidate>(jetH,itJet-jetH->begin()).castTo<pat::JetRef>();

        sumMom2=0; sumDzMom2=0;
        const reco::Track *thisTk;
        for(size_t i=0;i<clone.getPFConstituents().size();++i) {

            if(      clone.getPFConstituent(i)->gsfTrackRef().isNonnull() ) thisTk = &*clone.getPFConstituent(i)->gsfTrackRef();
            else if( clone.getPFConstituent(i)->trackRef().isNonnull()    ) thisTk = &*clone.getPFConstituent(i)->trackRef();
            else continue;

            mom2 = pow(clone.getPFConstituent(i)->pt(),2);
            dz   = thisTk->dz(pos);

            sumDzMom2 += mom2 * dz;
            sumMom2   += mom2;

        }

        if (sumMom2 != 0) {
            clone.addUserFloat("dz",sumDzMom2/sumMom2);
            clone.addUserFloat("mom2",sumMom2);
            clone.addUserFloat("dzMom2",sumDzMom2);
        } else {
            clone.addUserFloat("dz",-9999.9);
            clone.addUserFloat("mom2",-9999.9);
            clone.addUserFloat("dzMom2",-9999.9);
        }

        if ( storeJetId_ ) { 
           clone.addUserInt("jetId",(*JetIdH)[ jetRef ]);
           clone.addUserFloat("jetMva",(*JetMvaH)[ jetRef ]);
        } else {
          clone.addUserInt("jetId",-9);
          clone.addUserFloat("jetMva",-9999.9);
        }

        float ptd = -9999.9;
        ptd = clone.constituentPtDistribution();
        clone.addUserFloat("ptd",ptd);

        ///---- for QG discrimination ... and who knows what else
        ///---- see http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/tomc/QuarkGluonTagger/EightTeV/src/QGTagger.cc
//         Jet::Constituents constituents = clone->getJetConstituents();
        std::vector<reco::PFCandidatePtr> constituents = clone.getPFConstituents();
        float sum_pt2 = 0.;
        float sum_pt  = 0.;
        float sum_deta  = 0.;
        float sum_dphi  = 0.;
        float sum_deta2  = 0.;
        float sum_dphi2  = 0.;
        float sum_detadphi  = 0.;
        float sum_dR2  = 0.;
        float max_pt   = 0.;

        float QC_sum_pt2 = 0.;
        float QC_sum_pt  = 0.;
        float QC_sum_deta  = 0.;
        float QC_sum_dphi  = 0.;
        float QC_sum_deta2  = 0.;
        float QC_sum_dphi2  = 0.;
        float QC_sum_detadphi  = 0.;
        float QC_sum_dR2  = 0.;
        float QC_max_pt   = 0.;

        Int_t nChg_QC = 0, nChg_ptCut = 0, nNeutral_ptCut = 0;

        for( unsigned iConst=0; iConst<constituents.size(); ++iConst ) {
         if(!constituents[iConst].isNonnull()) continue;

         reco::TrackRef itrk = constituents[iConst]->trackRef();;

         bool trkForAxis = false;
         if(itrk.isNonnull()){                       //Track exists --> charged particle
          if(constituents[iConst]->pt() > 1.0) nChg_ptCut++;

          //Search for closest vertex to track
          reco::VertexCollection::const_iterator vtxClose = vtxH->begin();
          for(reco::VertexCollection::const_iterator vtx = vtxH->begin(); vtx != vtxH->end(); ++vtx){
           if(fabs(itrk->dz(vtx->position())) < fabs(itrk->dz(vtxClose->position()))) vtxClose = vtx;
          }

          if(vtxClose == vtxLead){
           Float_t dz = itrk->dz(vtxClose->position());
           Float_t dz_sigma = sqrt(pow(itrk->dzError(),2) + pow(vtxClose->zError(),2));

           if(itrk->quality(reco::TrackBase::qualityByName("highPurity")) && fabs(dz/dz_sigma) < 5.){
            trkForAxis = true;
            Float_t d0 = itrk->dxy(vtxClose->position());
            Float_t d0_sigma = sqrt(pow(itrk->d0Error(),2) + pow(vtxClose->xError(),2) + pow(vtxClose->yError(),2));
            if(fabs(d0/d0_sigma) < 5.) nChg_QC++;
           }
          }
         } else {                                //No track --> neutral constituents[iConst]icle
          if(constituents[iConst]->pt() > 1.0) nNeutral_ptCut++;
          trkForAxis = true;
         }


         float pt = constituents[iConst]->p4().Pt();
         float pt2 = pt*pt;

         if(trkForAxis){                   //If quality cuts, only use when trkForAxis
          QC_sum_pt += pt;
          QC_sum_pt2 += pt2;

          QC_sum_deta  += ((constituents[iConst]->eta() - clone.eta())*pt*pt);
          QC_sum_dphi  += ((2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*pt*pt);

          QC_sum_deta2 += ((constituents[iConst]->eta() - clone.eta())*(constituents[iConst]->eta() - clone.eta())*pt*pt);
          QC_sum_dphi2 += ((2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*(2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*pt*pt);

          QC_sum_detadphi += ((2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*(constituents[iConst]->eta() - clone.eta())*pt*pt);

          QC_sum_dR2    +=  ( (
            (constituents[iConst]->eta() - clone.eta())*(constituents[iConst]->eta() - clone.eta()) +
            (2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*(2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2))) 
                              )    *pt*pt);

          if (QC_max_pt < pt) QC_max_pt = pt;

         }

         //---- normal, without requirement of quality cuts
         sum_pt += pt;
         sum_pt2 += pt2;

         sum_deta  += ((constituents[iConst]->eta() - clone.eta())*pt*pt);
         sum_dphi  += ((2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*pt*pt);

         sum_deta2 += ((constituents[iConst]->eta() - clone.eta())*(constituents[iConst]->eta() - clone.eta())*pt*pt);
         sum_dphi2 += ((2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*(2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*pt*pt);

         sum_detadphi += ((2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*(constituents[iConst]->eta() - clone.eta())*pt*pt);

         sum_dR2    +=  ( (
           (constituents[iConst]->eta() - clone.eta())*(constituents[iConst]->eta() - clone.eta()) +
           (2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2)))*(2*atan(tan(((constituents[iConst]->phi()-clone.phi()))/2))) 
                          )    *pt*pt);

         if (max_pt < pt) max_pt = pt;

        }

        Float_t a = 0., b = 0., c = 0.;
        Float_t ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
        float axis1, axis2;

        //---- Quality cut
        if(QC_sum_pt2 > 0){
         clone.addUserFloat("QCptD",sqrt(QC_sum_pt2)/QC_sum_pt);
         clone.addUserFloat("QCRMScand",sqrt(QC_sum_dR2/QC_sum_pt2));
         clone.addUserFloat("QCRmax",QC_max_pt/QC_sum_pt);

         ave_deta = QC_sum_deta/QC_sum_pt2;
         ave_dphi = QC_sum_dphi/QC_sum_pt2;
         ave_deta2 = QC_sum_deta2/QC_sum_pt2;
         ave_dphi2 = QC_sum_dphi2/QC_sum_pt2;
         a = ave_deta2 - ave_deta*ave_deta;
         b = ave_dphi2 - ave_dphi*ave_dphi;
         c = -(QC_sum_detadphi/QC_sum_pt2 - ave_deta*ave_dphi);

         Float_t delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
         if(a+b-delta > 0) {
          axis1 = sqrt(0.5*(a+b+delta));
          axis2 = sqrt(0.5*(a+b-delta));
         }
         else {
          axis1 = 0.;
          axis2 = 0.;
         }
        }
        else {
         clone.addUserFloat("QCptD",0);
         clone.addUserFloat("QCRMScand",0);
         clone.addUserFloat("QCRmax",0);
         axis1 = 0.;
         axis2 = 0.;
        }

        clone.addUserFloat("QCaxis1",axis1);
        clone.addUserFloat("QCaxis2",axis2);

        //---- normal, without requirement of quality cuts
        if(sum_pt2 > 0){
         clone.addUserFloat("ptD",sqrt(sum_pt2)/sum_pt);
         clone.addUserFloat("RMScand",sqrt(sum_dR2/sum_pt2));
         clone.addUserFloat("Rmax",max_pt/sum_pt);

         ave_deta = sum_deta/sum_pt2;
         ave_dphi = sum_dphi/sum_pt2;
         ave_deta2 = sum_deta2/sum_pt2;
         ave_dphi2 = sum_dphi2/sum_pt2;
         a = ave_deta2 - ave_deta*ave_deta;
         b = ave_dphi2 - ave_dphi*ave_dphi;
         c = -(sum_detadphi/sum_pt2 - ave_deta*ave_dphi);

         Float_t delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
         if(a+b-delta > 0) {
          axis1 = sqrt(0.5*(a+b+delta));
          axis2 = sqrt(0.5*(a+b-delta));
         }
         else {
          axis1 = 0.;
          axis2 = 0.;
         }
        }
        else {
         clone.addUserFloat("ptD",0);
         clone.addUserFloat("RMScand",0);
         clone.addUserFloat("Rmax",0);
         axis1 = 0.;
         axis2 = 0.;
        }

        clone.addUserFloat("axis1",axis1);
        clone.addUserFloat("axis2",axis2);

        clone.addUserFloat("nChgQC",nChg_QC);
        clone.addUserFloat("nChgptCut",nChg_ptCut);
        clone.addUserFloat("nNeutralptCut",nNeutral_ptCut);

        pOut->push_back(clone);

    }
    iEvent.put(pOut);
}

PatJetBooster::~PatJetBooster() { }
void PatJetBooster::beginJob() { }
void PatJetBooster::endJob() { }
DEFINE_FWK_MODULE(PatJetBooster);
