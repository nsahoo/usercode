// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/GaussianSumVertexFit/interface/GsfVertexTrackCompatibilityEstimator.h"
#include "WWAnalysis/AnalysisStep/interface/MySingleDeposit.h"

#include<vector>

#include "WWAnalysis/Tools/interface/VertexReProducer.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"

//for sim-reco IP bias
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


class PatMuonBooster : public edm::EDProducer {
    public:
        explicit PatMuonBooster(const edm::ParameterSet&);
        ~PatMuonBooster();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        template <class T> T findClosestVertex(const double zPos, 
                const std::vector<T>& vtxs);


        edm::InputTag muonTag_;
        edm::InputTag trackTag_;
        edm::InputTag vertexTag_;

};

PatMuonBooster::PatMuonBooster(const edm::ParameterSet& iConfig) :
        muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
        trackTag_(iConfig.getUntrackedParameter<edm::InputTag>("trackTag")),
        vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag"))
{
  produces<pat::MuonCollection>();  
}


PatMuonBooster::~PatMuonBooster() {

}


// ------------ method called to produce the data  ------------
void PatMuonBooster::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    ESHandle<TransientTrackBuilder> theTTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    edm::Handle<edm::View<reco::Candidate> > muons;
    iEvent.getByLabel(muonTag_,muons);

    edm::Handle<reco::TrackCollection > tkH;
    iEvent.getByLabel(trackTag_,tkH);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(vertexTag_,vertices);

    edm::Handle<reco::BeamSpot> bs;
    iEvent.getByLabel(edm::InputTag("offlineBeamSpot"),bs);


    reco::Vertex vertexYesB;
    reco::Vertex vertexNoB;

    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);

    // here I set the biased PV 
    if(vertices->empty()) 
        vertexYesB = reco::Vertex(reco::Vertex::Point(bs->position().x(),bs->position().y(),bs->position().z()),
                reco::Vertex::Error());

    VertexReProducer revertex(vertices, iEvent);
    Handle<reco::BeamSpot>        pvbeamspot; 
    iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);


    // ----- here is the real loop over the muons ----
    for(edm::View<reco::Candidate>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){    
        const pat::MuonRef musRef = edm::RefToBase<reco::Candidate>(muons,mu-muons->begin()).castTo<pat::MuonRef>();
        pat::Muon clone = *edm::RefToBase<reco::Candidate>(muons,mu-muons->begin()).castTo<pat::MuonRef>();
        if(clone.type() == 8) continue;
        reco::TransientTrack tt = theTTBuilder->build(musRef->reco::Muon::innerTrack());

        double zPos = tt.track().vz();
        //if(!vertices->empty()) vertexYesB = findClosestVertex<reco::Vertex>(zPos,*vertices);
	if(!vertices->empty()) vertexYesB = vertices->front(); //take the first in the list


        // -- add info wrt YesBias vertex
        Measurement1D ip = IPTools::absoluteTransverseImpactParameter(tt,vertexYesB).second;
        Measurement1D ip3D = IPTools::absoluteImpactParameter3D(tt,vertexYesB).second;
        clone.addUserFloat(std::string("tip"),ip.value());
        clone.addUserFloat(std::string("tipErr"),ip.error());
        clone.addUserFloat(std::string("ip"),ip3D.value());
        clone.addUserFloat(std::string("ipErr"),ip3D.error());


        // ------- here I add the information about the IP significance wrt the PV
        reco::TrackCollection newTkCollection;
        bool foundMatch(false);
        for(reco::Vertex::trackRef_iterator itk = vertexYesB.tracks_begin(); itk!= vertexYesB.tracks_end(); itk++){
            bool refMatching = (itk->get() == &*(musRef->reco::Muon::innerTrack()) );
            if(refMatching){
                foundMatch = true;
            }else{
                newTkCollection.push_back(*itk->get());
            }
        }//track collection for vertexNoB is set

	//cout << "checking mu matching" << endl;
        if(!foundMatch) {
	  //cout << "WARNING: no muon matching found" << endl;
	  vertexNoB = vertexYesB;
        }else{      
            vector<TransientVertex> pvs = revertex.makeVertices(newTkCollection, *pvbeamspot, iSetup) ;
            if(pvs.empty()) {
                vertexNoB = reco::Vertex(reco::Vertex::Point(bs->position().x(),bs->position().y(),bs->position().z()),
                        reco::Vertex::Error());
            } else {
	      //vertexNoB = findClosestVertex<TransientVertex>(zPos,pvs);
	      vertexNoB = pvs.front(); //take the first in the list
            }
        }


        Measurement1D ip_2 = IPTools::absoluteTransverseImpactParameter(tt,vertexNoB).second;
        Measurement1D ip3D_2 = IPTools::absoluteImpactParameter3D(tt,vertexNoB).second;

        clone.addUserFloat(std::string("tip2"),ip_2.value());
        clone.addUserFloat(std::string("tipErr2"),ip_2.error());
        clone.addUserFloat(std::string("ip2"),ip3D_2.value());
        clone.addUserFloat(std::string("ipErr2"),ip3D_2.error());


        // ------- OLD style information (for backward compatibility)
        clone.addUserFloat( "dxyPV",clone.track()->dxy(vertexYesB.position()) );
        clone.addUserFloat( "dzPV",clone.track()->dz(vertexYesB.position()) );
        clone.addUserFloat( "dxyPV2",clone.track()->dxy(vertexNoB.position()) );
        clone.addUserFloat( "dzPV2",clone.track()->dz(vertexNoB.position()) );

        const reco::CandidateBaseRef musRef2(muons,mu-muons->begin());

        //Beta Calculation
        float num0003 = 0, num0004 = 0, num0703 = 0, num0704 = 0;
        float den0003 = 0, den0004 = 0, den0703 = 0, den0704 = 0;
        for(size_t k=0;k<tkH->size();++k) {
            reco::TrackRef tk(tkH,k);
            if( ROOT::Math::VectorUtil::DeltaR(tk->momentum(),musRef->track()->momentum()) < 0.01 ) 
                continue;
            //All tracks 0.3
            if( ROOT::Math::VectorUtil::DeltaR(tk->momentum(),musRef->track()->momentum()) < 0.3 ) {
                den0003 += tk->pt();
                if( fabs( musRef->track()->dz(tk->vertex()) ) < 0.2 ) 
                    num0003 += tk->pt();
            }
            //All tracks 0.4
            if( ROOT::Math::VectorUtil::DeltaR(tk->momentum(),musRef->track()->momentum()) < 0.4 ) {
                den0004 += tk->pt();
                if( fabs( musRef->track()->dz(tk->vertex()) ) < 0.2 ) 
                    num0004 += tk->pt();
            }
            if(tk->pt() < 0.7) continue;
            //pt>0.7 tracks 0.3
            if( ROOT::Math::VectorUtil::DeltaR(tk->momentum(),musRef->track()->momentum()) < 0.3 ) {
                den0703 += tk->pt();
                if( fabs( musRef->track()->dz(tk->vertex()) ) < 0.2 ) 
                    num0703 += tk->pt();
            }
            //pt>0.7 tracks 0.4
            if( ROOT::Math::VectorUtil::DeltaR(tk->momentum(),musRef->track()->momentum()) < 0.4 ) {
                den0704 += tk->pt();
                if( fabs( musRef->track()->dz(tk->vertex()) ) < 0.2 ) 
                    num0704 += tk->pt();
            }
        }
        clone.addUserFloat("beta0003",((den0003==0)?1:num0003/den0003));
        clone.addUserFloat("beta0004",((den0004==0)?1:num0004/den0004));
        clone.addUserFloat("beta0703",((den0703==0)?1:num0703/den0703));
        clone.addUserFloat("beta0704",((den0704==0)?1:num0704/den0704));

//         clone.addUserFloat("sigma",*sigmaH);

        pOut->push_back(clone);
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatMuonBooster::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatMuonBooster::endJob() { } 

template <class T> T PatMuonBooster::findClosestVertex(const double zPos, 
        const std::vector<T>& vtxs){
    double dist(99999);
    T returnVertex;
    //unsigned int size = vtxs.size();
    typename std::vector<T>::const_iterator vtx;
    for(vtx = vtxs.begin(); vtx != vtxs.end(); vtx++){
        double tmpDist = fabs(zPos-vtx->position().z()); 
        if(tmpDist<dist){
            dist=tmpDist;
            returnVertex = *vtx;
        }
    }
    return returnVertex;
}


//define this as a plug-in
DEFINE_FWK_MODULE(PatMuonBooster);
