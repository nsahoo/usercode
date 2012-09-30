#include <vector>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"

#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "WWAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

//MELA stuff
//#include "WWAnalysis/AnalysisStep/interface/HZZ4lMelaDiscriminator.h"
#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "ZZMatrixElement/MELA/interface/PseudoMELA.h"
#include "ZZMatrixElement/MELA/interface/SpinOneEvenMELA.h"
#include "ZZMatrixElement/MELA/interface/SpinOneOddMELA.h"
#include "ZZMatrixElement/MELA/interface/SpinTwoMinimalMELA.h"


#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"
#include "TMVA/MethodBDT.h"

class SkimEvent4LProducer : public edm::EDProducer {
    public:
        SkimEvent4LProducer(const edm::ParameterSet &iConfig) ;
        ~SkimEvent4LProducer() ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag src_;

        edm::InputTag reducedPFCands_;
        edm::InputTag jets_;
        edm::InputTag pfMet_;
        edm::InputTag vertices_;
        edm::InputTag rho_;

        bool          isMC_;
        edm::InputTag gensTag_;
        std::string higgsmassweightfile_;
        HiggsMassWeightProvider* hmwp;
        bool          isSignal_;
        bool doswap;
        edm::InputTag mcMatch_;
        bool doMELA_;
        bool doAnglesWithFSR_;
        bool doMassRes_;
        bool doBDT_;

        std::auto_ptr<Mela> mela_;
        std::auto_ptr<PseudoMELA> pseudoMela_;
        std::auto_ptr<SpinOneEvenMELA> spinOneEvenMela_;
        std::auto_ptr<SpinOneOddMELA> spinOneOddMela_;
        std::auto_ptr<SpinTwoMinimalMELA> spinTwoMinimalMela_;
        CompositeCandMassResolution massRes_;

        //BDT input variables
        TMVA::Reader *ScalarVsBkgBDTReader;
        std::string weightfileScalarVsBkg_;
        float MVAInputVar_costheta1;
        float MVAInputVar_costheta2;
        float MVAInputVar_costhetastar;
        float MVAInputVar_Phi;
        float MVAInputVar_Phi1;
        float MVAInputVar_mZ1;
        float MVAInputVar_mZ2;
        float MVAInputVar_ZZpt;
        float MVAInputVar_zzdotz1;
        float MVAInputVar_zzdotz2;
        float MVAInputVar_ZZptZ1ptCosDphi;
        float MVAInputVar_ZZptZ2ptCosDphi;
        float MVAInputVar_reduced_Z1pt;
        float MVAInputVar_reduced_Z2pt;
        float MVAInputVar_ZZy;
        float MVAInputVar_m4l;
};

SkimEvent4LProducer::SkimEvent4LProducer(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    reducedPFCands_(iConfig.getParameter<edm::InputTag>("reducedPFCands")),
    jets_(iConfig.getParameter<edm::InputTag>("jets")),
    pfMet_(iConfig.getParameter<edm::InputTag>("pfMet")),
    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    isMC_(iConfig.getParameter<bool>("isMC")),
    gensTag_(isMC_ ? (iConfig.existsAs<bool>("gensTag") ? iConfig.getParameter<edm::InputTag>("gensTag") : edm::InputTag("prunedGen")) : edm::InputTag("NOGENHERE")),
    higgsmassweightfile_(iConfig.existsAs<std::string>("higgsMassWeightFile")?iConfig.getParameter<std::string>("higgsMassWeightFile"):""),
    isSignal_(iConfig.existsAs<bool>("isSignal")?iConfig.getParameter<bool>("isSignal"):false),
    doswap(iConfig.existsAs<bool>("doswap")?iConfig.getParameter<bool>("doswap"):true),
    mcMatch_(isSignal_ ? iConfig.getParameter<edm::InputTag>("mcMatch") : edm::InputTag("FAKE")),
    doMELA_(iConfig.existsAs<bool>("doMELA")?iConfig.getParameter<bool>("doMELA"):false),
    doAnglesWithFSR_(iConfig.existsAs<bool>("doAnglesWithFSR")?iConfig.getParameter<bool>("doAnglesWithFSR"):true),
    doMassRes_(iConfig.existsAs<bool>("doMassRes")?iConfig.getParameter<bool>("doMassRes"):false),
    doBDT_(iConfig.existsAs<bool>("doBDT")?iConfig.getParameter<bool>("doBDT"):false),
    weightfileScalarVsBkg_(iConfig.existsAs<std::string>("weightfile_ScalarVsBkgBDT")?iConfig.getParameter<std::string>("weightfile_ScalarVsBkgBDT"):"")
{
    if (doMELA_) {
      //mela_.reset(new Mela(true)); //use true to set the configuration as for ICHEP
      mela_.reset(new Mela(false));
      pseudoMela_.reset(new PseudoMELA);
      spinOneEvenMela_.reset(new SpinOneEvenMELA());
      spinOneOddMela_.reset(new SpinOneOddMELA());
      spinTwoMinimalMela_.reset(new SpinTwoMinimalMELA());
    }
    if (doBDT_) {
      ScalarVsBkgBDTReader = new TMVA::Reader( "V" );

      ScalarVsBkgBDTReader->AddVariable( "costheta1",		&MVAInputVar_costheta1);
      ScalarVsBkgBDTReader->AddVariable( "costheta2",		&MVAInputVar_costheta2);
      ScalarVsBkgBDTReader->AddVariable( "costhetastar",      &MVAInputVar_costhetastar);
      ScalarVsBkgBDTReader->AddVariable( "Phi",		&MVAInputVar_Phi);
      ScalarVsBkgBDTReader->AddVariable( "Phi1",		&MVAInputVar_Phi1);
      ScalarVsBkgBDTReader->AddVariable( "mZ1",		&MVAInputVar_mZ1);
      ScalarVsBkgBDTReader->AddVariable( "mZ2",		&MVAInputVar_mZ2);
      ScalarVsBkgBDTReader->AddVariable( "ZZpt/m4l",		&MVAInputVar_ZZpt);
      ScalarVsBkgBDTReader->AddVariable( "ZZdotZ1/(m4l*mZ1)",	&MVAInputVar_zzdotz1);
      ScalarVsBkgBDTReader->AddVariable( "ZZdotZ2/(m4l*mZ2)",	&MVAInputVar_zzdotz2);
      ScalarVsBkgBDTReader->AddVariable( "ZZptCosDphiZ1pt",	&MVAInputVar_ZZptZ1ptCosDphi);
      ScalarVsBkgBDTReader->AddVariable( "ZZptCosDphiZ2pt",	&MVAInputVar_ZZptZ2ptCosDphi);
      ScalarVsBkgBDTReader->AddVariable( "Z1pt/m4l",		&MVAInputVar_reduced_Z1pt);
      ScalarVsBkgBDTReader->AddVariable( "Z2pt/m4l",		&MVAInputVar_reduced_Z2pt);
      ScalarVsBkgBDTReader->AddVariable( "ZZy",		&MVAInputVar_ZZy);
      // add spectators
      ScalarVsBkgBDTReader->AddSpectator("m4l",               &MVAInputVar_m4l);
      // initialize
      ScalarVsBkgBDTReader->BookMVA("BDTG", edm::FileInPath(weightfileScalarVsBkg_.c_str()).fullPath());
      
    }

    const char *higgsweightbase=getenv("CMSSW_BASE");
    std::string higgsweightpath(higgsweightbase);
    higgsweightpath += "/src/";
    if (higgsmassweightfile_ != "") higgsweightpath += higgsmassweightfile_;
    else higgsweightpath = "";

    hmwp = new HiggsMassWeightProvider(higgsweightpath);

    produces<std::vector<reco::SkimEvent4L> >();
}

SkimEvent4LProducer::~SkimEvent4LProducer() {
    delete hmwp;
}

void
SkimEvent4LProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    if (doMassRes_) massRes_.init(iSetup);

    edm::Handle<reco::CandidateView> src; 
    iEvent.getByLabel(src_, src);

    edm::Handle<reco::VertexCollection> vertices; 
    iEvent.getByLabel(vertices_, vertices);
    edm::Handle<pat::JetCollection> jets; 
    iEvent.getByLabel(jets_, jets);
    edm::Handle<reco::PFMETCollection> pfMet; 
    iEvent.getByLabel(pfMet_, pfMet);
    edm::Handle<std::vector<reco::LeafCandidate> > pfleaves; 
    iEvent.getByLabel(reducedPFCands_, pfleaves);

    edm::Handle<edm::Association<reco::GenParticleCollection> > mcMatch; 
    if (isSignal_) iEvent.getByLabel(mcMatch_, mcMatch);

    edm::Handle<std::vector<PileupSummaryInfo> > puH; 
    if (isMC_) iEvent.getByLabel("addPileupInfo", puH);

    edm::Handle<reco::GenParticleCollection> gensH;
    if (isMC_) iEvent.getByLabel(gensTag_, gensH);

    float genhiggsmass = 0.0;
    if (isMC_) {
        reco::GenParticleCollection gens = *gensH;
        for (std::size_t i = 0; i < gens.size(); i++) {
            if (gens[i].pdgId() == 25 && gens[i].status() == 3) {
                genhiggsmass = gens[i].mass();
            }
        }
    }

    std::auto_ptr<std::vector<reco::SkimEvent4L> > out(new std::vector<reco::SkimEvent4L>());
    for (reco::CandidateView::const_iterator it = src->begin(), ed= src->end(); it != ed; ++it) {
        const reco::CompositeCandidate *srczz = dynamic_cast<const reco::CompositeCandidate *>(&*it);
        if (srczz == 0) throw cms::Exception("CorruptData") << "Found something that is not a composite candidate" << std::endl;
        const pat::CompositeCandidate *patzz = dynamic_cast<const pat::CompositeCandidate *>(srczz);
        if (patzz == 0) {
            const reco::Candidate *d0 = srczz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            const reco::Candidate *d00 = srczz->daughter(0)->daughter(0);
            if (d00 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)->daughter(0)\n";
            out->push_back( reco::SkimEvent4L(*srczz, doswap) );
        } else {
            const reco::Candidate *d0 = patzz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            const reco::Candidate *d00 = patzz->daughter(0)->daughter(0);
            if (d00 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)->daughter(0)\n";
            out->push_back( reco::SkimEvent4L(*patzz, doswap) );
        }
        reco::SkimEvent4L &zz = out->back();

        zz.setVertex(vertices);
        zz.setPFMet(pfMet);
        zz.setJets(jets);
        zz.setPFLeaves(pfleaves);
        zz.setNumRecoVertices(vertices);

        zz.setGenHiggsMass(genhiggsmass);
        zz.setHiggsMassWeight(hmwp);

        if (doMELA_) {
	  reco::Particle::LorentzVector lp4[2][2];
	  int lIds[2][2];
	  bool includeFSR(true);
	  for (unsigned int i = 0; i < 2; ++i) {
	    const reco::Candidate * l1 = & zz.l(i,0), * l2 = & zz.l(i,1);
	    if (l1->charge() > 0) std::swap(l1,l2);
	    lp4[i][0] = l1->p4();
	    lp4[i][1] = l2->p4();
	    lIds[i][0] = l1->pdgId();
	    lIds[i][1] = l2->pdgId();
	    if (includeFSR) {
	      for (unsigned int j = 2, n = zz.daughter(i)->numberOfDaughters(); j < n; ++j) {
                const reco::Candidate &pho = zz.l(i,j);
                int jnear = (reco::deltaR2(pho, *l1) <= reco::deltaR2(pho, *l2) ? 0 : 1);
                lp4[i][jnear] += pho.p4();
	      }
	    }
	  }
	  
	  TLorentzVector thep4M11(lp4[0][0].X(), lp4[0][0].Y(), lp4[0][0].Z(), lp4[0][0].R());
	  TLorentzVector thep4M12(lp4[0][1].X(), lp4[0][1].Y(), lp4[0][1].Z(), lp4[0][1].R());
	  TLorentzVector thep4M21(lp4[1][0].X(), lp4[1][0].Y(), lp4[1][0].Z(), lp4[1][0].R());
	  TLorentzVector thep4M22(lp4[1][1].X(), lp4[1][1].Y(), lp4[1][1].Z(), lp4[1][1].R());

	  float costhetastar,costheta1,costheta2,phi,phistar1,kd,psig,pbkg;
	  bool withPt(false),withY(false);
	  int LHCsqrt(8);
	  mela_->computeKD(thep4M11,lIds[0][0],
			   thep4M12,lIds[0][1],
			   thep4M21,lIds[1][0],
			   thep4M22,lIds[1][1],
			   costhetastar,costheta1,costheta2,phi,phistar1,
			   kd,psig,pbkg,
			   withPt,withY,LHCsqrt);
	  zz.setCosThetaStar(costhetastar);
	  zz.setCosTheta1(costheta1);
	  zz.setCosTheta2(costheta2);
	  zz.setPhi(phi);
	  zz.setPhiStar1(phistar1);

	  if(kd<0 || kd>1){
	    /* ------ DEBUG ------
	    cout << "mass: " << (thep4M11+thep4M12+thep4M21+thep4M22).M() << endl;
	    cout << "kd: " << kd << endl;
	    cout << "l11,id11,ld12,id12: "
		 << thep4M11.Px() << " , " 
		 << thep4M11.Py() << " , " 
		 << thep4M11.Pz() << " , " 
		 << thep4M11.E() << " , " 
		 << lIds[0][0] << " , " 
		 << thep4M12.Px() << " , " 
		 << thep4M12.Py() << " , " 
		 << thep4M12.Pz() << " , " 
		 << thep4M12.E() << " , " 
		 << lIds[0][1] << endl;
	    cout << "l21,id21,ld22,id22: "
		 << thep4M21.Px() << " , " 
		 << thep4M21.Py() << " , " 
		 << thep4M21.Pz() << " , " 
		 << thep4M21.E() << " , " 
		 << lIds[1][0] << " , "
		 << thep4M22.Px() << " , " 
		 << thep4M22.Py() << " , " 
		 << thep4M22.Pz() << " , " 
		 << thep4M22.E() << " , " 		 
		 << lIds[1][1] << endl;
	    */
	    kd = -1;

	  }

	  zz.addUserFloat("mela",kd); 	  
	  zz.addUserFloat("melaSig",psig); 	  
	  zz.addUserFloat("melaBkg",pbkg);
	  
	  // adding to other special mela values:
	  pseudoMela_->computeKD(thep4M11,lIds[0][0],
				 thep4M12,lIds[0][1],
				 thep4M21,lIds[1][0],
				 thep4M22,lIds[1][1],
				 kd,psig,pbkg);
	  zz.addUserFloat("melaPSLD",kd); 	  

	  spinOneEvenMela_->computeKD(thep4M11,lIds[0][0],
				      thep4M12,lIds[0][1],
				      thep4M21,lIds[1][0],
				      thep4M22,lIds[1][1],
				      kd,psig,pbkg);
	  zz.addUserFloat("melaSpinOneEven",kd); 	  

	  spinOneOddMela_->computeKD(thep4M11,lIds[0][0],
				     thep4M12,lIds[0][1],
				     thep4M21,lIds[1][0],
				     thep4M22,lIds[1][1],
				     kd,psig,pbkg);
	  zz.addUserFloat("melaSpinOneOdd",kd); 	  

	  spinTwoMinimalMela_->computeKD(thep4M11,lIds[0][0],
					 thep4M12,lIds[0][1],
					 thep4M21,lIds[1][0],
					 thep4M22,lIds[1][1],
					 kd,psig,pbkg);
	  zz.addUserFloat("melaSpinTwoMinimal",kd); 	  

        }

        if (doBDT_) {
          MVAInputVar_costheta1 = zz.getCosTheta1();
          MVAInputVar_costheta2 = zz.getCosTheta2();
          MVAInputVar_costhetastar =  zz.getCosThetaStar();
          MVAInputVar_Phi = zz.getPhi();
          //MVAInputVar_Phi1 = zz.getPhi1();
          MVAInputVar_Phi1 = zz.getPhiStar1();
          MVAInputVar_mZ1 = zz.mz(0);
          MVAInputVar_mZ2 = zz.mz(1);
          MVAInputVar_ZZpt = zz.pt();
          MVAInputVar_zzdotz1 = zz.p4().Dot(zz.z(0).p4());
          MVAInputVar_zzdotz2 = zz.p4().Dot(zz.z(1).p4());


          TLorentzVector Z1vec; Z1vec.SetPtEtaPhiM( zz.z(0).p4().Pt(), zz.z(0).p4().Eta(), zz.z(0).p4().Phi(), zz.z(0).p4().M() );
          TLorentzVector Z2vec; Z2vec.SetPtEtaPhiM( zz.z(1).p4().Pt(), zz.z(1).p4().Eta(), zz.z(1).p4().Phi(), zz.z(1).p4().M() );
          TLorentzVector ZZvec = Z1vec + Z2vec;
          TVector3 zzvec3 = Z1vec.Vect();
          TVector3 z1vec3 = Z2vec.Vect();
          TVector3 z2vec3 = ZZvec.Vect();
          TVector2 zzvecxy = zzvec3.XYvector();
          TVector2 z1vecxy = z1vec3.XYvector();
          TVector2 z2vecxy = z2vec3.XYvector();
          MVAInputVar_ZZptZ1ptCosDphi = TMath::Cos(zzvecxy.DeltaPhi( z1vecxy));
          MVAInputVar_ZZptZ2ptCosDphi = TMath::Cos(zzvecxy.DeltaPhi( z2vecxy));
          MVAInputVar_reduced_Z1pt = zz.z(0).p4().Pt() / zz.mass();
          MVAInputVar_reduced_Z2pt = zz.z(1).p4().Pt() / zz.mass();
          MVAInputVar_ZZy = zz.p4().Rapidity();
          MVAInputVar_m4l = zz.mass();
          zz.addUserFloat("BDT_ScalarVsBkg_125", ScalarVsBkgBDTReader->EvaluateMVA("BDTG"));
        }

        if (doMassRes_) zz.addUserFloat("massErr", massRes_.getMassResolution(zz));
        if (isMC_) zz.setPileupInfo(*puH);
        if (isSignal_) zz.setGenMatches(*mcMatch);
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent4LProducer);
