#include "NtupleMakers/ElectronAnalyzer/interface/NewElectrons.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenParticle.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronAlgoA.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

/// Constructor
NewElectrons::NewElectrons(const ParameterSet& pset) {
  fileName = pset.getParameter<std::string>("RootFileName");
  baselineEleCollName =  pset.getParameter<std::string>("BaselineEleCollName");
  customEleCollName   =  pset.getParameter<std::string>("CustomEleCollName");
}

NewElectrons::~NewElectrons() {}

void NewElectrons::beginJob(const EventSetup& eventSetup) {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event","Event data");

  tree->Branch("mc_n", &nMC, "mc_n/I");
  tree->Branch("mc_pt", mc_pt, "mc_pt[mc_n]/F");
  tree->Branch("mc_eta", mc_eta, "mc_eta[mc_n]/F");
  tree->Branch("mc_phi", mc_phi, "mc_phi[mc_n]/F");
  tree->Branch("mc_id", mc_id, "mc_id[mc_n]/I");
  tree->Branch("mc_e", mc_e, "mc_e[mc_n]/F");
  tree->Branch("mc_mother", mc_mother, "mc_mother[mc_n]/I");
  tree->Branch("mc_crack", mc_crack, "mc_crack[mc_n]/I"); // 0 no crack, 1 crack
 
  tree->Branch("sc_e", sc_e, "sc_e[mc_n]/F");
  tree->Branch("sc_et", sc_et, "sc_et[mc_n]/F");
  tree->Branch("sc_eta", sc_eta, "sc_eta[mc_n]/F");
  tree->Branch("sc_phi", sc_phi, "sc_phi[mc_n]/F");
  tree->Branch("sc_dr", sc_dr, "sc_dr[mc_n]/F");
  tree->Branch("sc_type", sc_type, "sc_type[mc_n]/I"); // 0 barrel, 1 endcap

  tree->Branch("tk_pt", tk_pt, "tk_e[mc_n]/F");
  tree->Branch("tk_eta", tk_eta, "tk_eta[mc_n]/F");
  tree->Branch("tk_phi", tk_phi, "tk_phi[mc_n]/F");
  tree->Branch("tk_dr", tk_dr, "tk_dr[mc_n]/F");
  tree->Branch("tk_nhit", tk_nhit, "tk_nhit[mc_n]/I");
  tree->Branch("tk_layer", tk_layer, "tk_layer[mc_n]/I");
  tree->Branch("tk_subdet", tk_subdet, "tk_subdet[mc_n]/I");

  tree->Branch("el_pt", el_pt, "el_pt[mc_n]/F");
  tree->Branch("el_e", el_e, "el_e[mc_n]/F");
  tree->Branch("el_eta", el_eta, "el_eta[mc_n]/F");
  tree->Branch("el_phi", el_phi, "el_phi[mc_n]/F");
  tree->Branch("el_dr", el_dr, "el_dr[mc_n]/F");
  tree->Branch("el_eopin", el_eopin, "el_eopin[mc_n]/F");
  tree->Branch("el_eopout", el_eopout, "el_eopout[mc_n]/F");
  tree->Branch("el_pout", el_pout, "el_pout[mc_n]/F");
  tree->Branch("el_fbrem", el_fbrem, "el_fbrem[mc_n]/F");
  tree->Branch("el_hoe", el_hoe, "el_hoe[mc_n]/F");
  tree->Branch("el_detain", el_detain, "el_detain[mc_n]/F");
  tree->Branch("el_dphiin", el_dphiin, "el_dphiin[mc_n]/F");
  tree->Branch("el_detaout", el_detaout, "el_detaout[mc_n]/F");
  tree->Branch("el_dphiout", el_dphiout, "el_dphiout[mc_n]/F");
  tree->Branch("el_e3x3", el_e3x3, "el_e3x3[mc_n]/F");
  tree->Branch("el_e5x5", el_e5x5, "el_e5x5[mc_n]/F");
  tree->Branch("el_eseed", el_eseed, "el_eseed[mc_n]/F");
  tree->Branch("el_spp", el_spp, "el_spp[mc_n]/F");
  tree->Branch("el_see", el_see, "el_see[mc_n]/F");
  tree->Branch("el_class", el_class, "el_class[mc_n]/I");

  tree->Branch("el1_pt", el1_pt, "el1_pt[mc_n]/F");
  tree->Branch("el1_e", el1_e, "el1_e[mc_n]/F");
  tree->Branch("el1_eta", el1_eta, "el1_eta[mc_n]/F");
  tree->Branch("el1_phi", el1_phi, "el1_phi[mc_n]/F");
  tree->Branch("el1_dr", el1_dr, "el1_dr[mc_n]/F");
  tree->Branch("el1_eopin", el1_eopin, "el1_eopin[mc_n]/F");
  tree->Branch("el1_eopout", el1_eopout, "el1_eopout[mc_n]/F");
  tree->Branch("el1_pout", el1_pout, "el1_pout[mc_n]/F");
  tree->Branch("el1_fbrem", el1_fbrem, "el1_fbrem[mc_n]/F");
  tree->Branch("el1_hoe", el1_hoe, "el1_hoe[mc_n]/F");
  tree->Branch("el1_detain", el1_detain, "el1_detain[mc_n]/F");
  tree->Branch("el1_dphiin", el1_dphiin, "el1_dphiin[mc_n]/F");
  tree->Branch("el1_detaout", el1_detaout, "el1_detaout[mc_n]/F");
  tree->Branch("el1_dphiout", el1_dphiout, "el1_dphiout[mc_n]/F");
  tree->Branch("el1_e3x3", el1_e3x3, "el1_e3x3[mc_n]/F");
  tree->Branch("el1_e5x5", el1_e5x5, "el1_e5x5[mc_n]/F");
  tree->Branch("el1_eseed", el1_eseed, "el1_eseed[mc_n]/F");
  tree->Branch("el1_spp", el1_spp, "el1_spp[mc_n]/F");
  tree->Branch("el1_see", el1_see, "el1_see[mc_n]/F");
}

void NewElectrons::endJob() {

  file->Write();
  file->Close();
}

void NewElectrons::analyze(const Event & event, const EventSetup& eventSetup) {

  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;

  Handle<HepMCProduct> evt;
  event.getByLabel("VtxSmeared", evt);
  myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));
    cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  Handle<PixelMatchGsfElectronCollection> elh;
  event.getByLabel(baselineEleCollName, elh);
  const PixelMatchGsfElectronCollection* copy = elh.product();
  PixelMatchGsfElectronCollection::const_iterator ite; 
    cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  Handle<SuperClusterCollection> sch1;
  event.getByLabel("hybridSuperClusters", sch1);
  const SuperClusterCollection* scb = sch1.product();
  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  Handle<SuperClusterCollection> sch2;
  event.getByLabel("islandSuperClusters", "islandEndcapSuperClusters", sch2);
  const SuperClusterCollection* sce = sch2.product();
  SuperClusterCollection::const_iterator itscb, itsce;
  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  Handle<TrackCollection> tkh;
  event.getByLabel("ctfWithMaterialTracks", tkh);
  const TrackCollection* tracks = tkh.product();
  TrackCollection::const_iterator itt;

  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  //Handle<GlobalCtfElectronCollection> elh1;
  Handle<PixelMatchGsfElectronCollection> elh1;
  event.getByLabel(customEleCollName, elh1);
  //const GlobalCtfElectronCollection*  electrons1 = elh1.product();
  //GlobalCtfElectronCollection::const_iterator ite1;
  const PixelMatchGsfElectronCollection*  electrons1 = elh1.product();
  PixelMatchGsfElectronCollection::const_iterator ite1;
  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  nMC = 0;
  
  for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) { 

    if ((abs((*it)->pdg_id()) == 11) && ((*it)->status() != 3)) {      
      if (((*it)->momentum().perp() > 5.) && (fabs((*it)->momentum().eta()) < 2.5)) {

        mc_mother[nMC] = mother(*it);
        mc_pt[nMC] = (*it)->momentum().perp();
        mc_eta[nMC] = (*it)->momentum().eta();
        mc_phi[nMC] = (*it)->momentum().phi();
        mc_e[nMC] = (*it)->momentum().e();
        mc_id[nMC] = (*it)->pdg_id();

        // check if it is in a crack
        if (inCrack(fabs((*it)->momentum().eta())))
          mc_crack[nMC] = 1;
        else
          mc_crack[nMC] = 0;

        math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
        int type = 0;
        float theta = 1.;
        double dR, dRmin = 0.2;
        SuperClusterCollection::const_iterator nearSC;
        
        // loop over barrel superclusters
        for(itscb = scb->begin(); itscb != scb->end(); ++itscb) {
          math::XYZVector scv(itscb->x(), itscb->y(), itscb->z());
          
          dR = ROOT::Math::VectorUtil::DeltaR(scv, mcv);
          
          if (dR < dRmin) {
            dRmin = dR;
            nearSC = itscb;
            type = 0;
            theta = scv.Theta();
          }
        }
       
        for(itsce = sce->begin(); itsce != sce->end(); ++itsce) {
          math::XYZVector scv(itsce->x(), itsce->y(), itsce->z());
          dR = ROOT::Math::VectorUtil::DeltaR(scv, mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearSC = itsce;
            type = 1;
            theta = scv.Theta();
          }
        }
        // store info about SC	  
        if (dRmin < 0.2) {
          sc_e[nMC] = nearSC->energy(); 
          sc_et[nMC] = sin(theta)*nearSC->energy();
          sc_eta[nMC] = nearSC->eta();
          sc_phi[nMC] = nearSC->phi(); 
          sc_dr[nMC] = dRmin;
          sc_type[nMC] = type;
        } else {
          sc_e[nMC] = 0.;
          sc_et[nMC] = 0.;
          sc_eta[nMC] = 0.;
          sc_phi[nMC] = 0.;
          sc_type[nMC] = -1;
          sc_dr[nMC] = 0.2; 
        }
  
        // remove duplicate electrons
        PixelMatchGsfElectronCollection electrons;
        PixelMatchGsfElectronCollection::const_iterator it1, it2;

        for(it1=copy->begin(); it1!=copy->end(); ++it1) {

          bool isRemoved = false;
          for(it2=copy->begin(); it2!=copy->end(); ++it2) {
            if (it1 == it2)
              continue;
            if (((*it1).superCluster().id() == (*it2).superCluster().id()) &&
                ((*it1).superCluster().index() == (*it2).superCluster().index())) {
              
              float deltaEp1 = fabs((*it1).eSuperClusterOverP() - 1.);
              float deltaEp2 = fabs((*it2).eSuperClusterOverP() - 1.);
              if (deltaEp1 > deltaEp2) {
                isRemoved = true;
                break;
              }
            }
          }
          
          if (!isRemoved)
            electrons.push_back(*it1);
        }

        dRmin = 0.1;
        PixelMatchGsfElectronCollection::const_iterator nearElectron;	
        for(ite = electrons.begin(); ite != electrons.end(); ++ite) {
          dR = ROOT::Math::VectorUtil::DeltaR(ite->p4(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearElectron = ite;
          }
        }

        // strore info about Ele
        if (dRmin < 0.1) {
          el_pt[nMC] = nearElectron->trackMomentumAtVtx().R(); 
          el_eta[nMC] = nearElectron->eta(); 
          el_e[nMC] = nearElectron->energy();
          el_phi[nMC] = nearElectron->phi(); 
          el_dr[nMC] = dRmin; 
          el_eopin[nMC] = nearElectron->eSuperClusterOverP();
          el_eopout[nMC] = nearElectron->eSeedClusterOverPout();
          el_hoe[nMC] = nearElectron->hadronicOverEm();
          el_dphiin[nMC] = nearElectron->deltaPhiSuperClusterTrackAtVtx();
          el_detain[nMC] = nearElectron->deltaEtaSuperClusterTrackAtVtx();
          el_dphiout[nMC] = nearElectron->deltaPhiSeedClusterTrackAtCalo();
          el_detaout[nMC] = nearElectron->deltaEtaSeedClusterTrackAtCalo();
          float pin  = nearElectron->trackMomentumAtVtx().R();
          float pout = nearElectron->trackMomentumOut().R();
          el_pout[nMC] = pout;
          el_fbrem[nMC] = (pin-pout)/pin;
          el_class[nMC] = nearElectron->classification();
          R9_25_gsf(event, &(*nearElectron), el_eseed[nMC], el_e3x3[nMC], el_e5x5[nMC], el_spp[nMC], el_see[nMC]);
        } else {
          el_pt[nMC] = 0.;
          el_e[nMC] = 0.;
          el_eta[nMC] = 0.;
          el_phi[nMC] = 0.;
          el_dr[nMC] = 0.1;
          el_eopin[nMC] = 0.;
          el_eopout[nMC] = 0.;
          el_pout[nMC] = 0;
          el_hoe[nMC] = 0.;
          el_dphiin[nMC] = 0.;
          el_detain[nMC] = 0.;
          el_dphiout[nMC] = 0.;
          el_detaout[nMC] = 0.;
          el_fbrem[nMC] = 0.;
          el_eseed[nMC] = 0.;
          el_e3x3[nMC] = 0.;
          el_e5x5[nMC] = 0.;
          el_spp[nMC] = 0.;
          el_see[nMC] = 0.;
          el_class[nMC] = -1;
        }
          cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
        // new electrons collection
        //GlobalCtfElectronCollection::const_iterator nearElectron1;
        PixelMatchGsfElectronCollection::const_iterator nearElectron1;
        dRmin = 0.1;
        for(ite1 = electrons1->begin(); ite1 != electrons1->end(); ++ite1) {
          dR = ROOT::Math::VectorUtil::DeltaR(ite1->p4(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearElectron1 = ite1;
          }
        }
          cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
        // strore info about Ele
        if (dRmin < 0.1) {
          el1_pt[nMC] = nearElectron1->trackMomentumAtVtx().R(); 
          el1_eta[nMC] = nearElectron1->eta(); 
          el1_e[nMC] = nearElectron1->energy();
          el1_phi[nMC] = nearElectron1->phi(); 
          el1_dr[nMC] = dRmin; 
          el1_eopin[nMC] = nearElectron1->eSuperClusterOverP();
          el1_eopout[nMC] = nearElectron1->eSeedClusterOverPout();
          el1_hoe[nMC] = nearElectron1->hadronicOverEm();
          el1_dphiin[nMC] = nearElectron1->deltaPhiSuperClusterTrackAtVtx();
          el1_detain[nMC] = nearElectron1->deltaEtaSuperClusterTrackAtVtx();
          el1_dphiout[nMC] = nearElectron1->deltaPhiSeedClusterTrackAtCalo();
          el1_detaout[nMC] = nearElectron1->deltaEtaSeedClusterTrackAtCalo();
          float pin  = nearElectron1->trackMomentumAtVtx().R();
          float pout = nearElectron1->trackMomentumOut().R();
          el1_pout[nMC] = pout;
          el1_fbrem[nMC] = (pin-pout)/pin;
          el1_class[nMC] = nearElectron1->classification();
          R9_25_gsf(event, &(*nearElectron1), el1_eseed[nMC], el1_e3x3[nMC], el1_e5x5[nMC], el1_spp[nMC], el1_see[nMC]);
        } else {
          el1_pt[nMC] = 0.;
          el1_e[nMC] = 0.;
          el1_eta[nMC] = 0.;
          el1_phi[nMC] = 0.;
          el1_dr[nMC] = 0.1; 
          el1_eopin[nMC] = 0.;
          el1_eopout[nMC] = 0.;
          el1_pout[nMC] = 0;
          el1_hoe[nMC] = 0.;
          el1_dphiin[nMC] = 0.;
          el1_detain[nMC] = 0.;
          el1_dphiout[nMC] = 0.;
          el1_detaout[nMC] = 0.;
          el1_fbrem[nMC] = 0.;
          el1_eseed[nMC] = 0.;
          el1_e3x3[nMC] = 0.;
          el1_e5x5[nMC] = 0.;
          el1_spp[nMC] = 0.;
          el1_see[nMC] = 0.;
          el1_class[nMC] = -1;
        }
        
        // loop over combinatorial track finder tracks
        dRmin = 0.1;
        TrackCollection::const_iterator nearTk;
        for(itt = tracks->begin(); itt != tracks->end(); ++itt) {
          dR = ROOT::Math::VectorUtil::DeltaR(itt->momentum(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearTk = itt;
          }
        }
         // strore info about Tk
        if (dRmin < 0.1) {
          tk_pt[nMC] = nearTk->pt();
          tk_nhit[nMC] = nearTk->found();
          
          // check subdetector
          TrackingRecHitRef hit = nearTk->recHit(0);
          if (hit->isValid()) {
            DetId detid(hit->geographicalId());       
            unsigned int subdetId = static_cast<unsigned int>(detid.subdetId());
            switch (subdetId) {
            case StripSubdetector::TIB:
              {
                TIBDetId theTIBDetId(detid.rawId());
                tk_layer[nMC] = theTIBDetId.layer();
                tk_subdet[nMC] = 3;
                break;
              }
            case StripSubdetector::TID:
              {
                TIDDetId theTIDDetId(detid.rawId());
                tk_layer[nMC] = theTIDDetId.wheel();
                tk_subdet[nMC] = 4;
                break;
              }
            case StripSubdetector::TOB:
              {
                TOBDetId theTOBDetId(detid.rawId());
                tk_layer[nMC] = theTOBDetId.layer();
                tk_subdet[nMC] = 5;
                break;
              }
            case StripSubdetector::TEC:
              {
                TECDetId theTECDetId(detid.rawId());             
                tk_layer[nMC] = theTECDetId.wheel();
                tk_subdet[nMC] = 6;
                break;
              }
            }
            
          }else {
            tk_layer[nMC] = -1;
            tk_subdet[nMC] = -1;
          }

          tk_eta[nMC] = nearTk->eta(); 
          tk_phi[nMC] = nearTk->phi(); 
          tk_dr[nMC] = dRmin;
        } else {
          tk_pt[nMC] = 0.;
          tk_nhit[nMC] = 0;
          tk_subdet[nMC] = -1;
          tk_layer[nMC] = -1;
          tk_eta[nMC] = 0.;
          tk_phi[nMC] = 0.;
          tk_dr[nMC] = 0.2; 
        }
        
        nMC++;
      } else {
        cout << (*it)->momentum().perp() << "  " << (*it)->momentum().eta() << endl;
      }
    }
  }

  if (nMC != 0)
    tree->Fill();
}

bool NewElectrons::inCrack(float eta) {

  return (eta < 0.018 ||
          (eta>0.423 && eta<0.461) ||
          (eta>0.770 && eta<0.806) ||
          (eta>1.127 && eta<1.163) ||
          (eta>1.460 && eta<1.558));
}

int NewElectrons::mother(HepMC::GenParticle *p) {
  
  HepMC::GenParticle *pMother = p->mother();
  if (pMother) {
    while(pMother->pdg_id() == p->pdg_id()) {
      pMother = pMother->mother();
    }
    return pMother->pdg_id();
  } else 
    return -1;
}

void NewElectrons::R9_25_gsf(const Event & event, const reco::PixelMatchGsfElectron* e,
                            float& eseed, float& e3x3, float& e5x5, float& spp, float& see) {
  
  reco::SuperClusterRef sclRef=e->superCluster();

  edm::Handle<reco::BasicClusterShapeAssociationCollection> bH, eH;
  event.getByLabel("hybridSuperClusters", "hybridShapeAssoc", bH);
  const reco::BasicClusterShapeAssociationCollection* barrelClShp = bH.product();
  event.getByLabel("islandBasicClusters", "islandEndcapShapeAssoc", eH);
  const reco::BasicClusterShapeAssociationCollection* endcapClShp = eH.product();

  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
  DetId id = sclRef->seed()->getHitsByDetId()[0];
  if (id.subdetId() == EcalBarrel) {
    seedShpItr = barrelClShp->find(sclRef->seed());
  } else {
    seedShpItr = endcapClShp->find(sclRef->seed());
  }

  // Get the ClusterShapeRef corresponding to the BasicCluster
  const reco::ClusterShapeRef& seedShapeRef = seedShpItr->val;

  eseed = sclRef->seed()->energy();
  e3x3 = seedShapeRef->e3x3();
  e5x5 = seedShapeRef->e5x5();
  spp = seedShapeRef->covPhiPhi();
  see = seedShapeRef->covEtaEta();
}

void NewElectrons::R9_25_ctf(const Event & event, const reco::GlobalCtfElectron* e,
                             float& eseed, float& e3x3, float& e5x5, float& spp, float& see) {

  reco::SuperClusterRef sclRef=e->superCluster();

  edm::Handle<reco::BasicClusterShapeAssociationCollection> bH, eH;
  event.getByLabel("hybridSuperClusters", "hybridShapeAssoc", bH);
  const reco::BasicClusterShapeAssociationCollection* barrelClShp = bH.product();
  event.getByLabel("islandBasicClusters", "islandEndcapShapeAssoc", eH);
  const reco::BasicClusterShapeAssociationCollection* endcapClShp = eH.product();

  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
  DetId id = sclRef->seed()->getHitsByDetId()[0];
  if (id.subdetId() == EcalBarrel) {
    seedShpItr = barrelClShp->find(sclRef->seed());
  } else {
    seedShpItr = endcapClShp->find(sclRef->seed());
  }

  // Get the ClusterShapeRef corresponding to the BasicCluster                                                                                
  const reco::ClusterShapeRef& seedShapeRef = seedShpItr->val;

  eseed = sclRef->seed()->energy();
  e3x3 = seedShapeRef->e3x3();
  e5x5 = seedShapeRef->e5x5();
  spp = sqrt(seedShapeRef->covPhiPhi());
  see = sqrt(seedShapeRef->covEtaEta());
}
