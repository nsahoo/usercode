#include "WWAnalysis/AnalysisStep/interface/RegressionEnergyPatElectronProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <iostream>

using namespace edm ;
using namespace std ;
using namespace reco ;
using namespace pat ;

RegressionEnergyPatElectronProducer::RegressionEnergyPatElectronProducer( const edm::ParameterSet & cfg )
{

  produces<ElectronCollection>();

  inputPatElectrons = cfg.getParameter<edm::InputTag>("inputPatElectronsTag");
  energyRegressionType = cfg.getParameter<uint32_t>("energyRegressionType");
  regressionInputFile = cfg.getParameter<std::string>("regressionInputFile");
  debug = cfg.getParameter<bool>("debug");


  //****************************************************************************************
  //set up regression calculator
  //****************************************************************************************
  regressionEvaluator = new ElectronEnergyRegressionEvaluate();

  //set regression type
  ElectronEnergyRegressionEvaluate::ElectronEnergyRegressionType type = ElectronEnergyRegressionEvaluate::kNoTrkVar;
  if (energyRegressionType == 1) type = ElectronEnergyRegressionEvaluate::kNoTrkVar;
  else if (energyRegressionType == 2) type = ElectronEnergyRegressionEvaluate::kWithTrkVar;
  else if (energyRegressionType == 3) type = ElectronEnergyRegressionEvaluate::kNoTrkVarTwoPtBins;
  else if (energyRegressionType == 4) type = ElectronEnergyRegressionEvaluate::kWithTrkVarTwoPtBins;

  //load weights and initialize
  regressionEvaluator->initialize(edm::FileInPath(regressionInputFile.c_str()).fullPath(),type);

}
 
RegressionEnergyPatElectronProducer::~RegressionEnergyPatElectronProducer()
{
  delete regressionEvaluator;
}

void RegressionEnergyPatElectronProducer::produce( edm::Event & event, const edm::EventSetup & setup )
{

  assert(regressionEvaluator->isInitialized());


  //**************************************************************************
  //Get Number of Vertices
  //**************************************************************************
  Handle<reco::VertexCollection> hVertexProduct;
  event.getByLabel(edm::InputTag("goodPrimaryVertices"),hVertexProduct);
  const reco::VertexCollection inVertices = *(hVertexProduct.product());  

  // loop through all vertices
  Int_t nvertices = 0;
  for (reco::VertexCollection::const_iterator inV = inVertices.begin(); 
       inV != inVertices.end(); ++inV) {

    //pass these vertex cuts
    if (inV->ndof() >= 4
        && inV->position().Rho() <= 2.0
        && fabs(inV->z()) <= 24.0
      ) {
      nvertices++;
    }
  }

  //**************************************************************************
  //Get Rho
  //**************************************************************************
  double rho = 0;
  Handle<double> hRhoKt6PFJets;
  event.getByLabel(edm::InputTag("kt6PFJets","rho"), hRhoKt6PFJets);
  rho = (*hRhoKt6PFJets);

  edm::Handle<edm::View<reco::Candidate> > oldElectrons ;
  event.getByLabel(inputPatElectrons,oldElectrons) ;

  std::auto_ptr<ElectronCollection> electrons( new ElectronCollection ) ;

  ElectronCollection::const_iterator electron ;
  ElectronCollection::iterator ele ;
  // first clone the initial collection
  for(edm::View<reco::Candidate>::const_iterator ele=oldElectrons->begin(); ele!=oldElectrons->end(); ++ele) {    
    const pat::ElectronRef elecsRef = edm::RefToBase<reco::Candidate>(oldElectrons,ele-oldElectrons->begin()).castTo<pat::ElectronRef>();
    pat::Electron clone = *edm::RefToBase<reco::Candidate>(oldElectrons,ele-oldElectrons->begin()).castTo<pat::ElectronRef>();
    electrons->push_back(clone);
  }

  for
    ( ele = electrons->begin() ;
      ele != electrons->end() ;
      ++ele ) 
  {     


    bool printDebug = debug;
    if (printDebug) {
      std::cout << "***********************************************************************\n";
      std::cout << "Run Lumi Event: " << event.id().run() << " " << event.luminosityBlock() << " " << event.id().event() << "\n";
      std::cout << "Pat Electron : " << ele->pt() << " " << ele->eta() << " " << ele->phi() << "\n";
    }

    // apply regression energy
    Double_t RegressionMomentum = 0;

    //compute spp and sep
    double spp = 0;
    if (!isnan(ele->userFloat("covIPhiIPhi"))) spp = sqrt (ele->userFloat("covIPhiIPhi"));
    double sep;
    if (ele->sigmaIetaIeta()*spp > 0) {
      sep = ele->userFloat("covIEtaIPhi") / (ele->sigmaIetaIeta()*spp);
    } else if (ele->userFloat("covIEtaIPhi") > 0) {
      sep = 1.0;
    } else {
      sep = -1.0;
    }

    if (energyRegressionType == 1 || energyRegressionType == 3) {
      RegressionMomentum = regressionEvaluator->regressionValueNoTrkVar( ele->p(),
                                                                         ele->userFloat("rawEnergy"),
                                                                         ele->userFloat("eta"),
                                                                         ele->userFloat("phi"),
                                                                         ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                         ele->userFloat("etaWidth"),
                                                                         ele->userFloat("phiWidth"),
                                                                         ele->userInt("nBC"),
                                                                         ele->hadronicOverEm(),
                                                                         rho, 
                                                                         nvertices, 
                                                                         ele->userFloat("seedClusterEta"),
                                                                         ele->userFloat("seedClusterPhi"),
                                                                         ele->userFloat("seedClusterEnergy"),
                                                                         ele->userFloat("e3x3"),
                                                                         ele->userFloat("e5x5"),
                                                                         ele->sigmaIetaIeta(),
                                                                         spp,
                                                                         sep,
                                                                         ele->userFloat("eMax"),
                                                                         ele->userFloat("e2nd"),
                                                                         ele->userFloat("eTop"),
                                                                         ele->userFloat("eBottom"),
                                                                         ele->userFloat("eLeft"),
                                                                         ele->userFloat("eRight"),
                                                                         ele->userFloat("e2x5Max"),
                                                                         ele->userFloat("e2x5Top"),
                                                                         ele->userFloat("e2x5Bottom"),
                                                                         ele->userFloat("e2x5Left"),
                                                                         ele->userFloat("e2x5Right"),
                                                                         ele->pt(),
                                                                         ele->userFloat("seedieta"),
                                                                         ele->userFloat("seediphi"),
                                                                         ele->userFloat("etacryseed"),
                                                                         ele->userFloat("phicryseed"),
                                                                         ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                         printDebug);

    } else if (energyRegressionType == 2 || energyRegressionType == 4) {
      RegressionMomentum = regressionEvaluator->regressionValueWithTrkVar(ele->p(),
                                                                          ele->userFloat("rawEnergy"),
                                                                          ele->userFloat("eta"),
                                                                          ele->userFloat("phi"),
                                                                          ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                          ele->userFloat("etaWidth"),
                                                                          ele->userFloat("phiWidth"),
                                                                          ele->userInt("nBC"),
                                                                          ele->hadronicOverEm(),
                                                                          rho, 
                                                                          nvertices, 
                                                                          ele->userFloat("seedClusterEta"),
                                                                          ele->userFloat("seedClusterPhi"),
                                                                          ele->userFloat("seedClusterEnergy"),
                                                                          ele->userFloat("e3x3"),
                                                                          ele->userFloat("e5x5"),
                                                                          ele->sigmaIetaIeta(),
                                                                          spp,
                                                                          sep,
                                                                          ele->userFloat("eMax"),
                                                                          ele->userFloat("e2nd"),
                                                                          ele->userFloat("eTop"),
                                                                          ele->userFloat("eBottom"),
                                                                          ele->userFloat("eLeft"),
                                                                          ele->userFloat("eRight"),
                                                                          ele->userFloat("e2x5Max"),
                                                                          ele->userFloat("e2x5Top"),
                                                                          ele->userFloat("e2x5Bottom"),
                                                                          ele->userFloat("e2x5Left"),
                                                                          ele->userFloat("e2x5Right"),
                                                                          ele->pt(),
                                                                          ele->trackMomentumAtVtx().R(),
                                                                          ele->fbrem(),
                                                                          ele->charge(),
                                                                          ele->eSuperClusterOverP(),
                                                                          ele->userFloat("seedieta"),
                                                                          ele->userFloat("seediphi"),
                                                                          ele->userFloat("etacryseed"),
                                                                          ele->userFloat("phicryseed"),
                                                                          ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                          printDebug);
    } else {
      cout << "Error: RegressionType = " << energyRegressionType << " is not supported.\n";
    }
    
    math::XYZTLorentzVector oldMomentum = ele->p4();
    math::XYZTLorentzVector newMomentum = math::XYZTLorentzVector
      ( oldMomentum.x()*RegressionMomentum/oldMomentum.t(),
        oldMomentum.y()*RegressionMomentum/oldMomentum.t(),
        oldMomentum.z()*RegressionMomentum/oldMomentum.t(),
        RegressionMomentum ) ;

    ele->correctMomentum(newMomentum,ele->trackMomentumError(),ele->p4Error(reco::GsfElectron::P4_COMBINATION));
  }
  event.put(electrons) ;
     
}
  
  
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
DEFINE_FWK_MODULE(RegressionEnergyPatElectronProducer);
