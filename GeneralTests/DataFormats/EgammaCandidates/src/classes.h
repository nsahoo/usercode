//
// $Id: classes.h,v 1.15 2007/01/05 00:19:30 wmtan Exp $
//
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConvertedPhoton.h"
#include "DataFormats/EgammaCandidates/interface/SiStripElectron.h"
#include "DataFormats/EgammaCandidates/interface/PhotonIsolationAssociation.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/AssociationMap.h"


namespace {
  namespace {
    reco::PhotonCollection v1;
    edm::Wrapper<reco::PhotonCollection> w1;
    edm::Ref<reco::PhotonCollection> r1;
    edm::RefProd<reco::PhotonCollection> rp1;
    edm::RefVector<reco::PhotonCollection> rv1;

    reco::ElectronCollection v2;
    edm::Wrapper<reco::ElectronCollection> w2;
    edm::Ref<reco::ElectronCollection> r2;
    edm::RefProd<reco::ElectronCollection> rp2;
    edm::RefVector<reco::ElectronCollection> rv2;

    reco::PixelMatchElectronCollection v3;
    edm::Wrapper<reco::PixelMatchElectronCollection> w3;
    edm::Ref<reco::PixelMatchElectronCollection> r3;
    edm::RefProd<reco::PixelMatchElectronCollection> rp3;
    edm::RefVector<reco::PixelMatchElectronCollection> rv3;

    reco::PixelMatchGsfElectronCollection v4;
    edm::Wrapper<reco::PixelMatchGsfElectronCollection> w4;
    edm::Ref<reco::PixelMatchGsfElectronCollection> r4;
    edm::RefProd<reco::PixelMatchGsfElectronCollection> rp4;
    edm::RefVector<reco::PixelMatchGsfElectronCollection> rv4;

    reco::SiStripElectronCollection v5;
    edm::Wrapper<reco::SiStripElectronCollection> w5;
    edm::Ref<reco::SiStripElectronCollection> r5;
    edm::RefProd<reco::SiStripElectronCollection> rp5;
    edm::RefVector<reco::SiStripElectronCollection> rv5;

    reco::ConvertedPhotonCollection v6;
    edm::Wrapper<reco::ConvertedPhotonCollection> w6;
    edm::Ref<reco::ConvertedPhotonCollection> r6;
    edm::RefProd<reco::ConvertedPhotonCollection> rp6;
    edm::RefVector<reco::ConvertedPhotonCollection> rv6;

    reco::PhotonIsolationMap v66;
    edm::Wrapper<reco::PhotonIsolationMap> w66;
    edm::helpers::Key<edm::RefProd<reco::PhotonCollection > > h66;

    reco::ElectronIsolationMap v7;
    edm::Wrapper<reco::ElectronIsolationMap> w7;
    edm::helpers::Key<edm::RefProd<reco::ElectronCollection > > h7;

    reco::GlobalCtfElectronCollection v8;
    edm::Wrapper<reco::GlobalCtfElectronCollection> w8;
    edm::Ref<reco::GlobalCtfElectronCollection> r8;
    edm::RefProd<reco::GlobalCtfElectronCollection> rp8;
    edm::RefVector<reco::GlobalCtfElectronCollection> rv8;



    edm::reftobase::Holder<reco::Candidate, reco::ElectronRef> rb1;
    edm::reftobase::Holder<reco::Candidate, reco::PhotonRef> rb2;
    edm::reftobase::Holder<reco::Candidate, reco::SiStripElectronRef> rb3;
    edm::reftobase::Holder<reco::Candidate, reco::ConvertedPhotonRef> rb4;
  }
}
