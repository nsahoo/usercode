#include "WWAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h"

#include <cmath>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include <TMatrixD.h>

void CompositeCandMassResolution::init(const edm::EventSetup &iSetup) {
    iSetup.get<IdealMagneticFieldRecord>().get(magfield_);
}

void CompositeCandMassResolution::getLeaves(const reco::Candidate &c, std::vector<const reco::Candidate *> &out) const {
    if (c.numberOfDaughters() > 0) {
        for (int i = 0, n = c.numberOfDaughters(); i < n; ++i) {
            getLeaves(*c.daughter(i), out);
        }
    } else if (c.hasMasterClonePtr()) {
        getLeaves(*c.masterClonePtr(), out);
    } else if (c.hasMasterClone()) {
        getLeaves(*c.masterClone(), out);
    } else {
        out.push_back(&c);
    }
}

double CompositeCandMassResolution::getMassResolution(const reco::Candidate &c) const  {
    std::vector<const reco::Candidate *> leaves;
    getLeaves(c, leaves);
    int n = leaves.size(), ndim = n*3;
    TMatrixDSym bigCov(ndim);
    TMatrixD jacobian(1,ndim);
    for (int i = 0, o = 0; i < n; ++i, o += 3) {
        const reco::Candidate &ci = *leaves[i];
        fillP3Covariance(ci, bigCov, o);
        jacobian(0, o+0) = (c.energy()*(ci.px()/ci.energy()) - c.px())/c.mass();
        jacobian(0, o+1) = (c.energy()*(ci.py()/ci.energy()) - c.py())/c.mass();
        jacobian(0, o+2) = (c.energy()*(ci.pz()/ci.energy()) - c.pz())/c.mass();
    }
    /*static int debug_ = 0;
    if (++debug_ < 20) {
        std::cout << "Big matrix:   " << std::endl; bigCov.Print();
        std::cout << "Jacobian:     " << std::endl; jacobian.Print();
    }*/
    TMatrixDSym massCov = bigCov.Similarity(jacobian);
    /*if (debug_ < 20) {
        std::cout << "Final matrix: " << std::endl; massCov.Print();
    }*/
    double dm2 = massCov(0,0);
    return (dm2 > 0 ? std::sqrt(dm2) : 0.0);
}

void CompositeCandMassResolution::fillP3Covariance(const reco::Candidate &c, TMatrixDSym &bigCov, int offset) const {
    const reco::GsfElectron *gsf; const reco::Muon *mu;
    if ((gsf = dynamic_cast<const reco::GsfElectron *>(&c)) != 0) {
        fillP3Covariance(*gsf, bigCov, offset);
    } else if ((mu = dynamic_cast<const reco::Muon *>(&c)) != 0) {
        fillP3Covariance(*mu, bigCov, offset);
    } else {
        throw cms::Exception("Unknown type", typeid(c).name());
    }
}

void CompositeCandMassResolution::fillP3Covariance(const reco::Muon &c, TMatrixDSym &bigCov, int offset) const {
    fillP3Covariance(c, *c.track(), bigCov, offset);
}


void CompositeCandMassResolution::fillP3Covariance(const reco::GsfElectron &c, TMatrixDSym &bigCov, int offset) const {
    if (c.ecalDriven()) {
          // In order to produce a 4x4 matrix, we need a jacobian from (p) to (px,py,pz), i.e.
          //            [ Px/P  ]                
          //  C_(3x3) = [ Py/P  ] * sigma^2(P) * [ Px/P Py/P Pz/P  ]
          //            [ Pz/P  ]                
          double dp = c.p4Error(reco::GsfElectron::P4_COMBINATION);
          AlgebraicMatrix31 ptop3;
          ptop3(0,0) = c.px()/c.p();
          ptop3(1,0) = c.py()/c.p();
          ptop3(2,0) = c.pz()/c.p();
          AlgebraicSymMatrix33 mat = ROOT::Math::Similarity(ptop3, AlgebraicSymMatrix11(dp*dp) );
          for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) {
              bigCov(offset+i,offset+j) = mat(i,j);
          } } 
    } else {
        fillP3Covariance(c, *c.gsfTrack(), bigCov, offset);
    }
}

void CompositeCandMassResolution::fillP3Covariance(const reco::Candidate &c, const reco::Track &t, TMatrixDSym &bigCov, int offset) const {
      GlobalTrajectoryParameters gp(GlobalPoint(t.vx(), t.vy(),  t.vz()),
                    GlobalVector(t.px(),t.py(),t.pz()),
                    t.charge(),
                    magfield_.product());
      JacobianCurvilinearToCartesian curv2cart(gp);
      CartesianTrajectoryError cartErr= ROOT::Math::Similarity(curv2cart.jacobian(), t.covariance());
      const AlgebraicSymMatrix66 mat = cartErr.matrix();
      for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) {
            bigCov(offset+i,offset+j) = mat(i+3,j+3);
      } } 
}

