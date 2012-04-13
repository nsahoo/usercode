#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include <algorithm>
#include <cmath>

namespace utils2L {
    typedef StringObjectFunction<reco::Candidate,true>    Function;
    std::map<std::string,Function> functionCache2L_;
}

reco::SkimEvent2L::~SkimEvent2L() {}

reco::SkimEvent2L::SkimEvent2L(const reco::CompositeCandidate &src) :
    pat::CompositeCandidate(reco::CompositeCandidate(src.charge(), src.p4(), src.vertex(), src.pdgId(), src.status()))
{
    const reco::Candidate *srcl1 = src.daughter(0);
    const reco::Candidate *srcl2 = src.daughter(1);
    if (srcl1->pt() > srcl2->pt()) std::swap(srcl1, srcl2);
    addDaughter(*srcl1);
    addDaughter(*srcl2);
}

reco::SkimEvent2L::SkimEvent2L(const pat::CompositeCandidate &src) :
    pat::CompositeCandidate(src)
{
    clearDaughters();
    const reco::Candidate *srcl1 = src.daughter(0);
    const reco::Candidate *srcl2 = src.daughter(1);
    if (srcl1->pt() > srcl2->pt()) std::swap(srcl1, srcl2);
    addDaughter(*srcl1);
    addDaughter(*srcl2);
}

void reco::SkimEvent2L::setVertices(const edm::Handle<reco::VertexCollection> &vtxColl) 
{
    vertices_ = *vtxColl;
}

void reco::SkimEvent2L::setPFMet(const edm::Handle<reco::PFMETCollection> &pfMetColl)
{
    pfMet_ = pfMetColl->front();
}

void reco::SkimEvent2L::setPileupInfo(const edm::Handle<std::vector<PileupSummaryInfo> >& pu)
{
    for (size_t i = 0; i < pu->size(); i++) {
        if ((*pu)[i].getBunchCrossing() == 0) numsimvertices_ = (*pu)[i].getPU_NumInteractions();
    }
}

float reco::SkimEvent2L::lval(unsigned int il, const std::string &muExpr, const std::string &eleExpr) const 
{
    std::map<std::string,utils2L::Function>::const_iterator itMu = utils2L::functionCache2L_.find(muExpr);
    if (itMu == utils2L::functionCache2L_.end()) {
       itMu = utils2L::functionCache2L_.insert(std::make_pair(muExpr, utils2L::Function(muExpr))).first;
    }
    std::map<std::string,utils2L::Function>::const_iterator itEl = utils2L::functionCache2L_.find(eleExpr);
    if (itEl == utils2L::functionCache2L_.end()) {
       itEl = utils2L::functionCache2L_.insert(std::make_pair(eleExpr, utils2L::Function(eleExpr))).first;
    }
    const reco::Candidate &c = l(il);
    if      (typeid(c) == typeid(pat::Muon))     return itMu->second(c);
    else if (typeid(c) == typeid(pat::Electron)) return itEl->second(c);
    else throw cms::Exception("WrongType") << "Lepton " << il << " is of type " << typeid(c).name() << "\n";
}

float reco::SkimEvent2L::lval(unsigned int il, const std::string &expr) const
{
     std::map<std::string,utils2L::Function>::const_iterator itExpr = utils2L::functionCache2L_.find(expr);
     if (itExpr == utils2L::functionCache2L_.end()) {
        itExpr = utils2L::functionCache2L_.insert(std::make_pair(expr, utils2L::Function(expr))).first;
     }
     return itExpr->second(l(il));
}

