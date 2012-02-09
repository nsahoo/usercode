#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <iostream>
#include <algorithm>
using std::min;
using std::max;
using std::abs;

namespace utils {
    typedef StringObjectFunction<reco::Candidate,true>    Function;
    typedef StringCutObjectSelector<reco::Candidate,true> Selector;
    std::map<std::string,Function> functionCache_;
    std::map<std::string,Selector> selectorCache_;
} 

reco::SkimEvent4L::~SkimEvent4L() {}

reco::SkimEvent4L::SkimEvent4L(const reco::CompositeCandidate &src) :
    pat::CompositeCandidate(reco::CompositeCandidate(src.charge(), src.p4(), src.vertex(), src.pdgId(), src.status()))
{
    const reco::Candidate *srcz1 = src.daughter(0);
    const reco::Candidate *srcz2 = src.daughter(1);
    if (fabs(srcz1->mass()-91.188) > fabs(srcz2->mass()-91.188)) std::swap(srcz1, srcz2);
    addDaughter(*srcz1);
    addDaughter(*srcz2);
    init();
}

reco::SkimEvent4L::SkimEvent4L(const pat::CompositeCandidate &src) :
    pat::CompositeCandidate(src)
{
    if (numberOfDaughters() != src.numberOfDaughters()) throw cms::Exception("LogicError", "Didn't copy??");
    init();
}

void reco::SkimEvent4L::init() {
    if (abs(daughter(0)->pdgId()) == 11) {
        hypo_ = abs(daughter(1)->pdgId()) == 11 ? ZZ4EL : ZZ2EL2MU;
    } else {
        hypo_ = abs(daughter(1)->pdgId()) == 13 ? ZZ4MU : ZZ2MU2EL;
    }
}

void
reco::SkimEvent4L::setVertex(const edm::Handle<reco::VertexCollection> &vtxColl) 
{
    vtx_ = reco::VertexRef(vtxColl, 0);
    setVertex(vtx_->position());
}

void
reco::SkimEvent4L::setPFMet(const edm::Handle<reco::PFMETCollection> &pfMetColl)
{
    pfMet_ = reco::PFMETRef(pfMetColl,0);
}

void
reco::SkimEvent4L::setJets(const edm::Handle<pat::JetCollection> &jetColl, double ptMin)
{
    jets_.clear();
    for (unsigned int i = 0, n = jetColl->size(); i < n; ++i) {
        jets_.push_back(pat::JetRef(jetColl,i));
    }
}

void
reco::SkimEvent4L::setGenMatches(const edm::Association<reco::GenParticleCollection> &genMatch) 
{
    for (int i = 0; i < 2; ++i) {
        const reco::Candidate &recz = *daughter(i);
        std::vector<reco::GenParticleRef> & zGenMatches = (i == 0 ? z1GenMatches_ : z2GenMatches_);
        zGenMatches = std::vector<reco::GenParticleRef>();
        for (int j = 0, nj = recz.numberOfDaughters(); j < nj; ++j) {
            const reco::Candidate &proxy = *recz.daughter(j);
            reco::GenParticleRef match;
            if (proxy.hasMasterClonePtr())   match = genMatch[proxy.masterClonePtr()];
            else if (proxy.hasMasterClone()) match = genMatch[proxy.masterClone()];
            zGenMatches.push_back(match); 
        }
    }
}

const reco::Candidate & 
reco::SkimEvent4L::zByMass(unsigned int rank) const {
    bool z1first = (daughter(0)->mass() >= daughter(1)->mass());
    return z(z1first ? rank : 1-rank);
}

const reco::Candidate & 
reco::SkimEvent4L::lByPt(unsigned int rank) const {
    static std::pair<float,int> table_[4];
    table_[0].first = -daughter(0)->daughter(0)->pt(); table_[0].second = 0;
    table_[1].first = -daughter(0)->daughter(1)->pt(); table_[1].second = 1;
    table_[2].first = -daughter(1)->daughter(0)->pt(); table_[2].second = 2;
    table_[3].first = -daughter(1)->daughter(1)->pt(); table_[3].second = 3;
    std::sort(&table_[0], &table_[4]);
    return l(table_[rank].second / 2, table_[rank].second % 2);
}

float reco::SkimEvent4L::luserFloat(unsigned int iz, unsigned int il, const char *label)  const 
{
    return luserFloat(iz, il, std::string(label));
}

float reco::SkimEvent4L::luserFloat(unsigned int iz, unsigned int il, const std::string &label)  const 
{
    const reco::Candidate &c = l(iz,il);
    if      (typeid(c) == typeid(pat::Muon))     return (dynamic_cast<const pat::Muon &>(c)).userFloat(label);
    else if (typeid(c) == typeid(pat::Electron)) return (dynamic_cast<const pat::Electron &>(c)).userFloat(label);
    else throw cms::Exception("WrongType") << "Lepton " << iz << ", " << il << " is of type " << typeid(c).name() << "\n";
}


bool reco::SkimEvent4L::lgood(unsigned int iz, unsigned int il, const char *muId, const char *eleId) const 
{
    return lgood(iz,il, std::string(muId), std::string(eleId));
}

bool reco::SkimEvent4L::lgood(unsigned int iz, unsigned int il, const std::string &muId, const std::string &eleId) const 
{
     std::map<std::string,utils::Selector>::const_iterator itMuId = utils::selectorCache_.find(muId);
     if (itMuId == utils::selectorCache_.end()) {
        itMuId = utils::selectorCache_.insert(std::make_pair(muId, utils::Selector(muId))).first;
     }
     std::map<std::string,utils::Selector>::const_iterator itElId = utils::selectorCache_.find(eleId);
     if (itElId == utils::selectorCache_.end()) {
        itElId = utils::selectorCache_.insert(std::make_pair(eleId, utils::Selector(eleId))).first;
     }
     const reco::Candidate &c = l(iz,il);
     if      (typeid(c) == typeid(pat::Muon))     return itMuId->second(c);
     else if (typeid(c) == typeid(pat::Electron)) return itElId->second(c);
     else throw cms::Exception("WrongType") << "Lepton " << iz << ", " << il << " is of type " << typeid(c).name() << "\n";
}

float reco::SkimEvent4L::lval(unsigned int iz, unsigned int il, const char *expr) const 
{
    return lval(iz,il, std::string(expr));
}

float reco::SkimEvent4L::lval(unsigned int iz, unsigned int il, const std::string &expr) const 
{
     std::map<std::string,utils::Function>::const_iterator itExpr = utils::functionCache_.find(expr);
     if (itExpr == utils::functionCache_.end()) {
        itExpr = utils::functionCache_.insert(std::make_pair(expr, utils::Function(expr))).first;
     }
     return itExpr->second(l(iz,il));
}

float reco::SkimEvent4L::lval(unsigned int iz, unsigned int il, const char *muExpr, const char *eleExpr) const 
{
    return lval(iz,il, std::string(muExpr), std::string(eleExpr));
}

float reco::SkimEvent4L::lval(unsigned int iz, unsigned int il, const std::string &muExpr, const std::string &eleExpr) const 
{
     std::map<std::string,utils::Function>::const_iterator itMu = utils::functionCache_.find(muExpr);
     if (itMu == utils::functionCache_.end()) {
        itMu = utils::functionCache_.insert(std::make_pair(muExpr, utils::Function(muExpr))).first;
     }
     std::map<std::string,utils::Function>::const_iterator itEl = utils::functionCache_.find(eleExpr);
     if (itEl == utils::functionCache_.end()) {
        itEl = utils::functionCache_.insert(std::make_pair(eleExpr, utils::Function(eleExpr))).first;
     }
     const reco::Candidate &c = l(iz,il);
     if      (typeid(c) == typeid(pat::Muon))     return itMu->second(c);
     else if (typeid(c) == typeid(pat::Electron)) return itEl->second(c);
     else throw cms::Exception("WrongType") << "Lepton " << iz << ", " << il << " is of type " << typeid(c).name() << "\n";
}

reco::GenParticleRef reco::SkimEvent4L::genz(int iz) const 
{
    reco::GenParticleRef gl1 = genl(iz,0), gl2 = genl(iz,1);
    // search for mothers = not leptons
    while (gl1.isNonnull() && gl1->numberOfMothers() != 0 && 11 <= abs(gl1->pdgId()) && abs(gl1->pdgId()) <= 15) {
        gl1 = gl1->motherRef(0);
    }
    while (gl2.isNonnull() && gl2->numberOfMothers() != 0 && 11 <= abs(gl2->pdgId()) && abs(gl2->pdgId()) <= 15) {
        gl2 = gl2->motherRef(0);
    }
    if (gl1.isNonnull() && gl2.isNonnull() && gl1 == gl2) {
        return gl1;
    } else {
        return reco::GenParticleRef();
    }
}


reco::GenParticleRef reco::SkimEvent4L::genl(int iz, int il) const 
{
    const std::vector<reco::GenParticleRef> & zGenMatches = (iz == 0 ? z1GenMatches_ : z2GenMatches_);
    return (unsigned(il) < zGenMatches.size() ? zGenMatches[il] : reco::GenParticleRef());
}


const int reco::SkimEvent4L::channel() const { 
    switch(hypo_) {
        case ZZ4MU:    return 0;
        case ZZ4EL:    return 1;
        case ZZ2EL2MU: return 2;
        case ZZ2MU2EL: return 3;
        default: return -1;
    }
}

const float reco::SkimEvent4L::elePtMin() const 
{
    switch(hypo_) {
        case ZZ4MU:    return 999;
        case ZZ4EL:    return min(min(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(1)->pt()), min(daughter(1)->daughter(0)->pt(), daughter(1)->daughter(1)->pt()));
        case ZZ2EL2MU: return min(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(1)->pt());
        case ZZ2MU2EL: return min(daughter(1)->daughter(0)->pt(), daughter(1)->daughter(1)->pt());
        default: return 0;
    }
}

const float reco::SkimEvent4L::muPtMin() const 
{
    switch(hypo_) {
        case ZZ4EL:    return 999;
        case ZZ4MU:    return min(min(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(1)->pt()), min(daughter(1)->daughter(0)->pt(), daughter(1)->daughter(1)->pt()));
        case ZZ2MU2EL: return min(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(1)->pt());
        case ZZ2EL2MU: return min(daughter(1)->daughter(0)->pt(), daughter(1)->daughter(1)->pt());
        default: return 0;
    }
}

const float reco::SkimEvent4L::lPtMin() const 
{
    return min(min(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(1)->pt()), min(daughter(1)->daughter(0)->pt(), daughter(1)->daughter(1)->pt()));
}

const float reco::SkimEvent4L::eleEtaMax() const 
{
    switch(hypo_) {
        case ZZ4MU:    return 0;
        case ZZ4EL:    return max(max(fabs(daughter(0)->daughter(0)->eta()), fabs(daughter(0)->daughter(1)->eta())), max(fabs(daughter(1)->daughter(0)->eta()), fabs(daughter(1)->daughter(1)->eta())));
        case ZZ2EL2MU: return max(fabs(daughter(0)->daughter(0)->eta()), fabs(daughter(0)->daughter(1)->eta()));
        case ZZ2MU2EL: return max(fabs(daughter(1)->daughter(0)->eta()), fabs(daughter(1)->daughter(1)->eta()));
        default: return 999;
    }
}

const float reco::SkimEvent4L::muEtaMax() const 
{
    switch(hypo_) {
        case ZZ4EL:    return 0;
        case ZZ4MU:    return max(max(fabs(daughter(0)->daughter(0)->eta()), fabs(daughter(0)->daughter(1)->eta())), max(fabs(daughter(1)->daughter(0)->eta()), fabs(daughter(1)->daughter(1)->eta())));
        case ZZ2MU2EL: return max(fabs(daughter(0)->daughter(0)->eta()), fabs(daughter(0)->daughter(1)->eta()));
        case ZZ2EL2MU: return max(fabs(daughter(1)->daughter(0)->eta()), fabs(daughter(1)->daughter(1)->eta()));
        default: return 999;
    }
}

const float reco::SkimEvent4L::lEtaMax() const 
{
    return max(max(fabs(daughter(0)->daughter(0)->eta()), fabs(daughter(0)->daughter(1)->eta())), max(fabs(daughter(1)->daughter(0)->eta()), fabs(daughter(1)->daughter(1)->eta())));
}

const int reco::SkimEvent4L::nGoodLeptons(const char *muId, const char *eleId) const 
{
    return nGoodLeptons(std::string(muId), std::string(eleId));
}

const int reco::SkimEvent4L::nGoodLeptons(const std::string &muId, const std::string &eleId) const 
{
     int ngood = 0;

     std::map<std::string,utils::Selector>::const_iterator itMuId = utils::selectorCache_.find(muId);
     if (itMuId == utils::selectorCache_.end()) {
        itMuId = utils::selectorCache_.insert(std::make_pair(muId, utils::Selector(muId))).first;
     }
     std::map<std::string,utils::Selector>::const_iterator itElId = utils::selectorCache_.find(eleId);
     if (itElId == utils::selectorCache_.end()) {
        itElId = utils::selectorCache_.insert(std::make_pair(eleId, utils::Selector(eleId))).first;
     }

     for (int i = 0; i < 4; ++i) {
         const reco::Candidate &c = l(i/2,i%2);
         if      (typeid(c) == typeid(pat::Muon))     ngood += itMuId->second(c);
         else if (typeid(c) == typeid(pat::Electron)) ngood += itElId->second(c);
         else throw cms::Exception("WrongType") << "Lepton " << (i/2) << ", " << (i%2) << " is of type " << typeid(c).name() << "\n";
    }
}


