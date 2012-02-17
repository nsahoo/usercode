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
    clearDaughters();
    const reco::Candidate *srcz1 = src.daughter(0);
    const reco::Candidate *srcz2 = src.daughter(1);
    if (fabs(srcz1->mass()-91.188) > fabs(srcz2->mass()-91.188)) std::swap(srcz1, srcz2);
    addDaughter(*srcz1);
    addDaughter(*srcz2);
    init();
}

void reco::SkimEvent4L::init() {
    if (numberOfDaughters() != 2) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: I should have two daughters!\n";
    if (daughter(0) == 0) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z1 is null\n";
    if (daughter(1) == 0) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z2 is null\n";
    if (daughter(0)->numberOfDaughters() < 2) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z1 should have at least two daughters!\n";
    if (daughter(1)->numberOfDaughters() < 2) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z2 should have at least two daughters!\n";
    if (lproxy(0,0) == 0) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z1 daughter 0 is null\n";
    if (lproxy(0,1) == 0) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z1 daughter 1 is null\n";
    if (lproxy(1,0) == 0) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z2 daughter 0 is null\n";
    if (lproxy(1,1) == 0) throw cms::Exception("CorruptData") << "reco::SkimEvent4L::init: Z2 daughter 1 is null\n";
    if (abs(lproxy(0,0)->pdgId()) == 11) {
        hypo_ = abs(lproxy(1,0)->pdgId()) == 11 ? ZZ4EL : ZZ2EL2MU;
    } else {
        hypo_ = abs(lproxy(1,0)->pdgId()) == 13 ? ZZ4MU : ZZ2MU2EL;
    }
}

void
reco::SkimEvent4L::setEventInfo(edm::Event& e)
{
    event_ = e.id().event();
    lumi_ = e.id().luminosityBlock();
    run_ = e.id().run();
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
    table_[0].first = -lproxy(0,0)->pt(); table_[0].second = 0;
    table_[1].first = -lproxy(0,1)->pt(); table_[1].second = 1;
    table_[2].first = -lproxy(1,0)->pt(); table_[2].second = 2;
    table_[3].first = -lproxy(1,1)->pt(); table_[3].second = 3;
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

float reco::SkimEvent4L::worsePairCombRelIsoBaseline() const {
    float vals[4];
    vals[0] = lisoCombRelBaseline(0,0);
    vals[1] = lisoCombRelBaseline(0,1);
    vals[2] = lisoCombRelBaseline(1,0);
    vals[3] = lisoCombRelBaseline(1,1);
    std::sort(&vals[0], &vals[4]);
    return vals[3]+vals[2];
}

float reco::SkimEvent4L::lisoPf(unsigned int iz, unsigned int il, const char *name) const {
    char buff[255];
    sprintf(buff,"pfiso4l%s_%d_%d", name, iz, il);
    return userFloat(buff);
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
        case ZZ4EL:    return min(min(lproxy(0,0)->pt(), lproxy(0,1)->pt()), min(lproxy(1,0)->pt(), lproxy(1,1)->pt()));
        case ZZ2EL2MU: return min(lproxy(0,0)->pt(), lproxy(0,1)->pt());
        case ZZ2MU2EL: return min(lproxy(1,0)->pt(), lproxy(1,1)->pt());
        default: return 0;
    }
}

const float reco::SkimEvent4L::muPtMin() const 
{
    switch(hypo_) {
        case ZZ4EL:    return 999;
        case ZZ4MU:    return min(min(lproxy(0,0)->pt(), lproxy(0,1)->pt()), min(lproxy(1,0)->pt(), lproxy(1,1)->pt()));
        case ZZ2MU2EL: return min(lproxy(0,0)->pt(), lproxy(0,1)->pt());
        case ZZ2EL2MU: return min(lproxy(1,0)->pt(), lproxy(1,1)->pt());
        default: return 0;
    }
}

const float reco::SkimEvent4L::lPtMin() const 
{
    return min(min(lproxy(0,0)->pt(), lproxy(0,1)->pt()), min(lproxy(1,0)->pt(), lproxy(1,1)->pt()));
}

const float reco::SkimEvent4L::eleEtaMax() const 
{
    switch(hypo_) {
        case ZZ4MU:    return 0;
        case ZZ4EL:    return max(max(fabs(lproxy(0,0)->eta()), fabs(lproxy(0,1)->eta())), max(fabs(lproxy(1,0)->eta()), fabs(lproxy(1,1)->eta())));
        case ZZ2EL2MU: return max(fabs(lproxy(0,0)->eta()), fabs(lproxy(0,1)->eta()));
        case ZZ2MU2EL: return max(fabs(lproxy(1,0)->eta()), fabs(lproxy(1,1)->eta()));
        default: return 999;
    }
}

const float reco::SkimEvent4L::muEtaMax() const 
{
    switch(hypo_) {
        case ZZ4EL:    return 0;
        case ZZ4MU:    return max(max(fabs(lproxy(0,0)->eta()), fabs(lproxy(0,1)->eta())), max(fabs(lproxy(1,0)->eta()), fabs(lproxy(1,1)->eta())));
        case ZZ2MU2EL: return max(fabs(lproxy(0,0)->eta()), fabs(lproxy(0,1)->eta()));
        case ZZ2EL2MU: return max(fabs(lproxy(1,0)->eta()), fabs(lproxy(1,1)->eta()));
        default: return 999;
    }
}

const float reco::SkimEvent4L::lEtaMax() const 
{
    return max(max(fabs(lproxy(0,0)->eta()), fabs(lproxy(0,1)->eta())), max(fabs(lproxy(1,0)->eta()), fabs(lproxy(1,1)->eta())));
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

    return ngood;
}

const int reco::SkimEvent4L::nGoodPairs(const char *pairCut, bool anySign) const {
    return nGoodPairs(std::string(pairCut), anySign);
}

const int reco::SkimEvent4L::nGoodPairs(const std::string &pairCut, int anySign) const {
    std::map<std::string,utils::Selector>::const_iterator itPairCut = utils::selectorCache_.find(pairCut);
    if (itPairCut == utils::selectorCache_.end()) {
        itPairCut = utils::selectorCache_.insert(std::make_pair(pairCut, utils::Selector(pairCut))).first;
    }
    int nGood = 0;
    Pair pair;
    for (int i = 0; i < 4; ++i) { 
        const reco::Candidate &li = l(i/2, i%2);
        for (int j = i+1; j < 4; ++j) {
            const reco::Candidate &lj = l(j/2, j%2);
            if (!anySign && (li.charge() + lj.charge() != 0)) continue;
            pair.clear();
            pair.addDaughter(&li); 
            pair.addDaughter(&lj); 
            if (itPairCut->second(pair)) nGood++;
        }
    }
    return nGood;
}

