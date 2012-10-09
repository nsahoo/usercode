#ifndef YIELDMAKER_H
#define YIELDMAKER_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <TH1D.h>
#include <TH2D.h>

#include "FitMaker.h"
#include "FakeRateCalculator.h"

using namespace RooFit;


class YieldMaker {

    protected :

        RooRealVar rrchannel  ;
        RooRealVar rrz1mass   ;
        RooRealVar rrz2mass   ;
        RooRealVar rrmass     ;
        RooRealVar rrmela     ;
        RooRealVar rrweight   ;
        RooRealVar rrweighterr;
        RooArgSet argset      ;
        RooDataSet dataset    ;

    public :        

        YieldMaker():
            rrchannel  (RooRealVar("channel",   "channel",   0., 10.)), 
            rrz1mass   (RooRealVar("z1mass",    "z1mass",    0., 2000000.)), 
            rrz2mass   (RooRealVar("z2mass",    "z2mass",    0., 2000000.)), 
            rrmass     (RooRealVar("mass",      "mass",      0., 10000000.)),
            rrmela     (RooRealVar("mela",      "mela",      0., 1.)),
            rrweight   (RooRealVar("weight",    "weight",    -10., 10.)),
            rrweighterr(RooRealVar("weighterr", "weighterr", 0., 10000000.)),
            argset(RooArgSet(rrchannel, rrz1mass, rrz2mass, rrmass, rrmela, rrweight, rrweighterr, "argset")),
            dataset(RooDataSet("dataset", "dataset", argset))
        {}

        float getYield(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut) {

            float yield = 0.0;

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass = dataset.get(i)->getRealValue("z1mass");
                float z2mass = dataset.get(i)->getRealValue("z2mass");
                float mass   = dataset.get(i)->getRealValue("mass");
                float mela   = dataset.get(i)->getRealValue("mela");
                float weight = dataset.get(i)->getRealValue("weight");
                float ch     = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) yield += weight;
            }

            return yield;

        }

        float getYieldError(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut) {

            float yield    = 0.0;
            float yield_up = 0.0;
            float yield_dn = 0.0;

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float weight    = dataset.get(i)->getRealValue("weight");
                float weighterr = dataset.get(i)->getRealValue("weighterr");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                    yield    += weight;
                    yield_up += weight + weighterr;
                    yield_dn += weight - weighterr;
                }
            }

            return std::max<float>(yield_up - yield, yield - yield_dn);


        }

        float getYieldStatError(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut) {

            float yielderr = 0.0; 

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                    yielderr += weight*weight;
                }
            }

            return sqrt(yielderr);


        }

        float getCount(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut) {

            float yield = 0.0;

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass = dataset.get(i)->getRealValue("z1mass");
                float z2mass = dataset.get(i)->getRealValue("z2mass");
                float mass   = dataset.get(i)->getRealValue("mass");
                float mela   = dataset.get(i)->getRealValue("mela");
                float ch     = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) yield += 1.0;
            }

            return yield;

        }

        void get1DHist(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, TH1* hist) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(mass, weight);
            }

        }

        void get2DHist(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, TH2* hist) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(mass, mela, weight);
            }

        }

        RooDataSet getFitDataSet(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut) {
            RooRealVar m("mass",   "mass",   100, 100, 1000, "GeV/c^{2}");
            RooRealVar w("weight", "weight", 0.,  -10.,  10.);
            RooArgSet aset(m, w, "aset");
            RooDataSet dset("dataset","", aset);

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                    aset.setRealValue("mass", mass);
                    aset.setRealValue("weight", weight);
                    dset.add(aset);
                }
            }

            return dset;

        }

};


class DataYieldMaker : public YieldMaker {

    private :
        std::vector<std::pair<int, int> > runeventinfo;


    public:

        DataYieldMaker():YieldMaker(){}

        void fill(std::string filepath) {
            TFile* file = new TFile(filepath.c_str());
            TTree* tree = (TTree*)file->Get("zz4lTree/probe_tree");

            TBranch *bchannel   = tree->GetBranch("channel");
            TBranch *bz1mass    = tree->GetBranch("z1mass");
            TBranch *bz2mass    = tree->GetBranch("z2mass");
            TBranch *bmass      = tree->GetBranch("mass");
            TBranch *bmela      = tree->GetBranch("melaLD");
            TBranch *bevent     = tree->GetBranch("event");
            TBranch *brun       = tree->GetBranch("run");

            float channel   = 0.0;
            float z1mass    = 0.0;
            float z2mass    = 0.0;
            float mass      = 0.0;
            float mela      = 0.0;
            int   event     = 0;
            int   run       = 0;

            bchannel   ->SetAddress(&channel);
            bz1mass    ->SetAddress(&z1mass);
            bz2mass    ->SetAddress(&z2mass);
            bmass      ->SetAddress(&mass);
            bmela      ->SetAddress(&mela);
            bevent     ->SetAddress(&event);
            brun       ->SetAddress(&run);

            for (int i = 0; i < tree->GetEntries(); i++) {
                bchannel   ->GetEvent(i);
                bz1mass    ->GetEvent(i);
                bz2mass    ->GetEvent(i);
                bmass      ->GetEvent(i);
                bevent     ->GetEvent(i);
                brun       ->GetEvent(i);
                bmela      ->GetEvent(i);

                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
                }

                if (!existsAlready) {
                    argset.setRealValue("z1mass",    z1mass);
                    argset.setRealValue("z2mass",    z2mass);
                    argset.setRealValue("mass",      mass);
                    argset.setRealValue("mela",      mela);
                    argset.setRealValue("channel",   channel);
                    argset.setRealValue("weight",    1.0);
                    argset.setRealValue("weighterr", 0.0);
                    runeventinfo.push_back(std::pair<int, int>(run, event));
                    if (mass>140. && mass<300.) dataset.add(argset);
                }
            }
        }

        void getDataSet1D(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, RooDataSet& dset, RooRealVar& m) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                    m.setVal(mass);
                    RooArgSet aset(m, "argset_obs");
                    dset.add(aset);
                }
            }
        }


        void getDataSet2D(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, RooDataSet& dset, RooRealVar& m, RooRealVar& D) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                    m.setVal(mass);
                    D.setVal(mela);
                    RooArgSet aset(m, D, "argset_obs");
                    dset.add(aset);
                }
            }
        }
};


class ZXYieldMaker : public YieldMaker {

    private :
        std::vector<std::pair<int, int> > runeventinfo;    

        RooDataSet data2p2f;
        RooDataSet data3p2f;

        

    public:
       
        ZXYieldMaker():YieldMaker(){}
 
        void fill(std::string filepath, float wgt, FakeRateCalculator FR, bool doSS, bool isMC=false, bool doZXWgt=true, int PUWgtMode=1) {
        
            TFile* file = new TFile(filepath.c_str());
            TTree* tree = (TTree*)file->Get("zxTree/probe_tree");
            
            TBranch *bnumsim;
            if (isMC)bnumsim    = tree->GetBranch("numTrueInteractions"); 
            TBranch *bchannel   = tree->GetBranch("channel");
            TBranch *bz1mass    = tree->GetBranch("z1mass");
            TBranch *bz2mass    = tree->GetBranch("z2mass");
            TBranch *bmass      = tree->GetBranch("mass");
            TBranch *bl3pt      = tree->GetBranch("l3pt");
            TBranch *bl3eta     = tree->GetBranch("l3eta");
            TBranch *bl3phi     = tree->GetBranch("l3phi");
            TBranch *bl3id      = tree->GetBranch("l3idNew");
            TBranch *bl3iso     = tree->GetBranch("l3pfIsoComb04EACorr");
            TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
            TBranch *bl4pt      = tree->GetBranch("l4pt");
            TBranch *bl4eta     = tree->GetBranch("l4eta");
            TBranch *bl4phi     = tree->GetBranch("l4phi");
            TBranch *bl4id      = tree->GetBranch("l4idNew");
            TBranch *bl4iso     = tree->GetBranch("l4pfIsoComb04EACorr");
            TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
            TBranch *bl1pt      = tree->GetBranch("l1pt");
            TBranch *bl1eta     = tree->GetBranch("l1eta");
            TBranch *bl1phi     = tree->GetBranch("l1phi");
            TBranch *bl1pdgId   = tree->GetBranch("l1pdgId");
            TBranch *bl2pt      = tree->GetBranch("l2pt");
            TBranch *bl2eta     = tree->GetBranch("l2eta");
            TBranch *bl2phi     = tree->GetBranch("l2phi");
            TBranch *bl2pdgId   = tree->GetBranch("l2pdgId");
            TBranch *bmela      = tree->GetBranch("melaLD");
            TBranch *bevent     = tree->GetBranch("event");
            TBranch *brun       = tree->GetBranch("run");
           
            float numsim    = 0.0; 
            float channel   = 0.0;
            float z1mass    = 0.0;
            float z2mass    = 0.0;
            float mass      = 0.0;
            float l3pt      = 0.0;
            float l3eta     = 0.0;
            float l3phi     = 0.0;
            float l3id      = 0.0;
            float l3iso     = 0.0;
            float l3pdgId   = 0.0;
            float l4pt      = 0.0;
            float l4eta     = 0.0;
            float l4phi     = 0.0;
            float l4id      = 0.0;
            float l4iso     = 0.0;
            float l4pdgId   = 0.0;
            float l1pt      = 0.0;
            float l1eta     = 0.0;
            float l1phi     = 0.0;
            float l1pdgId   = 0.0;
            float l2pt      = 0.0;
            float l2eta     = 0.0;
            float l2phi     = 0.0;
            float l2pdgId   = 0.0;
            float mela      = 0.0;
            int   event     = 0;
            int   run       = 0;
       
            if (isMC) bnumsim    ->SetAddress(&numsim); 
            bchannel   ->SetAddress(&channel); 
            bz1mass    ->SetAddress(&z1mass);
            bz2mass    ->SetAddress(&z2mass);
            bmass      ->SetAddress(&mass);
            bl3pt      ->SetAddress(&l3pt);
            bl3eta     ->SetAddress(&l3eta);
            bl3phi     ->SetAddress(&l3phi);
            bl3id      ->SetAddress(&l3id);
            bl3iso     ->SetAddress(&l3iso);
            bl3pdgId   ->SetAddress(&l3pdgId);
            bl4pt      ->SetAddress(&l4pt);
            bl4eta     ->SetAddress(&l4eta);
            bl4phi     ->SetAddress(&l4phi);
            bl4id      ->SetAddress(&l4id);
            bl4iso     ->SetAddress(&l4iso);
            bl4pdgId   ->SetAddress(&l4pdgId);
            bl1pt      ->SetAddress(&l1pt);
            bl1eta     ->SetAddress(&l1eta);
            bl1phi     ->SetAddress(&l1phi);
            bl1pdgId   ->SetAddress(&l1pdgId);
            bl2pt      ->SetAddress(&l2pt);
            bl2eta     ->SetAddress(&l2eta);
            bl2phi     ->SetAddress(&l2phi);
            bl2pdgId   ->SetAddress(&l2pdgId);
            bmela      ->SetAddress(&mela);
            bevent     ->SetAddress(&event);
            brun       ->SetAddress(&run);
 
            for (int i = 0; i < tree->GetEntries(); i++) {
                if (isMC) bnumsim    ->GetEvent(i);
                bchannel   ->GetEvent(i);
                bz1mass    ->GetEvent(i);
                bz2mass    ->GetEvent(i);
                bmass      ->GetEvent(i);
                bl3pt      ->GetEvent(i);
                bl3eta     ->GetEvent(i);
                bl3phi     ->GetEvent(i);
                bl3id      ->GetEvent(i);
                bl3iso     ->GetEvent(i);
                bl3pdgId   ->GetEvent(i);
                bl4pt      ->GetEvent(i);
                bl4eta     ->GetEvent(i);
                bl4phi     ->GetEvent(i);
                bl4id      ->GetEvent(i);
                bl4iso     ->GetEvent(i);
                bl4pdgId   ->GetEvent(i);
                bevent     ->GetEvent(i);
                brun       ->GetEvent(i);
                bl1pt      ->GetEvent(i);
                bl1eta     ->GetEvent(i);
                bl1phi     ->GetEvent(i);
                bl1pdgId   ->GetEvent(i);
                bl2pt      ->GetEvent(i);
                bl2eta     ->GetEvent(i);
                bl2phi     ->GetEvent(i);
                bl2pdgId   ->GetEvent(i);
                bmela      ->GetEvent(i);
       
 
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
                }
        
                if (!existsAlready) {
                    argset.setRealValue("z1mass", z1mass);
                    argset.setRealValue("z2mass", z2mass);
                    argset.setRealValue("mass",   mass);
                    argset.setRealValue("mela",   mela);
                    argset.setRealValue("channel",channel);
                    runeventinfo.push_back(std::pair<int, int>(run, event));

                    float weight = wgt;
                    if (isMC) {
                        weight *= getPUWeight((int)numsim, PUWgtMode);
                        weight *= getSF(l1pt, l1eta, l1pdgId);
                        weight *= getSF(l2pt, l2eta, l2pdgId);
                    }

                    if (!doSS && l3pdgId == -l4pdgId) { 
                        float f1    = FR.getFakeRate(l3pt, l3eta, l3pdgId);
                        float f2    = FR.getFakeRate(l4pt, l4eta, l4pdgId);
                        float p1    = FR.getPromptRate(l3pt, l3eta, l3pdgId);
                        float p2    = FR.getPromptRate(l4pt, l4eta, l4pdgId);

                        float eps1  = f1/(1.0-f1);
                        float eps2  = f2/(1.0-f2);
                        float eta1  = (1.0-p1)/p1;
                        float eta2  = (1.0-p2)/p2;
                        float deno = (1-(eps1*eta1))*(1-(eps2*eta2));

                        if (doZXWgt) {
                            if      ((l3id==0 || l3iso/l3pt>0.4) && (l4id==0 || l4iso/l4pt>0.4)) weight *= -eps1*eps2;
                            else if ((l3id==0 || l3iso/l3pt>0.4) && (l4id==1 && l4iso/l4pt<0.4)) weight *= eps1;
                            else if ((l3id==1 && l3iso/l3pt<0.4) && (l4id==0 || l4iso/l4pt>0.4)) weight *= eps2;
                            else if ((l3id==1 && l3iso/l3pt<0.4) && (l4id==1 && l4iso/l4pt<0.4)) weight *= (-eps1*eta1 - eps2*eta2 + eps1*eps2*eta1*eta2);
                            weight /= deno;
                        }

                        argset.setRealValue("weight", weight);
                        argset.setRealValue("weighterr", 0.0);
                        dataset.add(argset);

                    }


                    else if (doSS && l3pdgId == l4pdgId) {
                        float f1    = FR.getFakeRate(l3pt, l3eta, l3pdgId);
                        float f2    = FR.getFakeRate(l4pt, l4eta, l4pdgId);
                        float f1_up = FR.getFakeRate(l3pt, l3eta, l3pdgId) + FR.getFakeRateErr(l3pt, l3eta, l3pdgId);
                        float f2_up = FR.getFakeRate(l4pt, l4eta, l4pdgId) + FR.getFakeRateErr(l4pt, l4eta, l4pdgId);
                        float f1_dn = FR.getFakeRate(l3pt, l3eta, l3pdgId) - FR.getFakeRateErr(l3pt, l3eta, l3pdgId);
                        float f2_dn = FR.getFakeRate(l4pt, l4eta, l4pdgId) - FR.getFakeRateErr(l4pt, l4eta, l4pdgId);
                        float sf = f1*f2;
                        float osss_mm = 1.30;
                        float osss_ee = 1.04;
                        float osss_em = 0.96;
                        if (channel == 0) sf*= osss_mm;
                        if (channel == 1) sf*= osss_ee;
                        if (channel == 2 || channel == 3) sf*= osss_em;

                        float sf_up = f1_up*f2_up;
                        if (channel == 0) sf_up*= osss_mm;
                        if (channel == 1) sf_up*= osss_ee;
                        if (channel == 2 || channel == 3) sf_up*= osss_em;

                        float sf_dn = f1_dn*f2_dn;
                        if (channel == 0) sf_dn*= osss_mm;
                        if (channel == 1) sf_dn*= osss_ee;
                        if (channel == 2 || channel == 3) sf_dn*= osss_em;

                        float weight_err = 0.0;
                        if (doZXWgt) weight *= sf;
                        if (doZXWgt) {
                            float weight_err_up = fabs(((sf != 0.0) ? weight * sf_up / sf : weight) - weight);
                            float weight_err_dn = fabs(((sf != 0.0) ? weight * sf_dn / sf : weight) - weight);
                            weight_err = (weight_err_up > weight_err_dn) ? weight_err_up : weight_err_dn; 
                        }

                        argset.setRealValue("weight", weight);
                        argset.setRealValue("weighterr", weight_err);
                        dataset.add(argset);

                    }

                }
            }
        }
};


class ZZYieldMaker : public YieldMaker {

    private:
        struct RunLumiEventInfo {
            unsigned run;
            unsigned lumi;
            unsigned event;
        };

        std::vector<RunLumiEventInfo> runeventinfo;    

    public:
        
        ZZYieldMaker():YieldMaker(){}

        void fill(std::string filepath, float wgt, float wgterr, bool isSignal, int PUWgtMode=1) {
            if (runeventinfo.size()>0) runeventinfo.clear();       

            TFile* file = new TFile(filepath.c_str());
            TTree* tree = (TTree*)file->Get("zz4lTree/probe_tree");
            
            TBranch *bchannel   = tree->GetBranch("channel");
            TBranch *bz1mass    = tree->GetBranch("z1mass");
            TBranch *bz2mass    = tree->GetBranch("z2mass");
            TBranch *bmass      = tree->GetBranch("mass");
            TBranch *bnumsim    = tree->GetBranch("numTrueInteractions");
            TBranch *bl3pt      = tree->GetBranch("l3pt");
            TBranch *bl3eta     = tree->GetBranch("l3eta");
            TBranch *bl3id      = tree->GetBranch("l3idNew");
            TBranch *bl3iso     = tree->GetBranch("l3pfIsoComb04EACorr");
            TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
            TBranch *bl4pt      = tree->GetBranch("l4pt");
            TBranch *bl4eta     = tree->GetBranch("l4eta");
            TBranch *bl4id      = tree->GetBranch("l4idNew");
            TBranch *bl4iso     = tree->GetBranch("l4pfIsoComb04EACorr");
            TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
            TBranch *bl1pt      = tree->GetBranch("l1pt");
            TBranch *bl1eta     = tree->GetBranch("l1eta");
            TBranch *bl1pdgId   = tree->GetBranch("l1pdgId");
            TBranch *bl2pt      = tree->GetBranch("l2pt");
            TBranch *bl2eta     = tree->GetBranch("l2eta");
            TBranch *bl2pdgId   = tree->GetBranch("l2pdgId");
            TBranch *bmela      = tree->GetBranch("melaLD");
            TBranch *bevent     = tree->GetBranch("event");
            TBranch *brun       = tree->GetBranch("run");
            TBranch *blumi      = tree->GetBranch("lumi");
            TBranch *bhiggswgt  = tree->GetBranch("genhiggsmassweight");
            
            float channel   = 0.0;
            float z1mass    = 0.0;
            float z2mass    = 0.0;
            float mass      = 0.0;
            float numsim    = 0.0;
            float l3pt      = 0.0;
            float l3eta     = 0.0;
            float l3id      = 0.0;
            float l3iso     = 0.0;
            float l3pdgId   = 0.0;
            float l4pt      = 0.0;
            float l4eta     = 0.0;
            float l4id      = 0.0;
            float l4iso     = 0.0;
            float l4pdgId   = 0.0;
            float l1pt      = 0.0;
            float l1eta     = 0.0;
            float l1pdgId   = 0.0;
            float l2pt      = 0.0;
            float l2eta     = 0.0;
            float l2pdgId   = 0.0;
            float mela      = 0.0;
            float higgswgt  = 0.0;
            unsigned event  = 0;
            unsigned run    = 0;
            unsigned lumi   = 0;

            bchannel   ->SetAddress(&channel);
            bz1mass    ->SetAddress(&z1mass);
            bz2mass    ->SetAddress(&z2mass);
            bmass      ->SetAddress(&mass);
            bnumsim    ->SetAddress(&numsim);
            bl3pt      ->SetAddress(&l3pt);
            bl3eta     ->SetAddress(&l3eta);
            bl3id      ->SetAddress(&l3id);
            bl3iso     ->SetAddress(&l3iso);
            bl3pdgId   ->SetAddress(&l3pdgId);
            bl4pt      ->SetAddress(&l4pt);
            bl4eta     ->SetAddress(&l4eta);
            bl4id      ->SetAddress(&l4id);
            bl4iso     ->SetAddress(&l4iso);
            bl4pdgId   ->SetAddress(&l4pdgId);
            bl1pt      ->SetAddress(&l1pt);
            bl1eta     ->SetAddress(&l1eta);
            bl1pdgId   ->SetAddress(&l1pdgId);
            bl2pt      ->SetAddress(&l2pt);
            bl2eta     ->SetAddress(&l2eta);
            bl2pdgId   ->SetAddress(&l2pdgId);
            bmela      ->SetAddress(&mela);
            bevent     ->SetAddress(&event);
            brun       ->SetAddress(&run);
            blumi      ->SetAddress(&lumi);
            bhiggswgt  ->SetAddress(&higgswgt);
        
            for (int i = 0; i < tree->GetEntries(); i++) {
                bchannel   ->GetEvent(i);
                bz1mass    ->GetEvent(i);
                bz2mass    ->GetEvent(i);
                bmass      ->GetEvent(i);
                bnumsim    ->GetEvent(i);
                bl3pt      ->GetEvent(i);
                bl3eta     ->GetEvent(i);
                bl3id      ->GetEvent(i);
                bl3iso     ->GetEvent(i);
                bl3pdgId   ->GetEvent(i);
                bl4pt      ->GetEvent(i);
                bl4eta     ->GetEvent(i);
                bl4id      ->GetEvent(i);
                bl4iso     ->GetEvent(i);
                bl4pdgId   ->GetEvent(i);
                bevent     ->GetEvent(i);
                brun       ->GetEvent(i);
                blumi      ->GetEvent(i);
                bl1pt      ->GetEvent(i);
                bl1eta     ->GetEvent(i);
                bl1pdgId   ->GetEvent(i);
                bl2pt      ->GetEvent(i);
                bl2eta     ->GetEvent(i);
                bl2pdgId   ->GetEvent(i);
                bmela      ->GetEvent(i);
                bhiggswgt  ->GetEvent(i);
        
                bool existsAlready = false;
       
                if (existsAlready) {
                    std::cout << "Run : " << run << " Lumi : " << lumi << " Event : " << event << std::endl;
                } 
 
                if (!existsAlready) {
                    argset.setRealValue("z1mass", z1mass);
                    argset.setRealValue("z2mass", z2mass);
                    argset.setRealValue("mass",   mass);
                    argset.setRealValue("mela",   mela);
                    argset.setRealValue("channel",channel);

                    RunLumiEventInfo rlei;
                    rlei.run = run;
                    rlei.lumi = lumi;
                    rlei.event = event;
                    runeventinfo.push_back(rlei);

                    float weight = wgt;
                    float weighterr = wgterr;

                    weight  *= getPUWeight((int)numsim, PUWgtMode);
                    weighterr  *= getPUWeight((int)numsim, PUWgtMode);
                    weight *= getSF(l1pt, l1eta, l1pdgId);
                    weight *= getSF(l2pt, l2eta, l2pdgId);
                    weight *= getSF(l3pt, l3eta, l3pdgId);
                    weight *= getSF(l4pt, l4eta, l4pdgId);

                    weighterr *= getSF(l1pt, l1eta, l1pdgId);
                    weighterr *= getSF(l2pt, l2eta, l2pdgId);
                    weighterr *= getSF(l3pt, l3eta, l3pdgId);
                    weighterr *= getSF(l4pt, l4eta, l4pdgId);

                    if (isSignal) weight    *= higgswgt;
                    if (isSignal) weighterr *= higgswgt;

                    argset.setRealValue("weight", weight);
                    argset.setRealValue("weighterr", weighterr);
                    dataset.add(argset);
                }
            }
        }

};

#endif
