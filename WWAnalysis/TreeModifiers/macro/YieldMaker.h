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

#include "WWAnalysis/TreeModifiers/interface/HiggsMassWeightProvider.h"
#include "FitMaker.h"
#include "FakeRateCalculator.h"
#include "scales2.h"

using namespace RooFit;


class YieldMaker {

    protected :

        RooRealVar rrchannel            ;
        RooRealVar rrz1mass             ;
        RooRealVar rrz2mass             ;
        RooRealVar rrmass               ;
        RooRealVar rrmela               ;
        RooRealVar rrmelaPS             ;
        RooRealVar rrmelaSpinTwoMinimal ;
        RooRealVar rrweight             ;
        RooRealVar rrweighterr          ;
        RooRealVar rrnjets              ;
        RooRealVar rrjetdeta            ;
        RooRealVar rrjetmass            ;
        RooRealVar rrpt4l               ;
        RooArgSet  argset               ;
        RooDataSet dataset              ;

    public :        

        YieldMaker():
            rrchannel  (RooRealVar("channel",   "channel",   0., 10.)), 
            rrz1mass   (RooRealVar("z1mass",    "z1mass",    0., 2000000.)), 
            rrz2mass   (RooRealVar("z2mass",    "z2mass",    0., 2000000.)), 
            rrmass     (RooRealVar("mass",      "mass",      0., 10000000.)),
            rrmela     (RooRealVar("mela",      "mela",      0., 1.)),
            rrmelaPS   (RooRealVar("melaPS",    "melaPS",      0., 1.)),
            rrmelaSpinTwoMinimal  (RooRealVar("melaSpinTwoMinimal",      "melaSpinTwoMinimal",  0., 1.)),
            rrweight   (RooRealVar("weight",    "weight",    -10., 10.)),
            rrweighterr(RooRealVar("weighterr", "weighterr", 0., 10000000.)),
            rrnjets    (RooRealVar("njets",     "njets",     0., 10000.)),
            rrjetdeta  (RooRealVar("jetdeta",   "jetdeta",   0., 100.)),
            rrjetmass  (RooRealVar("jetmass",   "jetdeta",   0., 10000000.)),
            rrpt4l     (RooRealVar("pt4l",      "pt4l",      0., 10000000.)),
            argset(RooArgSet(rrchannel, rrz1mass, rrz2mass, rrmass, rrmela, rrmelaPS, rrmelaSpinTwoMinimal, rrweight, rrweighterr, "argset1"), RooArgSet(rrnjets, rrjetdeta, rrjetmass, rrpt4l, "argset2"), "argset"),
            dataset(RooDataSet("dataset", "dataset", argset))
        {}

        float getYield(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {


            float yield = 0.0;

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass = dataset.get(i)->getRealValue("z1mass");
                float z2mass = dataset.get(i)->getRealValue("z2mass");
                float mass   = dataset.get(i)->getRealValue("mass");
                float mela   = dataset.get(i)->getRealValue("mela");
                float weight = dataset.get(i)->getRealValue("weight");
                float ch     = dataset.get(i)->getRealValue("channel");
                float njets  = dataset.get(i)->getRealValue("njets");
                float pt4l   = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut && (pt4lup<0 || pt4l<pt4lup)) yield += weight;
            }

            return yield;

        }

        float getYieldError(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {

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
                float njets     = dataset.get(i)->getRealValue("njets");
                float pt4l      = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut && mela>melacut && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut  && (pt4lup<0 || pt4l<pt4lup)) {
                    yield    += weight;
                    yield_up += weight + weighterr;
                    yield_dn += weight - weighterr;
                }
            }

            return std::max<float>(yield_up - yield, yield - yield_dn);


        }

        float getYieldStatError(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {

            float yielderr = 0.0; 

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");
                float njets     = dataset.get(i)->getRealValue("njets");
                float pt4l      = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut && mela>melacut && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut  && (pt4lup<0 || pt4l<pt4lup)) {
                    yielderr += weight*weight;
                }
            }

            return sqrt(yielderr);


        }

        float getCount(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {

            float yield = 0.0;

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass = dataset.get(i)->getRealValue("z1mass");
                float z2mass = dataset.get(i)->getRealValue("z2mass");
                float mass   = dataset.get(i)->getRealValue("mass");
                float mela   = dataset.get(i)->getRealValue("mela");
                float ch     = dataset.get(i)->getRealValue("channel");
                float njets  = dataset.get(i)->getRealValue("njets");
                float pt4l   = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut  && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut && (pt4lup<0 || pt4l<pt4lup)) yield += 1.0;
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

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(mass, weight);
            }

        }

        void get2DHist(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, TH2* hist, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");
                float njets     = dataset.get(i)->getRealValue("njets");
                float pt4l      = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut && (pt4lup<0 || pt4l<pt4lup)) hist->Fill(mass, mela, weight);
            }

        }

        void fillMelaPSVsMass2DHist(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, TH2* hist) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float melaPS    = dataset.get(i)->getRealValue("melaPS");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                  hist->Fill(mass, melaPS, weight);
                }
            }

        }

        void fillMelaSpinTwoMinimalVsMass2DHist(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, TH2* hist) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float melaSpinTwoMinimal   = dataset.get(i)->getRealValue("melaSpinTwoMinimal");
                float weight    = dataset.get(i)->getRealValue("weight");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(mass, melaSpinTwoMinimal, weight);
            }

        }


        RooDataSet getFitDataSet(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {
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
                float njets     = dataset.get(i)->getRealValue("njets");
                float pt4l      = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut && (pt4lup<0 || pt4l<pt4lup)) {
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
            TFile* file = TFile::Open(filepath.c_str());
            TTree* tree = (TTree*)file->Get("zz4lTree/probe_tree");

            TBranch *bchannel              = tree->GetBranch("channel");
            TBranch *bz1mass               = tree->GetBranch("z1mass");
            TBranch *bz2mass               = tree->GetBranch("z2mass");
            TBranch *bmass                 = tree->GetBranch("mass");
            TBranch *bmela                 = tree->GetBranch("melaLD");
            TBranch *bmelaPS               = tree->GetBranch("melaPSLD");
            TBranch *bmelaSpinTwoMinimal   = tree->GetBranch("melaSpinTwoMinimal");
            TBranch *bnjets30              = tree->GetBranch("njets30");
            TBranch *bmjj                  = tree->GetBranch("mjj");
            TBranch *bj1eta                = tree->GetBranch("jet1eta");
            TBranch *bj2eta                = tree->GetBranch("jet2eta");
            TBranch *bl1pt                 = tree->GetBranch("l1pt");
            TBranch *bl2pt                 = tree->GetBranch("l2pt");
            TBranch *bl3pt                 = tree->GetBranch("l3pt");
            TBranch *bl4pt                 = tree->GetBranch("l4pt");
            TBranch *bl1eta                = tree->GetBranch("l1eta");
            TBranch *bl2eta                = tree->GetBranch("l2eta");
            TBranch *bl3eta                = tree->GetBranch("l3eta");
            TBranch *bl4eta                = tree->GetBranch("l4eta");
            TBranch *bl1phi                = tree->GetBranch("l1phi");
            TBranch *bl2phi                = tree->GetBranch("l2phi");
            TBranch *bl3phi                = tree->GetBranch("l3phi");
            TBranch *bl4phi                = tree->GetBranch("l4phi");
            TBranch *bevent                = tree->GetBranch("event");
            TBranch *brun                  = tree->GetBranch("run");

            float channel            = 0.0;
            float z1mass             = 0.0;
            float z2mass             = 0.0;
            float mass               = 0.0;
            float mela               = 0.0;
            float melaPS             = 0.0;
            float melaSpinTwoMinimal = 0.0;
            int   event              = 0;
            int   run                = 0;
            float njets30            = 0.0;
            float mjj                = 0.0;
            float j1eta              = 0.0;
            float j2eta              = 0.0;
            float l1pt               = 0.0;
            float l2pt               = 0.0;
            float l3pt               = 0.0;
            float l4pt               = 0.0;
            float l1eta              = 0.0;
            float l2eta              = 0.0;
            float l3eta              = 0.0;
            float l4eta              = 0.0;
            float l1phi              = 0.0;
            float l2phi              = 0.0;
            float l3phi              = 0.0;
            float l4phi              = 0.0;

            bchannel   ->SetAddress(&channel);
            bz1mass    ->SetAddress(&z1mass);
            bz2mass    ->SetAddress(&z2mass);
            bmass      ->SetAddress(&mass);
            bmela      ->SetAddress(&mela);
            bmelaPS    ->SetAddress(&melaPS);
            bmelaSpinTwoMinimal ->SetAddress(&melaSpinTwoMinimal);
            bevent     ->SetAddress(&event);
            brun       ->SetAddress(&run);
            bnjets30   ->SetAddress(&njets30);
            bmjj       ->SetAddress(&mjj);
            bj1eta     ->SetAddress(&j1eta);
            bj2eta     ->SetAddress(&j2eta);
            bl1pt      ->SetAddress(&l1pt);
            bl2pt      ->SetAddress(&l2pt);
            bl3pt      ->SetAddress(&l3pt);
            bl4pt      ->SetAddress(&l4pt);
            bl1eta     ->SetAddress(&l1eta);
            bl2eta     ->SetAddress(&l2eta);
            bl3eta     ->SetAddress(&l3eta);
            bl4eta     ->SetAddress(&l4eta);
            bl1phi     ->SetAddress(&l1phi);
            bl2phi     ->SetAddress(&l2phi);
            bl3phi     ->SetAddress(&l3phi);
            bl4phi     ->SetAddress(&l4phi);

            for (int i = 0; i < tree->GetEntries(); i++) {
                bchannel            ->GetEvent(i);
                bz1mass             ->GetEvent(i);
                bz2mass             ->GetEvent(i);
                bmass               ->GetEvent(i);
                bevent              ->GetEvent(i);
                brun                ->GetEvent(i);
                bmela               ->GetEvent(i);
                bmelaPS             ->GetEvent(i);
                bmelaSpinTwoMinimal ->GetEvent(i);
                bnjets30            ->GetEvent(i);
                bmjj                ->GetEvent(i);
                bj2eta              ->GetEvent(i);
                bj1eta              ->GetEvent(i);
                bl1pt               ->GetEvent(i);
                bl2pt               ->GetEvent(i);
                bl3pt               ->GetEvent(i);
                bl4pt               ->GetEvent(i);
                bl1eta              ->GetEvent(i);
                bl2eta              ->GetEvent(i);
                bl3eta              ->GetEvent(i);
                bl4eta              ->GetEvent(i);
                bl1phi              ->GetEvent(i);
                bl2phi              ->GetEvent(i);
                bl3phi              ->GetEvent(i);
                bl4phi              ->GetEvent(i);

                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
                }

                if (!existsAlready) {
                    argset.setRealValue("z1mass",             z1mass);
                    argset.setRealValue("z2mass",             z2mass);
                    argset.setRealValue("mass",               mass);
                    argset.setRealValue("mela",               mela);
                    argset.setRealValue("melaPS",             melaPS);
                    argset.setRealValue("melaSpinTwoMinimal", melaSpinTwoMinimal);
                    argset.setRealValue("channel",            channel);
                    argset.setRealValue("weight",             1.0);
                    argset.setRealValue("weighterr",          0.0);
                    argset.setRealValue("njets",              njets30);
                    argset.setRealValue("jetdeta",            fabs(j2eta-j1eta));
                    argset.setRealValue("jetmass",            mjj);

                    float p1x = l1pt*cos(l1phi);
                    float p2x = l2pt*cos(l2phi);
                    float p3x = l3pt*cos(l3phi);
                    float p4x = l4pt*cos(l4phi);

                    float p1y = l1pt*sin(l1phi);
                    float p2y = l2pt*sin(l2phi);
                    float p3y = l3pt*sin(l3phi);
                    float p4y = l4pt*sin(l4phi);

                    float px  = p1x+p2x+p3x+p4x;
                    float py  = p1y+p2y+p3y+p4y;
                    
                    argset.setRealValue("pt4l", sqrt(px*px + py*py));

                    //cout << "fill psmela: " << melaPS << " : " << mass << " : " << mela << endl;

                    runeventinfo.push_back(std::pair<int, int>(run, event));
                    //if (mass>140. && mass<300.) dataset.add(argset);
                    dataset.add(argset);
                }
            }
        }

        void getDataSet1D(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, RooDataSet& dset, RooRealVar& m, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float ch        = dataset.get(i)->getRealValue("channel");
                float njets     = dataset.get(i)->getRealValue("njets");
                float pt4l      = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut && (pt4lup<0 || pt4l<pt4lup)) {
                    m.setVal(mass);
                    RooArgSet aset(m, "argset_obs");
                    dset.add(aset);
                }
            }
        }


        void getDataSet2D(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, RooDataSet& dset, RooRealVar& m, RooRealVar& D, int njetcut=-1, bool veto=false, float pt4lcut=0.0, float pt4lup=-1) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float ch        = dataset.get(i)->getRealValue("channel");
                float njets     = dataset.get(i)->getRealValue("njets");
                float pt4l      = dataset.get(i)->getRealValue("pt4l");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut && (njetcut<0 || (!veto && njets==njetcut) || (veto && njets!=njetcut)) && pt4l>=pt4lcut && (pt4lup<0 || pt4l<pt4lup)) {
                    m.setVal(mass);
                    D.setVal(mela);
                    RooArgSet aset(m, D, "argset_obs");
                    dset.add(aset);
                }
            }
        }

        void getDataSetMassVsPSMela(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, RooDataSet& dset, RooRealVar& m, RooRealVar& D) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float melaPS    = dataset.get(i)->getRealValue("melaPS");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2mass<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                    m.setVal(mass);
                    D.setVal(melaPS);
                    RooArgSet aset(m, D, "argset_obs");
                    dset.add(aset);
                }
            }
        }

        void getDataSetMassVsSpinTwoMinimalMela(int channel, float z1min, float z2min, float m4lmin, float m4lmax, float melacut, RooDataSet& dset, RooRealVar& m, RooRealVar& D) {

            for (int i = 0; i < dataset.numEntries(); i++) {
                float z1mass    = dataset.get(i)->getRealValue("z1mass");
                float z2mass    = dataset.get(i)->getRealValue("z2mass");
                float mass      = dataset.get(i)->getRealValue("mass");
                float mela      = dataset.get(i)->getRealValue("mela");
                float melaSpinTwoMinimal  = dataset.get(i)->getRealValue("melaSpinTwoMinimal");
                float ch        = dataset.get(i)->getRealValue("channel");

                if (channel == (int)ch && z1mass>z1min && z1mass<120 && z2mass>z2min && z2min<120 && mass>m4lmin && mass<m4lmax && mela>melacut) {
                    m.setVal(mass);
                    D.setVal(melaSpinTwoMinimal);
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
        
            TFile* file = TFile::Open(filepath.c_str());
            TTree* tree = doSS ? (TTree*)file->Get("zxTreeSS/probe_tree") : (TTree*)file->Get("zxTreeOS/probe_tree");
            //TTree* tree = (TTree*)file->Get("zxTree/probe_tree");
            
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
            TBranch *bmelaPS               = tree->GetBranch("melaPSLD");
            TBranch *bmelaSpinTwoMinimal   = tree->GetBranch("melaSpinTwoMinimal");
            TBranch *bevent     = tree->GetBranch("event");
            TBranch *brun       = tree->GetBranch("run");
            TBranch *bnjets30   = tree->GetBranch("njets30");
            TBranch *bmjj       = tree->GetBranch("mjj");
            TBranch *bj1eta     = tree->GetBranch("jet1eta");
            TBranch *bj2eta     = tree->GetBranch("jet2eta");
           
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
            float njets30   = 0.0;
            float mjj       = 0.0;
            float j1eta     = 0.0;
            float j2eta     = 0.0;
            float mela      = 0.0;
            float melaPS             = 0.0;
            float melaSpinTwoMinimal = 0.0;
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
            bnjets30   ->SetAddress(&njets30);
            bmjj       ->SetAddress(&mjj);
            bj1eta     ->SetAddress(&j1eta);
            bj2eta     ->SetAddress(&j2eta);
            bmela      ->SetAddress(&mela);
            bmelaPS    ->SetAddress(&melaPS);
            bmelaSpinTwoMinimal ->SetAddress(&melaSpinTwoMinimal);
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
                bnjets30   ->GetEvent(i);
                bmjj       ->GetEvent(i);
                bj1eta     ->GetEvent(i);
                bj2eta     ->GetEvent(i);
                bmelaPS    ->GetEvent(i);
                bmelaSpinTwoMinimal ->GetEvent(i);

 
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
                }
        
                if (!existsAlready) {
                    argset.setRealValue("z1mass", z1mass);
                    argset.setRealValue("z2mass", z2mass);
                    argset.setRealValue("mass",   mass);
                    argset.setRealValue("mela",   mela);
                    argset.setRealValue("melaPS",             melaPS);
                    argset.setRealValue("melaSpinTwoMinimal", melaSpinTwoMinimal);
                    argset.setRealValue("channel",channel);
                    argset.setRealValue("njets",  njets30);
                    argset.setRealValue("jetdeta",fabs(j2eta-j1eta));
                    argset.setRealValue("jetmass",mjj);
                    runeventinfo.push_back(std::pair<int, int>(run, event));

                    float p1x = l1pt*cos(l1phi);
                    float p2x = l2pt*cos(l2phi);
                    float p3x = l3pt*cos(l3phi);
                    float p4x = l4pt*cos(l4phi);

                    float p1y = l1pt*sin(l1phi);
                    float p2y = l2pt*sin(l2phi);
                    float p3y = l3pt*sin(l3phi);
                    float p4y = l4pt*sin(l4phi);

                    float px  = p1x+p2x+p3x+p4x;
                    float py  = p1y+p2y+p3y+p4y;
                    
                    argset.setRealValue("pt4l", sqrt(px*px + py*py));


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

        void fill(std::string filepath, float wgt, float wgterr, bool isSignal, int PUWgtMode=1, int hmass=0, bool i7=true) {
            if (runeventinfo.size()>0) runeventinfo.clear();       

            TFile* file = TFile::Open(filepath.c_str());
            TTree* tree = (TTree*)file->Get("zz4lTree/probe_tree");
            
            TBranch *bchannel   = tree->GetBranch("channel");
            TBranch *bz1mass    = tree->GetBranch("z1mass");
            TBranch *bz2mass    = tree->GetBranch("z2mass");
            TBranch *bmass      = tree->GetBranch("mass");
            TBranch *bnumsim    = tree->GetBranch("numTrueInteractions");
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
            TBranch *bmelaPS               = tree->GetBranch("melaPSLD");
            TBranch *bmelaSpinTwoMinimal   = tree->GetBranch("melaSpinTwoMinimal");
            TBranch *bevent     = tree->GetBranch("event");
            TBranch *brun       = tree->GetBranch("run");
            TBranch *blumi      = tree->GetBranch("lumi");
            TBranch *bhiggswgt  = tree->GetBranch("genhiggsmass");
            TBranch *bnjets30   = tree->GetBranch("njets30");
            TBranch *bmjj       = tree->GetBranch("mjj");
            TBranch *bj1eta     = tree->GetBranch("jet1eta");
            TBranch *bj2eta     = tree->GetBranch("jet2eta");
            
            float channel   = 0.0;
            float z1mass    = 0.0;
            float z2mass    = 0.0;
            float mass      = 0.0;
            float numsim    = 0.0;
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
            float njets30   = 0.0;
            float mjj       = 0.0;
            float j1eta     = 0.0;
            float j2eta     = 0.0;
            float mela      = 0.0;
            float melaPS             = 0.0;
            float melaSpinTwoMinimal = 0.0;
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
            bmelaPS    ->SetAddress(&melaPS);
            bmelaSpinTwoMinimal ->SetAddress(&melaSpinTwoMinimal);
            bevent     ->SetAddress(&event);
            brun       ->SetAddress(&run);
            blumi      ->SetAddress(&lumi);
            bhiggswgt  ->SetAddress(&higgswgt);
            bnjets30   ->SetAddress(&njets30);
            bmjj       ->SetAddress(&mjj);
            bj1eta     ->SetAddress(&j1eta);
            bj2eta     ->SetAddress(&j2eta);

            std::stringstream weightss;
            if (isSignal) {
                if (hmass>=400 && (int(hmass))%50 == 0) {
                    if (i7)weightss <<"/home/avartak/CMS/Higgs/HCP/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs"<<hmass<<"_7TeV_Lineshape+Interference.txt";
                    else   weightss <<"/home/avartak/CMS/Higgs/HCP/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs"<<hmass<<"_8TeV_Lineshape+Interference.txt";
                }
                else weightss << "";
            }
            else weightss << "";
            HiggsMassWeightProvider hmwp(weightss.str());

 
            for (int i = 0; i < tree->GetEntries(); i++) {
                bchannel   ->GetEvent(i);
                bz1mass    ->GetEvent(i);
                bz2mass    ->GetEvent(i);
                bmass      ->GetEvent(i);
                bnumsim    ->GetEvent(i);
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
                blumi      ->GetEvent(i);
                bl1pt      ->GetEvent(i);
                bl1eta     ->GetEvent(i);
                bl1phi     ->GetEvent(i);
                bl1pdgId   ->GetEvent(i);
                bl2pt      ->GetEvent(i);
                bl2eta     ->GetEvent(i);
                bl2phi     ->GetEvent(i);
                bl2pdgId   ->GetEvent(i);
                bmela      ->GetEvent(i);
                bnjets30   ->GetEvent(i);
                bmjj       ->GetEvent(i);
                bj2eta     ->GetEvent(i);
                bj1eta     ->GetEvent(i);
                bmelaPS             ->GetEvent(i);
                bmelaSpinTwoMinimal ->GetEvent(i);
                bhiggswgt  ->GetEvent(i);
        
                bool existsAlready = false;
      
                /* 
                if (existsAlready) {
                    std::cout << "Run : " << run << " Lumi : " << lumi << " Event : " << event << std::endl;
                } 
                */
 
                if (!existsAlready) {
                    argset.setRealValue("z1mass", z1mass);
                    argset.setRealValue("z2mass", z2mass);
                    argset.setRealValue("mass",   mass);
                    argset.setRealValue("mela",   mela);
                    argset.setRealValue("melaPS", melaPS);
                    argset.setRealValue("melaSpinTwoMinimal", melaSpinTwoMinimal);
                    argset.setRealValue("channel",channel);
                    argset.setRealValue("njets",  njets30);
                    argset.setRealValue("jetdeta",fabs(j2eta-j1eta));
                    argset.setRealValue("jetmass",mjj);

                    float p1x = l1pt*cos(l1phi);
                    float p2x = l2pt*cos(l2phi);
                    float p3x = l3pt*cos(l3phi);
                    float p4x = l4pt*cos(l4phi);

                    float p1y = l1pt*sin(l1phi);
                    float p2y = l2pt*sin(l2phi);
                    float p3y = l3pt*sin(l3phi);
                    float p4y = l4pt*sin(l4phi);

                    float px  = p1x+p2x+p3x+p4x;
                    float py  = p1y+p2y+p3y+p4y;

                    argset.setRealValue("pt4l", sqrt(px*px + py*py));

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

                    if (isSignal) {
                        if (hmass>=400 && (int(hmass))%50 == 0) {
                            weight    *= hmwp.getWeight(higgswgt); 
                            weighterr *= hmwp.getWeight(higgswgt); 
                        }
                    }

                    argset.setRealValue("weight", weight);
                    argset.setRealValue("weighterr", weighterr);
                    dataset.add(argset);
                }
            }
        }

};

#endif
