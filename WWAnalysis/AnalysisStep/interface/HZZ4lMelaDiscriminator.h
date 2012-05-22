#ifndef HZZ4LMelaDiscriminator_h
#define HZZ4LMelaDiscriminator_h
#include <cmath>
class AngularPdfFactory;
class RooRealVar; 
class TH1F; class TH2F; class TH3F; class TDirectory;
class RooqqZZ_JHU;

class HZZ4LMelaDiscriminator {
    public:
        enum Model { SMHiggs, PSHiggs, QQZZ };
        HZZ4LMelaDiscriminator(Model model, const char *qqZZTemplateFileName) ;
        ~HZZ4LMelaDiscriminator() ;
        double get(double m4l, double mz1, double mz2, double cosThetaStar, double cosTheta1, double cosTheta2, double phi, double phi1) const ;
    private:
        struct QQZZ8DHistos {
            public:
              QQZZ8DHistos(TDirectory *dir) ;
              ~QQZZ8DHistos() ;
              double get(bool normalized,double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1) const ;
            private:
              TH1F *h_mzz;
              TH3F *h_mzzm1m2;
              TH2F *h_mzzcosthetastar;
              TH2F *h_mzzcostheta1;
              TH2F *h_mzzcostheta2;
              TH2F *h_mzzphi1;
              TH2F *h_mzzphi;
        };

        Model model_;
        RooRealVar *mzz_, *z1mass_, *z2mass_, *cosThetaStar_, *cosTheta1_, *cosTheta2_, *phi_, *phi1_;
        AngularPdfFactory *APF_;
        QQZZ8DHistos      *QQZZH_;
        RooqqZZ_JHU       *QQZZP_;

        void checkZorder(double& z1mass, double& z2mass,
                 double& costhetastar, double& costheta1,
                 double& costheta2, double& phi,
                 double& phistar1) const ;
};

#endif
