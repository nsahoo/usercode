#ifndef CARDTEMPLETE_H
#define CARDTEMPLETE_H

float getQQZZPDFUncertainty7TeV(float m4l) {
    return 1. + 0.0035*sqrt(m4l-30.);
}

float getGGZZPDFUncertainty7TeV(float m4l) {
    return 1. + 0.0066*sqrt(m4l-10.);
}

float getQQZZQCDScaleUncertainty7TeV(float m4l) {
    return 1. + 0.01*sqrt((m4l-20.)/13.);
}

float getGGZZQCDScaleUncertainty7TeV(float m4l) {
    return 1.04 + 0.1*sqrt((m4l+40.)/40.);
}

float getGGHPDFUncertaintyUp7TeV(float mass) {
    if (mass >= 100. && mass < 105.) return 1.078;
    else if (mass >= 105. && mass < 111.) return 1.077;
    else if (mass >= 111. && mass < 116.) return 1.076;
    else if (mass >= 116. && mass < 133.) return 1.075;
    else if (mass >= 133. && mass < 141.) return 1.074;
    else if (mass >= 141. && mass < 149.) return 1.073;
    else if (mass >= 149. && mass < 153.) return 1.074;
    else if (mass >= 153. && mass < 174.) return 1.075;
    else if (mass >= 174. && mass < 220.) return 1.074;
    else if (mass >= 220                ) return 1.073;
    else return 1.0;
}

float getGGHPDFUncertaintyDown7TeV(float mass) {
    if (mass >= 100.  && mass < 106. ) return 0.935;
    else if (mass >= 106.  && mass < 109. ) return 0.934;
    else if (mass >= 109.  && mass < 113.5) return 0.933;
    else if (mass >= 113.5 && mass < 120. ) return 0.932;
    else if (mass >= 120.  && mass < 131. ) return 0.931;
    else if (mass >= 131.  && mass < 138. ) return 0.930;
    else if (mass >= 138.  && mass < 148. ) return 0.931;
    else if (mass >= 148.  && mass < 152. ) return 0.930;
    else if (mass >= 152.  && mass < 162. ) return 0.929;
    else if (mass >= 162.  && mass < 165. ) return 0.928;
    else if (mass >= 165.  && mass < 170. ) return 0.927;
    else if (mass >= 170.  && mass < 178. ) return 0.926;
    else if (mass >= 178.  && mass < 192. ) return 0.925;
    else if (mass >= 192.  && mass < 198. ) return 0.924;
    else if (mass >= 198.  && mass < 208. ) return 0.923;
    else if (mass >= 208.  && mass < 212. ) return 0.922;
    else if (mass >= 212.  && mass < 218. ) return 0.923;
    else if (mass >= 218.                 ) return 0.924;
    else return 1.0;
}

float getGGHQCDScaleUncertaintyUp7TeV(float mass) {
    if      (mass >= 100.  && mass < 103. ) return 1.079;
    else if (mass >= 103.  && mass < 106. ) return 1.078;
    else if (mass >= 106.  && mass < 109. ) return 1.077;
    else if (mass >= 109.  && mass < 111.5) return 1.076;
    else if (mass >= 111.5 && mass < 115. ) return 1.075;
    else if (mass >= 115.  && mass < 118.5) return 1.074;
    else if (mass >= 118.5 && mass < 122.5) return 1.073;
    else if (mass >= 122.5 && mass < 131.5) return 1.071;
    else if (mass >= 131.5 && mass < 136. ) return 1.070;
    else if (mass >= 136.  && mass < 141. ) return 1.069;
    else if (mass >= 141.  && mass < 146. ) return 1.068;
    else if (mass >= 146.  && mass < 152. ) return 1.067;
    else if (mass >= 152.  && mass < 158. ) return 1.066;
    else if (mass >= 158.  && mass < 164. ) return 1.065;
    else if (mass >= 164.  && mass < 172. ) return 1.064;
    else if (mass >= 172.  && mass < 180. ) return 1.063;
    else if (mass >= 180.  && mass < 186. ) return 1.062;
    else if (mass >= 186.  && mass < 196. ) return 1.061;
    else if (mass >= 196.  && mass < 216. ) return 1.060;
    else if (mass >= 216.                 ) return 1.059;
    else return 1.0;
}


float getGGHQCDScaleUncertaintyDown7TeV(float mass) {
    if      (mass >= 100.  && mass < 104. ) return 0.916;
    else if (mass >= 104.  && mass < 108. ) return 0.917;
    else if (mass >= 108.  && mass < 111.5) return 0.918;
    else if (mass >= 111.5 && mass < 115.5) return 0.918;
    else if (mass >= 115.5 && mass < 120. ) return 0.920;
    else if (mass >= 120.  && mass < 125. ) return 0.921;
    else if (mass >= 125.  && mass < 130. ) return 0.922;
    else if (mass >= 130.  && mass < 135.5) return 0.923;
    else if (mass >= 135.5 && mass < 142. ) return 0.924;
    else if (mass >= 142.  && mass < 148. ) return 0.925;
    else if (mass >= 148.  && mass < 154. ) return 0.926;
    else if (mass >= 154.  && mass < 162. ) return 0.927;
    else if (mass >= 162.  && mass < 170. ) return 0.928;
    else if (mass >= 170.  && mass < 178. ) return 0.929;
    else if (mass >= 178.  && mass < 186. ) return 0.930;
    else if (mass >= 186.  && mass < 196. ) return 0.931;
    else if (mass >= 196.  && mass < 206. ) return 0.932;
    else if (mass >= 206.  && mass < 216. ) return 0.933;
    else if (mass >= 216.                 ) return 0.934;
    else return 1.0;
}

float getVBFPDFUncertaintyUp7TeV(float mass) {
    if      (mass >= 100.  && mass < 110. ) return 1.022;
    else if (mass >= 110.  && mass < 115.5) return 1.023;
    else if (mass >= 115.5 && mass < 123.5) return 1.024;
    else if (mass >= 123.5 && mass < 134. ) return 1.025;
    else if (mass >= 134.  && mass < 142. ) return 1.026;
    else if (mass >= 142.  && mass < 152. ) return 1.027;
    else if (mass >= 152.  && mass < 162. ) return 1.028;
    else if (mass >= 162.  && mass < 170. ) return 1.029;
    else if (mass >= 170.  && mass < 180. ) return 1.030;
    else if (mass >= 180.  && mass < 188. ) return 1.031;
    else if (mass >= 188.  && mass < 198. ) return 1.032;
    else if (mass >= 198.  && mass < 208. ) return 1.033;
    else if (mass >= 208.  && mass < 216. ) return 1.034;
    else if (mass >= 220.                 ) return 1.035;
    else return 1.0;
}

float getVBFPDFUncertaintyDown7TeV(float mass) {
    if      (mass >= 100.  && mass < 184. ) return 0.979;
    else if (mass >= 184.                 ) return 0.980;
    else return 1.0;
}

float getVBFQCDScaleUncertaintyUp7TeV(float mass) {
    if      (mass >= 100.  && mass < 105. ) return 1.004;
    else if (mass >= 105.  && mass < 110. ) return 1.003;
    else if (mass >= 110.  && mass < 111. ) return 1.005;
    else if (mass >= 111.  && mass < 112.5) return 1.004;
    else if (mass >= 112.5 && mass < 114.5) return 1.003;
    else if (mass >= 114.5 && mass < 117. ) return 1.002;
    else if (mass >= 117.  && mass < 131.5) return 1.003;
    else if (mass >= 131.5 && mass < 137.5) return 1.004;
    else if (mass >= 131.5 && mass < 139. ) return 1.003;
    else if (mass >= 136.  && mass < 142. ) return 1.002;
    else if (mass >= 142.  && mass < 144. ) return 1.003;
    else if (mass >= 144.  && mass < 147. ) return 1.004;
    else if (mass >= 147.  && mass < 150. ) return 1.003;
    else if (mass >= 150.  && mass < 152. ) return 1.002;
    else if (mass >= 152.  && mass < 157. ) return 1.003;
    else if (mass >= 157.  && mass < 159. ) return 1.002;
    else if (mass >= 159.  && mass < 164. ) return 1.001;
    else if (mass >= 164.  && mass < 178. ) return 1.002;
    else if (mass >= 178.  && mass < 184. ) return 1.001;
    else if (mass >= 184.  && mass < 188. ) return 1.003;
    else if (mass >= 188.                 ) return 1.001;
    else return 1.0;
}


float getVBFQCDScaleUncertaintyDown7TeV(float mass) {
    if      (mass >= 100.  && mass < 110. ) return 0.997;
    else if (mass >= 110.  && mass < 116.5) return 0.998;
    else if (mass >= 116.5 && mass < 127. ) return 0.997;
    else if (mass >= 127.  && mass < 133. ) return 0.998;
    else if (mass >= 133.  && mass < 139.5) return 0.999;
    else if (mass >= 139.5 && mass < 141. ) return 0.998;
    else if (mass >= 141.  && mass < 159. ) return 0.999;
    else if (mass >= 159.  && mass < 164. ) return 0.998;
    else if (mass >= 164.  && mass < 178. ) return 0.999;
    else if (mass >= 178.  && mass < 186. ) return 0.998;
    else if (mass >= 186.  && mass < 188. ) return 0.999;
    else if (mass >= 188.  && mass < 192. ) return 0.998;
    else if (mass >= 192.  && mass < 194. ) return 0.997;
    else if (mass >= 194.  && mass < 198. ) return 0.996;
    else if (mass >= 198.  && mass < 200. ) return 0.997;
    else if (mass >= 200.  && mass < 202. ) return 0.998;
    else if (mass >= 202.  && mass < 212. ) return 0.997;
    else if (mass >= 212.                 ) return 0.996;
    else return 1.0;
}


std::string createCardTemplate(int channel, bool doShapeAnalysis, std::string workspacefilename) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
    if (doShapeAnalysis && channel == 0) {
        card += std::string("shapes * * ") + workspacefilename + " w_4mu:$PROCESS\n";
    }
    if (doShapeAnalysis && channel == 1) {
        card += std::string("shapes * * ") + workspacefilename + " w_4e:$PROCESS\n";
    }
    if (doShapeAnalysis && channel == 2) {
        card += std::string("shapes * * ") + workspacefilename + " w_2e2mu:$PROCESS\n";
    }
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN           BIN            BIN             BIN             BIN\n";
        card += "process SIG_GGH_NAME  SIG_VBF_NAME   BKG_QQZZ_NAME   BKG_GGZZ_NAME   BKG_ZJETS_NAME\n";
        card += "process -2            -1             1               2               3\n";
        card += "rate    SIG_GGH_YIELD SIG_VBF_YIELD  BKG_QQZZ_YIELD  BKG_GGZZ_YIELD  BKG_ZJETS_YIELD\n";
        card += "------------\n";
        card += "lumi_7TeV                 lnN        1.022   1.022   1.022    1.022    -\n";
        card += "pdf_gg                    lnN        GGH_PDF -       -        GGZZ_PDF -\n";
        card += "pdf_qqbar                 lnN         -      VBF_PDF QQZZ_PDF -        -\n";
        card += "pdf_hzz4l_accept          lnN        1.02    1.02    -        -        -\n";
        card += "QCDscale_ggH              lnN        GGH_QCD -       -        -        -\n";
        card += "QCDscale_qqH              lnN        -       VBF_PDF -        -        -\n";
        card += "QCDscale_ggVV             lnN        -       -       -        GGZZ_QCD -\n";
        card += "QCDscale_VV               lnN        -       -       QQZZ_QCD -        -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02    1.02    -        -        -\n";
    if (channel == 0) {
        card += "CMS_eff_m                 lnN        1.0224  1.0224  1.0224   1.0224   -\n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -       -        -        0.5/1.6\n";
    }
    else if (channel == 1) {
        card += "CMS_eff_e                 lnN        1.0447  1.0447  1.0447   1.0447   -\n";
        card += "CMS_hzz4e_Zjets           lnN        -       -       -        -        0.5/1.6\n";
    }
    else if (channel == 2) {
        card += "CMS_eff_m                 lnN        1.0224  1.0224  1.0224   1.0224   -\n";
        card += "CMS_eff_e                 lnN        1.0447  1.0447  1.0447   1.0447   -\n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -        -        0.5/1.6\n";
    }

    return card;
}

#endif
