#ifndef CARDTEMPLETE_H
#define CARDTEMPLETE_H

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
        card += "bin     BIN      BIN        BIN       BIN       BIN\n";
        card += "process sig_ggH sig_VBF  bkg_qqzz bkg_ggzz bkg_zjets\n";
        card += "process -2     -1        1        2        3\n";
        card += "rate    SIG_GGH SIG_VBF  BKG_QQZZ BKG_GGZZ BKG_ZJETS\n";
        card += "------------\n";
        card += "lumi                      lnN        1.045  1.045  1.045   1.045  -\n";
        card += "pdf_gg                    lnN        1.0755 -      -       1.0708 -\n";
        card += "pdf_qqbar                 lnN         -     1.023  1.0341  -      -\n";
        card += "pdf_H4l_accept            lnN        1.02   1.02   -       -      -\n";
        card += "QCDscale_ggH              lnN        1.0975 -      -       -      -\n";
        card += "QCDscale_qqH              lnN        -      1.003  -       -      -\n";
        card += "QCDscale_VH               lnN        -      -      -       -      -\n";
        card += "theoryUncXS_HighMH        lnN        1      1      -       -      -\n";
        card += "QCDscale_ggVV             lnN        -      -      -       1.2431 -\n";
        card += "QCDscale_VV               lnN        -      -      1.0284  -      -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02   1.02   -       -      -\n";
    if (channel == 0 || channel == 2) {
        card +=  "CMS_eff_m                lnN        1.0224 1.0224 1.0224  1.0224 -\n";
        card +=  "CMS_scale_m              lnN        1.01   1.01   1.01    1.01   -\n";
    }
    if (channel == 1 || channel == 2) {
        card += "CMS_eff_e                 lnN        1.0447 1.0447 1.0447  1.0447 -\n";
        card += "CMS_scale_e               lnN        1.02   1.02   1.02    1.02   -\n";
    }
        card += "CMS_hzz4mu_Zjets_EXTRAP   lnN        -      -      -       -      2.00\n";

    return card;
}

#endif
