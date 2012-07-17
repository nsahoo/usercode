#ifndef CARDTEMPLETE_H
#define CARDTEMPLETE_H

std::string createCardTemplate(bool do7TeV, int channel, bool do1D, std::string workspacefilename) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
    if (channel == 0) {
        card += std::string("shapes * * ") + workspacefilename + " w_4mu:$PROCESS\n";
    }
    if (channel == 1) {
        card += std::string("shapes * * ") + workspacefilename + " w_4e:$PROCESS\n";
    }
    if (channel == 2) {
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
    if (do7TeV) {        
        card += "lumi_7TeV                 lnN        1.022   1.022   1.022    1.022    -\n";
    }
    else {        
        card += "lumi_8TeV                 lnN        1.05    1.05    1.05     1.05     -\n";
    }
        card += "pdf_gg                    lnN        GGH_PDF -       -        GGZZ_PDF -\n";
        card += "pdf_qqbar                 lnN         -      VBF_PDF QQZZ_PDF -        -\n";
        card += "pdf_hzz4l_accept          lnN        1.02    1.02    -        -        -\n";
        card += "QCDscale_ggH              lnN        GGH_QCD -       -        -        -\n";
        card += "QCDscale_qqH              lnN        -       VBF_PDF -        -        -\n";
        card += "QCDscale_ggVV             lnN        -       -       -        GGZZ_QCD -\n";
        card += "QCDscale_VV               lnN        -       -       QQZZ_QCD -        -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02    1.02    -        -        -\n";
    if (channel == 0) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015    -\n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -       -        -        0.5/1.6\n";
        if (do7TeV) {
        card += "sig_4mu_mean_err_7TeV     param      MEAN_0  MEAN_1                     \n";
        card += "sig_4mu_sigma_err_7TeV    param      SIGMA_0 SIGMA_1                    \n";
        }
        else {
        card += "sig_4mu_mean_err_8TeV     param      MEAN_0  MEAN_1                     \n";
        card += "sig_4mu_sigma_err_8TeV    param      SIGMA_0 SIGMA_1                    \n";
        }
    }
    else if (channel == 1) {
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -\n";
        card += "CMS_hzz4e_Zjets           lnN        -       -       -        -        0.5/1.6\n";
        if (do7TeV) {
        card += "sig_4e_mean_err_7TeV      param      MEAN_0  MEAN_1                     \n";
        card += "sig_4e_sigma_err_7TeV     param      SIGMA_0 SIGMA_1                    \n";
        }
        else {
        card += "sig_4e_mean_err_8TeV      param      MEAN_0  MEAN_1                     \n";
        card += "sig_4e_sigma_err_8TeV     param      SIGMA_0 SIGMA_1                    \n";
        }
    }
    else if (channel == 2) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -\n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -        -        0.5/1.6\n";
        if (do7TeV) {
        card += "sig_2e2mu_mean_err_7TeV   param      MEAN_0  MEAN_1                     \n";
        card += "sig_2e2mu_sigma_err_7TeV  param      SIGMA_0 SIGMA_1                    \n";
        }
        else {
        card += "sig_2e2mu_mean_err_8TeV   param      MEAN_0  MEAN_1                     \n";
        card += "sig_2e2mu_sigma_err_8TeV  param      SIGMA_0 SIGMA_1                    \n";
        }
    }
    if (!do1D) {
        card += "CMS_zz4l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    return card;
}

#endif
