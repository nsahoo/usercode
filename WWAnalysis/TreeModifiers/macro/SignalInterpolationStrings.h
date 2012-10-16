#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

std::string getSignalCBMeanString(float m, int ch, bool do7TeV) {
    stringstream fss;
    fss << "( ";  

    if (!doFFT) fss << "@0 + ";

    if (m<400.) {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    else {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    fss << " ) + @0*@1";
    return fss.str();
}

std::string getSignalCBSigmaString(float m, int ch, do7TeV) {
    stringstream fss;
    fss << "( ";  

    if (m<400.) {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    else {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    fss << " ) * (1+@1)";
    return fss.str();
}

std::string getSignalCBAlphaLString(float m, int ch, bool do7TeV) {
    stringstream fss;
   
    if (m<400.) {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    else {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }

    return fss.str();
}

std::string getSignalCBNLString(float m, int ch, bool do7TeV) {
    stringstream fss;

    if (m<400.) {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    else {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }


    return fss.str();
}

std::string getSignalCBAlphaRString(float m, int ch, bool do7TeV) {
    stringstream fss;

    if (m<400.) {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    else {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }


    return fss.str();
}

std::string getSignalCBNRString(float m, int ch, bool do7TeV) {
    stringstream fss;

    if (m<400.) {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }
    else {
        if (do7TeV) {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
        else {
            if (ch == 0) fss << "";
            if (ch == 1) fss << "";
            if (ch == 2) fss << "";
        }
    }


    return fss.str();
}


std::string getSignalBWGammaString(float m, int ch, bool do7TeV) {
    stringstream fss;

    if (do7TeV) {
        if (ch == 0) fss << "";
        if (ch == 1) fss << "";
        if (ch == 2) fss << "";
    }
    else {
        if (ch == 0) fss << "";
        if (ch == 1) fss << "";
        if (ch == 2) fss << "";
    }


    return fss.str();
}


std::string getYieldEfficiencyString(float m, int ch, bool do7TeV) {
    stringstream fss;

    if (do7TeV) {
        if      (ch == 0) return "(-2.52639) + (0.0553644)*@0 + (-0.000405676)*@0*@0 + (1.61242e-06)*@0*@0*@0 + (-3.78806e-09)*@0*@0*@0*@0 + (5.41746e-12)*@0*@0*@0*@0*@0 + (-4.63766e-15)*@0*@0*@0*@0*@0*@0 + (2.18648e-18)*@0*@0*@0*@0*@0*@0*@0 + (-4.37074e-22)*@0*@0*@0*@0*@0*@0*@0*@0";
        else if (ch == 1) return "(-2.59319) + (0.0514931)*@0 + (-0.000369965)*@0*@0 + (1.46651e-06)*@0*@0*@0 + (-3.49839e-09)*@0*@0*@0*@0 + (5.1664e-12)*@0*@0*@0*@0*@0 + (-4.63104e-15)*@0*@0*@0*@0*@0*@0 + (2.30959e-18)*@0*@0*@0*@0*@0*@0*@0 + (-4.91469e-22)*@0*@0*@0*@0*@0*@0*@0*@0";
        else              return "(-2.76074) + (0.0595281)*@0 + (-0.000455363)*@0*@0 + (1.90513e-06)*@0*@0*@0 + (-4.7343e-09)*@0*@0*@0*@0 + (7.17618e-12)*@0*@0*@0*@0*@0 + (-6.50896e-15)*@0*@0*@0*@0*@0*@0 + (3.24419e-18)*@0*@0*@0*@0*@0*@0*@0 + (-6.83026e-22)*@0*@0*@0*@0*@0*@0*@0*@0";
    }
    else {
        if      (ch == 0) return "(-3.28633) + (0.0682262)*@0 + (-0.000496505)*@0*@0 + (1.9491e-06)*@0*@0*@0 + (-4.50967e-09)*@0*@0*@0*@0 + (6.33717e-12)*@0*@0*@0*@0*@0 + (-5.31633e-15)*@0*@0*@0*@0*@0*@0 + (2.44814e-18)*@0*@0*@0*@0*@0*@0*@0 + (-4.76086e-22)*@0*@0*@0*@0*@0*@0*@0*@0";
        else if (ch == 1) return "(-2.59332) + (0.0490742)*@0 + (-0.000332895)*@0*@0 + (1.22995e-06)*@0*@0*@0 + (-2.70031e-09)*@0*@0*@0*@0 + (3.63235e-12)*@0*@0*@0*@0*@0 + (-2.94501e-15)*@0*@0*@0*@0*@0*@0 + (1.32399e-18)*@0*@0*@0*@0*@0*@0*@0 + (-2.53976e-22)*@0*@0*@0*@0*@0*@0*@0*@0";
        else              return "(-3.02899) + (0.0623266)*@0 + (-0.000461463)*@0*@0 + (1.86983e-06)*@0*@0*@0 + (-4.51599e-09)*@0*@0*@0*@0 + (6.68363e-12)*@0*@0*@0*@0*@0 + (-5.94745e-15)*@0*@0*@0*@0*@0*@0 + (2.92147e-18)*@0*@0*@0*@0*@0*@0*@0 + (-6.08659e-22)*@0*@0*@0*@0*@0*@0*@0*@0";
    }

    return fss.str();
}

