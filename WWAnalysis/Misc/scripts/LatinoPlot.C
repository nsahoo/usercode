#include<vector>

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFrame.h"
#endif

enum samp { iHWW, iWW, iZJets, iTop, iWZ, iWJets, nSamples };

float xPos[nSamples+1] = {0.22,0.22,0.22,0.39,0.39,0.39,0.39}; 
float yOff[nSamples+1] = {0,1,2,0,1,2,3};

const Float_t _tsize   = 0.03;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis*  axis,
        TString coordinate,
        TString title)
{
    axis->SetLabelFont  (   42);
    axis->SetLabelOffset(0.015);
    axis->SetLabelSize  (0.050);
    axis->SetNdivisions (  505);
    axis->SetTitleFont  (   42);
    axis->SetTitleOffset(  1.5);
    axis->SetTitleSize  (0.050);

    if (coordinate == "y") axis->SetTitleOffset(1.6);

    axis->SetTitle(title);
}


//------------------------------------------------------------------------------
// THStackAxisFonts
//------------------------------------------------------------------------------
void THStackAxisFonts(THStack* h,
        TString  coordinate,
        TString  title)
{
    TAxis* axis = NULL;

    if (coordinate.Contains("x")) axis = h->GetHistogram()->GetXaxis();
    if (coordinate.Contains("y")) axis = h->GetHistogram()->GetYaxis();

    AxisFonts(axis, coordinate, title);
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend(Float_t x1,
        Float_t y1,
        TH1F*   hist,
        TString label,
        TString option)
{
    TLegend* legend = new TLegend(x1,
            y1,
            x1 + _xoffset,
            y1 + _yoffset);

    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (_tsize);

    legend->AddEntry(hist, label.Data(), option.Data());

    legend->Draw();
}


class LatinoPlot {

    public: 
        LatinoPlot() { _hist.resize(nSamples,0); _data = 0; _breakdown = false;}
        void setMCHist   (const samp &s, TH1F * h)  { _hist[s]       = h;  } 
        void setDataHist (TH1F * h)                 { _data          = h;  } 
        void setHWWHist  (TH1F * h)                 { setMCHist(iHWW  ,h); } 
        void setWWHist   (TH1F * h)                 { setMCHist(iWW   ,h); } 
        void setZJetsHist(TH1F * h)                 { setMCHist(iZJets,h); } 
        void setTopHist  (TH1F * h)                 { setMCHist(iTop  ,h); } 
        void setWZHist   (TH1F * h)                 { setMCHist(iWZ   ,h); } 
        void setWJetsHist(TH1F * h)                 { setMCHist(iWJets,h); } 

        void Draw(const int &rebin=1) {

            Color_t _sampleColor[nSamples];
            _sampleColor[iHWW  ] = kRed+1;
            _sampleColor[iWW   ] = kAzure-9;
            _sampleColor[iZJets] = kGreen+2;
            _sampleColor[iTop  ] = kYellow;
            _sampleColor[iWZ   ] = kAzure-2;
            _sampleColor[iWJets] = kBlack;

            //setUpStyle();
            if(!gPad) new TCanvas();

            THStack* hstack = new THStack();
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;

                _hist[i]->Rebin(rebin);
                _hist[i]->SetLineColor(_sampleColor[i]);

                // signal gets overlaid
                if (i == iHWW) continue;

                _hist[i]->SetFillColor(_sampleColor[i]);
                _hist[i]->SetFillStyle(1001);

                hstack->Add(_hist[i]);
            }

            if(_hist[iHWW]) _hist[iHWW]->SetLineWidth(3);
            if(_data) _data->Rebin(rebin);
            if(_data) _data->SetLineColor  (kBlack);
            if(_data) _data->SetMarkerStyle(kFullCircle);

            hstack->Draw("hist");
            if(_hist[iHWW]) _hist[iHWW]->Draw("hist,same");
            if(_data) _data->Draw("ep,same");

            hstack->SetTitle("CMS preliminary");

            if( gPad->GetLogy() )
                hstack->SetMaximum(1000 * hstack->GetMaximum());
            else
                hstack->SetMaximum(1.55 * hstack->GetMaximum());

            if(_breakdown) {
                THStackAxisFonts(hstack, "y", "entries");
                hstack->GetHistogram()->LabelsOption("v");
            } else {
                THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                THStackAxisFonts(hstack, "y", TString::Format("entries / %.0f %s", _hist[iHWW]->GetBinWidth(0),_units.Data()));
            }

            // total mess to get it nice, should be redone
            size_t j=0;
            if(_data        ) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _data,         " data",  "lp"); j++; }
            if(_hist[iHWW  ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHWW  ], " HWW",    "l"); j++; }
            if(_hist[iWW   ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWW   ], " WW",     "f"); j++; }
            if(_hist[iZJets]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZJets], " Z+jets", "f"); j++; }
            if(_hist[iTop  ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iTop  ], " top",    "f"); j++; }
            if(_hist[iWZ   ]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWZ   ], " WW/WZ",  "f"); j++; }
            if(_hist[iWJets]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWJets], " W+jets", "f"); j++; }


            TLatex* luminosity = new TLatex(0.9, 0.815, TString::Format("L = %g pb^{-1}",_lumi));
            luminosity->SetNDC();
            luminosity->SetTextAlign(32);
            luminosity->SetTextFont(42);
            luminosity->SetTextSize(_tsize);
            luminosity->Draw("same");

        }
        void setLumi(const float &l) { _lumi = l; }
        void setLabel(const TString &s) { _xLabel = s; }
        void setUnits(const TString &s) { _units = s; }
        void setBreakdown(const bool &b = true) { _breakdown = b; }

    private: 
        std::vector<TH1F*> _hist;
        TH1F* _data;

        //MWL
        float _lumi;
        TString _xLabel;
        TString _units;
        bool _breakdown;


};


