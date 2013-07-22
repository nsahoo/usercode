#include "TCut.h"
#include "TH2D.h"
#include "TMarker.h"

TGraph* bestFit(TTree *t, TString x, TString y, TCut cut) {
    int nfind = t->Draw(y+":"+x, cut + "deltaNLL == 0");
    if (nfind == 0) {
        TGraph *gr0 = new TGraph(1);
        gr0->SetPoint(0,-999,-999);
        gr0->SetMarkerStyle(34); gr0->SetMarkerSize(2.0);
        return gr0;
    } else {
        TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
        gr0->SetMarkerStyle(34); gr0->SetMarkerSize(2.0);
        if (gr0->GetN() > 1) gr0->Set(1);
        return gr0;
    }
}

TH2 *treeToHist2D(TTree *t, TString x, TString y, TString name, TCut cut, double xmin, double xmax, double ymin, double ymax, int xbins, int ybins) {
    t->Draw(Form("2*deltaNLL:%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)", y.Data(), x.Data(), name.Data(), xbins, xmin, xmax, ybins, ymin, ymax), cut + "deltaNLL != 0", "PROF");
    TH2 *prof = (TH2*) gROOT->FindObject(name+"_prof");
    TH2D *h2d = new TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
    for (int ix = 1; ix <= xbins; ++ix) {
        for (int iy = 1; iy <= ybins; ++iy) {
             double z = prof->GetBinContent(ix,iy);
             if (z != z) z = (name.Contains("bayes") ? 0 : 999); // protect agains NANs
             h2d->SetBinContent(ix, iy, z);
        }
    }
    h2d->SetDirectory(0);
    return h2d;
}

TH2D* frameTH2D(TH2D *in, double threshold){
        // NEW LOGIC:
        //   - pretend that the center of the last bin is on the border if the frame
        //   - add one tiny frame with huge values
        double frameValue = 1000;
        if (TString(in->GetName()).Contains("bayes")) frameValue = 0.0;

	Double_t xw = in->GetXaxis()->GetBinWidth(1);
	Double_t yw = in->GetYaxis()->GetBinWidth(1);

	Int_t nx = in->GetNbinsX();
	Int_t ny = in->GetNbinsY();

	Double_t x0 = in->GetXaxis()->GetXmin();
	Double_t x1 = in->GetXaxis()->GetXmax();

	Double_t y0 = in->GetYaxis()->GetXmin();
	Double_t y1 = in->GetYaxis()->GetXmax();
        Double_t xbins[999], ybins[999]; 
        double eps = 0.1;

        xbins[0] = x0 - eps*xw - xw; xbins[1] = x0 + eps*xw - xw;
        for (int ix = 2; ix <= nx; ++ix) xbins[ix] = x0 + (ix-1)*xw;
        xbins[nx+1] = x1 - eps*xw + 0.5*xw; xbins[nx+2] = x1 + eps*xw + xw;

        ybins[0] = y0 - eps*yw - yw; ybins[1] = y0 + eps*yw - yw;
        for (int iy = 2; iy <= ny; ++iy) ybins[iy] = y0 + (iy-1)*yw;
        ybins[ny+1] = y1 - eps*yw + yw; ybins[ny+2] = y1 + eps*yw + yw;
        
	TH2D *framed = new TH2D(
			Form("%s framed",in->GetName()),
			Form("%s framed",in->GetTitle()),
			nx + 2, xbins,
			ny + 2, ybins 
			);

	//Copy over the contents
	for(int ix = 1; ix <= nx ; ix++){
		for(int iy = 1; iy <= ny ; iy++){
			framed->SetBinContent(1+ix, 1+iy, in->GetBinContent(ix,iy));
		}
	}
	//Frame with huge values
	nx = framed->GetNbinsX();
	ny = framed->GetNbinsY();
	for(int ix = 1; ix <= nx ; ix++){
		framed->SetBinContent(ix,  1, frameValue);
		framed->SetBinContent(ix, ny, frameValue);
	}
	for(int iy = 2; iy <= ny-1 ; iy++){
		framed->SetBinContent( 1, iy, frameValue);
		framed->SetBinContent(nx, iy, frameValue);
	}

	return framed;
}

TList* contourFromTH2(TH2 *h2in, double threshold, int minPoints=20) {
    std::cout << "Getting contour at threshold " << threshold << " from " << h2in->GetName() << std::endl;
    //http://root.cern.ch/root/html/tutorials/hist/ContourList.C.html
    Double_t contours[1];
    contours[0] = threshold;
    if (h2in->GetNbinsX() * h2in->GetNbinsY() > 10000) minPoints = 50;

    TH2D *h2 = frameTH2D((TH2D*)h2in,threshold);

    h2->SetContour(1, contours);

    // Draw contours as filled regions, and Save points
    h2->Draw("CONT Z LIST");
    gPad->Update(); // Needed to force the plotting and retrieve the contours in TGraphs


    // Get Contours
    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;

    if (conts == NULL || conts->GetSize() == 0){
        printf("*** No Contours Were Extracted!\n");
        return 0;
    }

    TList *ret = new TList();
    for(int i = 0; i < conts->GetSize(); i++){
        contLevel = (TList*)conts->At(i);
        printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
        for (int j = 0, n = contLevel->GetSize(); j < n; ++j) {
            TGraph *gr1 = (TGraph*) contLevel->At(j);
            if (gr1->GetN() > minPoints) ret->Add(gr1->Clone());
            //break;
        }
    }
    return ret;
}

void styleMultiGraph(TList *tmg, int lineColor, int lineWidth, int lineStyle) {
    for (int i = 0; i < tmg->GetSize(); ++i) {
        TGraph *g = (TGraph*) tmg->At(i);
        g->SetLineColor(lineColor); g->SetLineWidth(lineWidth); g->SetLineStyle(lineStyle);
    }
}
void styleMultiGraphMarker(TList *tmg, int markerColor, int markerSize, int markerStyle) {
    for (int i = 0; i < tmg->GetSize(); ++i) {
        TGraph *g = (TGraph*) tmg->At(i);
        g->SetMarkerColor(markerColor); g->SetMarkerSize(markerSize); g->SetMarkerStyle(markerStyle);
    }
}

void contour(TString xvar, float xmin, float xmax, TString yvar, float ymin, float ymax) {
    TTree *tree = (TTree*) gFile->Get("limit") ;
    TH2 *hist2d = treeToHist2D(tree, xvar, yvar, "h2d", "", xmin, xmax, ymin, ymax, 300, 300);
    TGraph *fit = bestFit(tree, xvar, yvar, "");
    TList *c68 = contourFromTH2(hist2d, 2.30);
    TList *c95 = contourFromTH2(hist2d, 5.99);
    styleMultiGraph(c68, /*color=*/1, /*width=*/3, /*style=*/1);
    styleMultiGraph(c95, /*color=*/1, /*width=*/3, /*style=*/9);
    hist2d->Draw("AXIS");
    hist2d->Draw("COLZ 9 SAME");
    c68->Draw("L SAME");
    c95->Draw("L SAME");
    fit->Draw("P SAME");
    double smx = 1.0, smy = 1.0; TMarker m;
    m.SetMarkerSize(3.0); m.SetMarkerColor(97); m.SetMarkerStyle(33); 
    m.DrawMarker(smx,smy);
    m.SetMarkerSize(1.8); m.SetMarkerColor(89); m.SetMarkerStyle(33); 
    m.DrawMarker(smx,smy);
}
void contours() {}
