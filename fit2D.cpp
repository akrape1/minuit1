#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TMath.h>
#include <iostream>

TH2F *hdata    = nullptr;
TH2F *hbkg     = nullptr;
TH2F *hbkgPDF  = nullptr;  // normalized background shape
double hbkgNorm = 0.0;


double signalKernel(double x, double y, const double *p)
{
    double A   = p[0]; // total signal yield (used later with normalization)
    double mu1 = p[1];
    double mu2 = p[2];
    double s1  = p[3];
    double s2  = p[4];

    // return only the SHAPE (no A here)
    double kx = TMath::Exp(-0.5 * (x - mu1)*(x - mu1) / (s1*s1));
    double ky = TMath::Exp(-0.5 * (y - mu2)*(y - mu2) / (s2*s2));
    return kx * ky;
}

// ------------------------
// Full model = A * S(x,y) + B * Bkg(x,y)
// where S and Bkg are normalized shapes that integrate to 1 over the histogram
// ------------------------
double model(double x, double y, const double *p, double sigNorm)
{
    double A = p[0]; // total signal yield
    double B = p[5]; // total background yield

    // signal shape at (x,y)
    double kern = signalKernel(x, y, p);
    double sShape = (sigNorm > 0.0) ? (kern / sigNorm) : 0.0;

    // background shape from hbkgPDF
    int ix = hbkgPDF->GetXaxis()->FindBin(x);
    int iy = hbkgPDF->GetYaxis()->FindBin(y);
    double bShape = hbkgPDF->GetBinContent(ix, iy); // already normalized

    // expected counts in this bin
    return A * sShape + B * bShape;
}

// ------------------------
// chi2 function 
// ------------------------
void fcn(int &npar, double *, double &fval, double *p, int)
{
    // First, compute signal normalization over all bins
    double sigNorm = 0.0;
    for (int ix = 1; ix <= hdata->GetNbinsX(); ix++) {
        double x = hdata->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= hdata->GetNbinsY(); iy++) {
            double y = hdata->GetYaxis()->GetBinCenter(iy);
            sigNorm += signalKernel(x, y, p);
        }
    }

    double chi2 = 0.0;

    for (int ix = 1; ix <= hdata->GetNbinsX(); ix++) {
        double x = hdata->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= hdata->GetNbinsY(); iy++) {
            double y = hdata->GetYaxis()->GetBinCenter(iy);

            double d = hdata->GetBinContent(ix, iy);
            if (d < 0) continue;  // should not happen, but be safe

            double m = model(x, y, p, sigNorm);

            double err = (d > 0.0) ? TMath::Sqrt(d) : 1.0;
            double pull = (d - m) / err;
            chi2 += pull * pull;
        }
    }

    fval = chi2;
}

void fit2D()
{
    // --------------------------------------
    // Load histograms
    // --------------------------------------
    TFile *f = TFile::Open("fitInputs.root");

    hdata = (TH2F*)f->Get("hdata");
    hbkg  = (TH2F*)f->Get("hbkg");

    double dataInt = hdata->Integral();
    hbkgNorm = hbkg->Integral();

    std::cout << "Integral hdata = " << dataInt  << std::endl;
    std::cout << "Integral hbkg  = " << hbkgNorm << std::endl;

    // --------------------------------------
    // Build normalized background PDF (shape only)
    // --------------------------------------
    hbkgPDF = (TH2F*)hbkg->Clone("hbkgPDF");
    hbkgPDF->Scale(1.0 / hbkgNorm);

    // --------------------------------------
    // Find peak position in the data for initial mu1, mu2
    // --------------------------------------
    int binx, biny, binz;
    hdata->GetMaximumBin(binx, biny, binz);
    double mu1_start = hdata->GetXaxis()->GetBinCenter(binx);
    double mu2_start = hdata->GetYaxis()->GetBinCenter(biny);

    // --------------------------------------
    // MINUIT setup
    // Parameters: 0:A, 1:mu1, 2:mu2, 3:sigma1, 4:sigma2, 5:B
    // A and B are total yields; A+B ~ dataInt
    // --------------------------------------
    TMinuit min(6);
    min.SetFCN(fcn);

    double A_start = 0.2 * dataInt;  // guess 20% signal
    double B_start = 0.8 * dataInt;  // guess 80% background

    min.DefineParameter(0, "A",      A_start,    10.0, 0.0, dataInt*2.0);
    min.DefineParameter(1, "mu1",    mu1_start,  0.1,  mu1_start - 10.0, mu1_start + 10.0);
    min.DefineParameter(2, "mu2",    mu2_start,  0.1,  mu2_start - 10.0, mu2_start + 10.0);
    min.DefineParameter(3, "sigma1", 5.0,        0.1,  1.0, 20.0);
    min.DefineParameter(4, "sigma2", 5.0,        0.1,  1.0, 20.0);
    min.DefineParameter(5, "B",      B_start,    10.0, 0.0, dataInt*2.0);


    min.Migrad();

    // --------------------------------------
    // Retrieve parameters and print results
    // --------------------------------------
    const char* parName[6] = {"A", "mu1", "mu2", "sigma1", "sigma2", "B"};
    double p[6], e[6];

    std::cout << "\n=== Fit Parameter Results ===\n";
    for (int i = 0; i < 6; i++) {
        min.GetParameter(i, p[i], e[i]);
        std::cout << parName[i] << " = " << p[i] << " ± " << e[i] << std::endl;
    }

    // --------------------------------------
    // Compute chi2 and reduced chi2
    // --------------------------------------
    Double_t fmin, edm, errdef;
    Int_t npari, nparx, istat;
    min.mnstat(fmin, edm, errdef, npari, nparx, istat);

    double chi2 = fmin;

    int nbinsUsed = 0;
    for (int ix = 1; ix <= hdata->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= hdata->GetNbinsY(); iy++) {
            if (hdata->GetBinContent(ix, iy) >= 0) nbinsUsed++;
        }
    }

    int ndf = nbinsUsed - npari;
    double chi2red = (ndf > 0) ? (chi2 / ndf) : 0.0;

    std::cout << "\nChi2 = " << chi2
              << "   NDF = " << ndf
              << "   Reduced Chi2 = " << chi2red << "\n\n";

    // --------------------------------------
    // Recompute signal normalization with final params for plotting
    // --------------------------------------
    double sigNorm = 0.0;
    for (int ix = 1; ix <= hdata->GetNbinsX(); ix++) {
        double x = hdata->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= hdata->GetNbinsY(); iy++) {
            double y = hdata->GetYaxis()->GetBinCenter(iy);
            sigNorm += signalKernel(x, y, p);
        }
    }

    // --------------------------------------
    // Build output histograms
    // --------------------------------------
    TH2F *hfit    = (TH2F*)hdata->Clone("hfit");    hfit->Reset();
    TH2F *hres    = (TH2F*)hdata->Clone("hres");    hres->Reset();
    TH2F *hsub    = (TH2F*)hdata->Clone("hsub");    hsub->Reset();
    TH2F *hsig    = (TH2F*)hdata->Clone("hsig");    hsig->Reset();
    TH2F *hbkgfit = (TH2F*)hdata->Clone("hbkgfit"); hbkgfit->Reset();

    for (int ix = 1; ix <= hdata->GetNbinsX(); ix++) {
        double x = hdata->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= hdata->GetNbinsY(); iy++) {
            double y = hdata->GetYaxis()->GetBinCenter(iy);

            double kern = signalKernel(x, y, p);
            double sShape = (sigNorm > 0.0) ? (kern / sigNorm) : 0.0;
            double sigVal = p[0] * sShape;

            double bShape = hbkgPDF->GetBinContent(ix, iy);
            double bkgVal = p[5] * bShape;

            double m = sigVal + bkgVal;
            double d = hdata->GetBinContent(ix, iy);

            hfit->SetBinContent(ix,  iy, m);
            hres->SetBinContent(ix,  iy, d - m);
            hsub->SetBinContent(ix,  iy, d - bkgVal);
            hsig->SetBinContent(ix,  iy, sigVal);
            hbkgfit->SetBinContent(ix,iy, bkgVal);
        }
    }

    // --------------------------------------
    // Make the 2×2 canvas
    // --------------------------------------
    TCanvas *c = new TCanvas("c", "2D Fit Results", 1200, 1000);
    c->Divide(2,2);

    c->cd(1);
    hdata->SetTitle("Data");
    hdata->Draw("lego2");

    c->cd(2);
    hfit->SetTitle("Fit Result (Signal + Background)");
    hfit->Draw("lego2");

    c->cd(3);
    hres->SetTitle("Residuals (Data - Fit)");
    hres->Draw("lego2");

    c->cd(4);
    hsub->SetTitle("Data - Best-Fit Background");
    hsub->Draw("lego2");

    c->SaveAs("fit2D_results.png");
}
