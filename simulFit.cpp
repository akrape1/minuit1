#include <TH1.h>
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMinuit.h>

TH1F *h1 = nullptr;
TH1F *h2 = nullptr;

// --- Global parameters ---
// p[0] = mean
// p[1] = sigma
// p[2] = A1 (signal norm exp1)
// p[3] = B1 (background norm exp1)
// p[4] = A2 (signal norm exp2)
// p[5] = B2 (background norm exp2)
// p[6] = n  (power exponent exp2)

double chi2_global(const double *p)
{
    double chi2 = 0.0;

    // ---------- Experiment 1 model ----------
    for (int i = 1; i <= h1->GetNbinsX(); i++)
    {
        double x = h1->GetBinCenter(i);
        double y = h1->GetBinContent(i);
        double err = sqrt(y);
        if (err < 1) err = 1;

        double sig = p[2] * exp(-0.5 * pow((x - p[0]) / p[1], 2));
        double bkg = p[3] * exp(-x / 20.0);
        double model = sig + bkg;

        chi2 += pow((y - model) / err, 2);
    }

    // ---------- Experiment 2 model ----------
    for (int i = 1; i <= h2->GetNbinsX(); i++)
    {
        double x = h2->GetBinCenter(i);
        double y = h2->GetBinContent(i);
        double err = sqrt(y);
        if (err < 1) err = 1;

        double sig = p[4] * exp(-0.5 * pow((x - p[0]) / p[1], 2));
        double bkg = p[5] * pow(x, p[6]);
        double model = sig + bkg;

        chi2 += pow((y - model) / err, 2);
    }

    return chi2;
}

// MINUIT objective
void fcn(int &npar, double *, double &f, double *par, int /*flag*/)
{
    f = chi2_global(par);
}

void simulFit()
{
    // --- Load histograms ---
    TFile *f = TFile::Open("experiments.root");
    h1 = (TH1F*)f->Get("hexp1");
    h2 = (TH1F*)f->Get("hexp2");

    // --- Minuit ---
    TMinuit min(7);
    min.SetFCN(fcn);

    // initial parameters
    min.DefineParameter(0, "mean",    75.0, 0.1, 60, 90);
    min.DefineParameter(1, "sigma",   4.5,  0.1, 1, 10);
    min.DefineParameter(2, "A1",      40,   1, 0, 500);
    min.DefineParameter(3, "B1",      50,   1, 0, 5000);
    min.DefineParameter(4, "A2",      40,   1, 0, 500);
    min.DefineParameter(5, "B2",   5e4,    50, 0, 1e7);
    min.DefineParameter(6, "n",       -2.2, 0.05, -10, 0);

    min.Migrad();

        // --- parameters ---
    double mean, emean, sigma, esigma;
    min.GetParameter(0, mean,  emean);
    min.GetParameter(1, sigma, esigma);

    // --- chi2 and ndof ---
    Double_t fmin, edm, errdef; 
    Int_t npari, nparx, istat;
    min.mnstat(fmin, edm, errdef, npari, nparx, istat);

    double chi2 = fmin;
    int npoints = h1->GetNbinsX() + h2->GetNbinsX();
    int ndf    = npoints - npari;
    double prob = TMath::Prob(chi2, ndf);

    // --- print summary ---
    std::cout << "Mean  = " << mean  << " ± " << emean  << std::endl;
    std::cout << "Sigma = " << sigma << " ± " << esigma << std::endl;
    std::cout << "chi2/ndf = " << chi2 << " / " << ndf
            << " = " << chi2/ndf << std::endl;
    std::cout << "Chi2 probability = " << prob << std::endl;

    double p[7], e[7];
    for (int i = 0; i < 7; i++)
        min.GetParameter(i, p[i], e[i]);

    TF1 *f1 = new TF1("fit_exp1",
                      [&](double *x, double *) {
                          double xx = x[0];
                          return p[2]*exp(-0.5*pow((xx-p[0])/p[1],2))
                                 + p[3]*exp(-xx/20.0);
                      },
                      h1->GetXaxis()->GetXmin(),
                      h1->GetXaxis()->GetXmax(), 0);

    TF1 *f2 = new TF1("fit_exp2",
                      [&](double *x, double *) {
                          double xx = x[0];
                          return p[4]*exp(-0.5*pow((xx-p[0])/p[1],2))
                                 + p[5]*pow(xx, p[6]);
                      },
                      h2->GetXaxis()->GetXmin(),
                      h2->GetXaxis()->GetXmax(), 0);

    f1->SetLineColor(kRed); f1->SetLineWidth(2);
    f2->SetLineColor(kRed); f2->SetLineWidth(2);

    TCanvas *c = new TCanvas("c","Simultaneous Fit",1200,500);
    c->Divide(2,1);

    c->cd(1);
    h1->Draw();
    f1->Draw("same");

    c->cd(2);
    h2->Draw();
    f2->Draw("same");

    c->SaveAs("simulFit_results.png");
}
