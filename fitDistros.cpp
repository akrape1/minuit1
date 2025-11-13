#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>

using namespace std;

TH1F *hdata = nullptr;
TF1 *fparam = nullptr;

// -------------------------------------------------------------------------
// Model 1: sum of two Gaussians, amplitudes are peak heights (counts/bin)
double twoGauss(double *x, double *par) {
  double xx = x[0];
  double A1 = par[0];   // height of first Gaussian
  double mu1 = par[1];
  double sigma1 = par[2];
  double A2 = par[3];   // height of second Gaussian
  double mu2 = par[4];
  double sigma2 = par[5];
  double g1 = A1 * exp(-0.5 * pow((xx - mu1)/sigma1, 2));
  double g2 = A2 * exp(-0.5 * pow((xx - mu2)/sigma2, 2));
  return g1 + g2;   // counts per bin height
}

// -------------------------------------------------------------------------
// Model 2: Gumbel, amplitude is peak height (counts/bin)
double gumbel(double *x, double *par) {
  double xx = x[0];
  double A = par[0];    // peak height
  double mu = par[1];   // location
  double beta = par[2]; // scale
  double z = (xx - mu)/beta;
  return A * exp(-z - exp(-z)); //A * (1.0/beta) * exp(-z - exp(-z))
}

// -------------------------------------------------------------------------
// Chi2 calculator (models return counts/bin)
double calcChi2(TH1F* h, TF1* f) {
  double chi2 = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double x = h->GetBinCenter(i);
    double n = h->GetBinContent(i);
    double mu = f->Eval(x);
    double err = h->GetBinError(i);
    if (err <= 0) err = (n > 0) ? sqrt(n) : 1.0;
    chi2 += pow((n - mu)/err, 2);
  }
  return chi2;
}

// -------------------------------------------------------------------------
void fcn_chi2(int &npar, double *deriv, double &f, double par[], int flag) {
  for (int i = 0; i < npar; i++) fparam->SetParameter(i, par[i]);
  f = calcChi2(hdata, fparam);
}

// -------------------------------------------------------------------------
void doChi2Fit(TH1F* h, TF1* model, double *startPar, double *step,
               int npar, double &chi2, int &ndf, double *outpar, double *err) {

  TMinuit minuit(npar);
  minuit.SetFCN(fcn_chi2);

  for (int i = 0; i < npar; ++i)
    minuit.DefineParameter(i, Form("p%d", i), startPar[i], step[i], 0, 0);

  hdata = h;
  fparam = model;
  minuit.Migrad();

  double fmin, fedm, errdef;
  int npari, nparx, istat;
  minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  chi2 = fmin;
  ndf = h->GetNbinsX() - npar;

  for (int i = 0; i < npar; i++) minuit.GetParameter(i, outpar[i], err[i]);
  model->SetParameters(outpar);
}

// -------------------------------------------------------------------------
void fitDistros() {
  gStyle->SetOptStat(0);

  TFile *f = TFile::Open("distros.root");
  if (!f || f->IsZombie()) { cerr << "Cannot open distros.root\n"; return; }

  TH1F *h1 = (TH1F*)f->Get("dist1");
  if (!h1) { cerr << "Cannot find dist1\n"; return; }
  h1->Sumw2();

  // --- Two Gaussian model ---
  TF1 *fTwoG = new TF1("fTwoG", twoGauss,
                       h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(), 6);

  double ymax = h1->GetMaximum();
  double xmean = h1->GetMean();
  double xstd = h1->GetRMS();

  double startG[6] = {ymax, xmean - xstd, xstd/2.0,
                      ymax/2.0, xmean + xstd, xstd/2.0};
  double stepG[6]  = {0.05*ymax, 0.1*xstd, 0.05*xstd,
                      0.05*ymax, 0.1*xstd, 0.05*xstd};
  double chi2G, outG[6], errG[6]; int ndfG;

  doChi2Fit(h1, fTwoG, startG, stepG, 6, chi2G, ndfG, outG, errG);
  double pvalG = TMath::Prob(chi2G, ndfG);

  // --- Gumbel model ---
  TF1 *fGumb = new TF1("fGumb", gumbel,
                       h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(), 3);

  double startGb[3] = {ymax, xmean, xstd/2.0};
  double stepGb[3]  = {0.05*ymax, 0.1*xstd, 0.05*xstd};
  double chi2Gb, outGb[3], errGb[3]; int ndfGb;

  doChi2Fit(h1, fGumb, startGb, stepGb, 3, chi2Gb, ndfGb, outGb, errGb);
  double pvalGb = TMath::Prob(chi2Gb, ndfGb);

  // --- Print ---
  cout << "\n========== Results ==========\n";
  cout << "Two-Gaussian Fit:\n";
  cout << "  Chi2 = " << chi2G << "  ndf = " << ndfG
       << "  Chi2/ndf = " << chi2G/ndfG
       << "  p-value = " << pvalG << endl;

  cout << "\nGumbel Fit:\n";
  cout << "  Chi2 = " << chi2Gb << "  ndf = " << ndfGb
       << "  Chi2/ndf = " << chi2Gb/ndfGb
       << "  p-value = " << pvalGb << endl;

  cout << "\nPreferred model: "
       << ((pvalG > pvalGb) ? "Two-Gaussian" : "Gumbel")
       << " (higher p-value)\n==============================\n";

  // --- Plot ---
  TCanvas *c1 = new TCanvas("c1", "Fits Comparison", 1200, 600);
  c1->Divide(2,1);

  c1->cd(1);
  h1->Draw("E");
  fTwoG->SetLineColor(kRed);
  fTwoG->SetLineWidth(2);
  fTwoG->Draw("same");
  TLatex l1; l1.SetNDC();
  l1.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.2f", chi2G/ndfG));

  c1->cd(2);
  h1->Draw("E");
  fGumb->SetLineColor(kBlue);
  fGumb->SetLineWidth(2);
  fGumb->Draw("same");
  TLatex l2; l2.SetNDC();
  l2.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.2f", chi2Gb/ndfGb));

  c1->SaveAs("fitDistros_results.png");
}
