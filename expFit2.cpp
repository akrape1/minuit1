// this is generated to be like the example one but executable thru ROOT. chi2 is done here

#include <iostream>
#include <cmath>
#include <TMinuit.h>
#include <TRandom2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TAxis.h>

using namespace std;

// ---------- GLOBALS ----------
TH1F *hdata = nullptr;
TF1 *fparam = nullptr;

//-------------------------------------------------------------------------
// Exponential PDF
double expPdf(double* xPtr, double par[]) {
  double x = *xPtr;
  double A = par[0];
  double lam = par[1];
  return A * exp(x / lam);
}

//-------------------------------------------------------------------------
// Negative log-likelihood
double calcNLL(TH1F* h, TF1* f) {
  double nll = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double x = h->GetBinCenter(i);
    int n = (int)(h->GetBinContent(i));
    double mu = f->Eval(x);
    if (mu < 1e-10) mu = 1e-10;
    nll -= n * TMath::Log(mu) - mu - TMath::LnGamma(n + 1);
  }
  return 2 * nll;
}

//-------------------------------------------------------------------------
// Chi2 calculation
double calcChi2(TH1F* h, TF1* f) {
  double chi2 = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double x = h->GetBinCenter(i);
    double n = h->GetBinContent(i);
    double mu = f->Eval(x);
    double err = h->GetBinError(i);
    if (err <= 0) err = (n > 0) ? sqrt(n) : 1.0;
    chi2 += pow((n - mu) / err, 2);
  }
  return chi2;
}

//-------------------------------------------------------------------------
// FCNs for Minuit
void fcn_nll(int& npar, double* deriv, double& f, double par[], int flag) {
  for (int i = 0; i < npar; i++) fparam->SetParameter(i, par[i]);
  f = calcNLL(hdata, fparam);
}

void fcn_chi2(int& npar, double* deriv, double& f, double par[], int flag) {
  for (int i = 0; i < npar; i++) fparam->SetParameter(i, par[i]);
  f = calcChi2(hdata, fparam);
}


void expFit2(bool useChi2 = false) {

  cout << "\n==========================\n";
  cout << (useChi2 ? "Performing Chi2 fit..." : "Performing NLL fit...") << endl;

  // Style
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "hxy");
  gStyle->SetLabelFont(42, "xyz");

  // Generate data
  TRandom2 r(0);
  double xmin = 0.0, xmax = 4.0, lambda = 3.14159;
  int nentries = 1000;
  TH1F* hexp = new TH1F("hexp", "Exponential distribution;x;# events", 50, xmin, xmax);
  for (int i = 0; i < nentries; i++) hexp->Fill(r.Exp(lambda));
  hexp->Sumw2();

  // Prepare model
  const int npar = 2;
  TF1* myfunc = new TF1("myfunc", expPdf, xmin, xmax, npar);

  TMinuit minuit(npar);
  minuit.SetFCN(useChi2 ? fcn_chi2 : fcn_nll);

  double par[npar] = {hexp->GetMaximum(), -2};
  double stepSize[npar] = {fabs(par[0]*0.1), fabs(par[1]*0.1)};
  TString parName[npar] = {"A", "lambda"};

  for (int i = 0; i < npar; i++)
    minuit.DefineParameter(i, parName[i].Data(), par[i], stepSize[i], 0, 0);

  hdata = hexp;
  fparam = myfunc;

  // Fit
  minuit.Migrad();

  // Extract results
  double outpar[npar], err[npar];
  for (int i = 0; i < npar; i++) minuit.GetParameter(i, outpar[i], err[i]);
  myfunc->SetParameters(outpar);

  // Draw
  TCanvas* c1 = new TCanvas("c1", "Fit Result", 800, 600);
  hexp->Draw("E");
  myfunc->SetLineWidth(2);
  myfunc->Draw("same");
  c1->SaveAs("./test.png");

  // Print summary
  double fmin, fedm, errdef;
  int npari, nparx, istat;
  minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);

  cout << "Fit minimum = " << fmin << endl;
  cout << "Fit status = " << istat << endl;
  for (int i = 0; i < npar; ++i)
    cout << parName[i] << " = " << outpar[i] << " +- " << err[i] << endl;

  if (useChi2) {
    int ndf = hexp->GetNbinsX() - npar;
    cout << "Chi2/ndf = " << fmin / ndf << endl;
  }

  cout << "==========================\n";
}
