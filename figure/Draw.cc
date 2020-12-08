#include "TFile.h"
#include "TTree.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <TH1D.h>
#include <TStyle.h>
#include "SetDrawNcuOpt.h"
#include "THStack.h"

void Time();

void Draw()
{
    //Time();
    //Velocity();
    Momentum();
}

void Time()
{
    auto c1 = new TCanvas("c", "BPRE", 10, 10, 800, 600);
    double mean1, rms1, mean2, rms2;
    TFile *file1 = new TFile("../output/tev5mm_pythia6_zprime5tev_ww_0199_hepsim.slcio.root");
    TFile *file2 = new TFile("../output/tev5mm_pythia6_zprime5tev_qq_0199_hepsim.slcio.root");

    //TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Leading"));
    //TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Leading"));
    TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Trailing"));
    TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Trailing"));

    mean1 = sig->GetMean();
    rms1 = sig->GetRMS();

    mean2 = bkg->GetMean();
    rms2 = bkg->GetRMS();

    SetDrawNcuOpt(sig, kRed, "Time of flight - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    SetDrawNcuOpt(bkg, kBlue, "Time of flight - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    SetRootStyle();

    sig->Scale(1.0 / sig->Integral());
    bkg->Scale(1.0 / bkg->Integral());

    sig->GetXaxis()->SetRangeUser(6.5, 24.5);
    sig->GetYaxis()->SetRangeUser(0, 0.3);
    bkg->GetXaxis()->SetRangeUser(6.5, 24.5);
    bkg->GetYaxis()->SetRangeUser(0, 0.3);

    bkg->DrawNormalized("hist");
    sig->DrawNormalized("same&&hist");

    TLegend *leg2 = new TLegend(0.5, 0.45, 0.89, 0.85);

    char *mean = Form("Sig mean = %.2f ", mean1);
    char *rms = Form("Sig rms = %.2f ", rms1);

    char *bmean = Form("Bkg mean = %.2f ", mean2);
    char *brms = Form("Bkg rms = %.2f ", rms2);

    leg2->SetBorderSize(0);
    leg2->SetHeader("FD group -SiFCC");
    leg2->AddEntry(sig, "Z'(5TeV)#rightarrowW^{+}W^{-} #rightarrow 2 subjet", "l");
    leg2->AddEntry((TObject *)0, mean, "");
    leg2->AddEntry((TObject *)0, rms, "");
    leg2->AddEntry((TObject *)0, "", "");
    leg2->AddEntry(bkg, "Z'(5TeV)#rightarrowq#bar{q} #rightarrow 1 subjet", "l");
    leg2->AddEntry((TObject *)0, bmean, "");
    leg2->AddEntry((TObject *)0, brms, "");

    leg2->Draw("same");
    c1->Print("LeadingTime.png");
}
void Velocity()
{
    auto c1 = new TCanvas("c", "BPRE", 10, 10, 800, 600);
    double mean1, rms1, mean2, rms2;
    TFile *file1 = new TFile("../output/tev5mm_pythia6_zprime5tev_ww_0199_hepsim.slcio.root");
    TFile *file2 = new TFile("../output/tev5mm_pythia6_zprime5tev_qq_0199_hepsim.slcio.root");

    //TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Leading"));
    //TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Leading"));
    TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Trailing_V"));
    TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Trailing_V"));

    mean1 = sig->GetMean();
    rms1 = sig->GetRMS();

    mean2 = bkg->GetMean();
    rms2 = bkg->GetRMS();

    SetDrawNcuOpt(sig, kRed, "Trailing Velocity  - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    SetDrawNcuOpt(bkg, kBlue, "Trailing Velocity - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    SetRootStyle();

    sig->Scale(1.0 / sig->Integral());
    bkg->Scale(1.0 / bkg->Integral());

    sig->GetXaxis()->SetRangeUser(6.5, 24.5);
    sig->GetYaxis()->SetRangeUser(0, 0.3);
    bkg->GetXaxis()->SetRangeUser(6.5, 24.5);
    bkg->GetYaxis()->SetRangeUser(0, 0.3);

    bkg->DrawNormalized("hist");
    sig->DrawNormalized("same&&hist");

    TLegend *leg2 = new TLegend(0.5, 0.45, 0.89, 0.85);

    char *mean = Form("Sig mean = %.2f ", mean1);
    char *rms = Form("Sig rms = %.2f ", rms1);

    char *bmean = Form("Bkg mean = %.2f ", mean2);
    char *brms = Form("Bkg rms = %.2f ", rms2);

    leg2->SetBorderSize(0);
    leg2->SetHeader("FD group -SiFCC");
    leg2->AddEntry(sig, "Z'(5TeV)#rightarrowW^{+}W^{-} #rightarrow 2 subjet", "l");
    leg2->AddEntry((TObject *)0, mean, "");
    leg2->AddEntry((TObject *)0, rms, "");
    leg2->AddEntry((TObject *)0, "", "");
    leg2->AddEntry(bkg, "Z'(5TeV)#rightarrowq#bar{q} #rightarrow 1 subjet", "l");
    leg2->AddEntry((TObject *)0, bmean, "");
    leg2->AddEntry((TObject *)0, brms, "");

    leg2->Draw("same");
    c1->Print("TrailingV.png");
}

void Momentum()
{
    auto c1 = new TCanvas("c", "BPRE", 10, 10, 800, 600);
    double mean1, rms1, mean2, rms2;
    TFile *file1 = new TFile("../output/tev5mm_pythia6_zprime5tev_ww_0199_hepsim.slcio.root");
    TFile *file2 = new TFile("../output/tev5mm_pythia6_zprime5tev_qq_0199_hepsim.slcio.root");

    //TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Leading"));
    //TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Leading"));
    TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Trailing_P"));
    TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Trailing_P"));

    mean1 = sig->GetMean();
    rms1 = sig->GetRMS();

    mean2 = bkg->GetMean();
    rms2 = bkg->GetRMS();

    SetDrawNcuOpt(sig, kRed, "Momentum  - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    SetDrawNcuOpt(bkg, kBlue, "Momentum - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    SetRootStyle();

    sig->Scale(1.0 / sig->Integral());
    bkg->Scale(1.0 / bkg->Integral());

    sig->GetXaxis()->SetRangeUser(6.5, 24.5);
    sig->GetYaxis()->SetRangeUser(0, 0.3);
    bkg->GetXaxis()->SetRangeUser(6.5, 24.5);
    bkg->GetYaxis()->SetRangeUser(0, 0.3);

    bkg->DrawNormalized("hist");
    sig->DrawNormalized("same&&hist");

    TLegend *leg2 = new TLegend(0.5, 0.45, 0.89, 0.85);

    char *mean = Form("Sig mean = %.2f ", mean1);
    char *rms = Form("Sig rms = %.2f ", rms1);

    char *bmean = Form("Bkg mean = %.2f ", mean2);
    char *brms = Form("Bkg rms = %.2f ", rms2);

    leg2->SetBorderSize(0);
    leg2->SetHeader("FD group -SiFCC");
    leg2->AddEntry(sig, "Z'(5TeV)#rightarrowW^{+}W^{-} #rightarrow 2 subjet", "l");
    leg2->AddEntry((TObject *)0, mean, "");
    leg2->AddEntry((TObject *)0, rms, "");
    leg2->AddEntry((TObject *)0, "", "");
    leg2->AddEntry(bkg, "Z'(5TeV)#rightarrowq#bar{q} #rightarrow 1 subjet", "l");
    leg2->AddEntry((TObject *)0, bmean, "");
    leg2->AddEntry((TObject *)0, brms, "");

    leg2->Draw("same");
    c1->Print("Momentum.png");
}
