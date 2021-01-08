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
void Momentum();
void Velocity();

void Draw()
{
    //Choose Which function you want to run
    //Time();
    //Velocity();
    Momentum();
}

void Time()
{
    //Create the TCanvas first
    auto c1 = new TCanvas("c", "BPRE", 10, 10, 800, 600);
    double mean1, rms1, mean2, rms2;

    //Input the root file
    TFile *file1 = new TFile("../output/tev5mm_pythia6_zprime5tev_ww_0199_hepsim.slcio.root");
    TFile *file2 = new TFile("../output/tev5mm_pythia6_zprime5tev_qq_0199_hepsim.slcio.root");

    //TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Leading"));
    //TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Leading"));
    TH1F *sig = ((TH1F *)file1->Get("Timing_detector_next_to_trailing"));
    TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_next_to_trailing"));

    mean1 = sig->GetMean();
    rms1 = sig->GetRMS();

    mean2 = bkg->GetMean();
    rms2 = bkg->GetRMS();

    SetDrawNcuOpt(sig, kRed, "Time of flight - collision point to HCAL(Next to Trailing)", "T [ns]", "Arbitrary number");
    SetDrawNcuOpt(bkg, kBlue, "Time of flight - collision point to HCAL(Next to Trailing)", "T [ns]", "Arbitrary number");
    //SetDrawNcuOpt(sig, kRed, "Time of flight - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    //SetDrawNcuOpt(bkg, kBlue, "Time of flight - collision point to HCAL(Trailing)", "T [ns]", "Arbitrary number");
    SetRootStyle();

    //Normalize signal and background
    sig->Scale(1.0 / sig->Integral());
    bkg->Scale(1.0 / bkg->Integral());

    //customized your signal and background range
    sig->GetXaxis()->SetRangeUser(6.5, 24.5);
    sig->GetYaxis()->SetRangeUser(0, 0.3);
    bkg->GetXaxis()->SetRangeUser(6.5, 24.5);
    bkg->GetYaxis()->SetRangeUser(0, 0.3);

    bkg->Draw("hist");
    sig->Draw("same&&hist");

    //customized your TLegend
    TLegend *leg2 = new TLegend(0.5, 0.45, 0.89, 0.85);

    char *mean = Form("Sig mean = %.2f ", mean1);
    char *rms = Form("Sig rms = %.2f ", rms1);

    char *bmean = Form("Bkg mean = %.2f ", mean2);
    char *brms = Form("Bkg rms = %.2f ", rms2);

    //Create the TLegend
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
    //c1->Print("TrailingTime.png");
    c1->Print("Next_to_TrailingTime.png");
}
void Velocity()
{
    auto c1 = new TCanvas("c", "BPRE", 10, 10, 800, 600);
    c1->SetLogx();

    double mean1, rms1, mean2, rms2;
    TFile *file1 = new TFile("../output/tev5mm_pythia6_zprime5tev_ww_0199_hepsim.slcio.root");
    TFile *file2 = new TFile("../output/tev5mm_pythia6_zprime5tev_qq_0199_hepsim.slcio.root");

    //TFile *file1 = new TFile("../output/testsig.root");
    //TFile *file2 = new TFile("../output/testbkg.root");

    //TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Leading"));
    //TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Leading"));
    TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Trailing_V"));
    TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Trailing_V"));

    mean1 = sig->GetMean();
    rms1 = sig->GetRMS();

    mean2 = bkg->GetMean();
    rms2 = bkg->GetRMS();

    SetDrawNcuOpt(sig, kRed, "Trailing Velocity  - collision point to HCAL(Trailing)", "V (log)", "Arbitrary number");
    SetDrawNcuOpt(bkg, kBlue, "Trailing Velocity - collision point to HCAL(Trailing)", "V (log)", "Arbitrary number");
    SetRootStyle();

    //gPad->SetLogx(1);
    sig->Scale(1.0 / sig->Integral());
    bkg->Scale(1.0 / bkg->Integral());

    //sig->GetXaxis()->SetRangeUser(8.4, 8.5);
    //sig->GetYaxis()->SetRangeUser(0, 0.8);
    //bkg->GetXaxis()->SetRangeUser(8.4, 8.5);
    //bkg->GetYaxis()->SetRangeUser(0, 0.8);

    //gPad->SetLogy();

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

    TH1F *sig = ((TH1F *)file1->Get("Timing_detector_Trailing_P"));
    TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_Trailing_P"));

    //TH1F *sig = ((TH1F *)file1->Get("Timing_detector_next_to_trailing_P"));
    //TH1F *bkg = ((TH1F *)file2->Get("Timing_detector_next_to_trailing_P"));

    mean1 = sig->GetMean();
    rms1 = sig->GetRMS();

    mean2 = bkg->GetMean();
    rms2 = bkg->GetRMS();

    SetDrawNcuOpt(sig, kRed, "Momentum  - collision point to HCAL(Trailing)", "P [GeV]", "Arbitrary number");
    SetDrawNcuOpt(bkg, kBlue, "Momentum - collision point to HCAL(Trailing)", "P [GeV]", "Arbitrary number");

    //SetDrawNcuOpt(sig, kRed, "Momentum  - collision point to HCAL(Next to Trailing)", "P [GeV]", "Arbitrary number");
    //SetDrawNcuOpt(bkg, kBlue, "Momentum - collision point to HCAL(Next to Trailing)", "P [GeV]", "Arbitrary number");
    SetRootStyle();

    sig->Scale(1.0 / sig->Integral());
    bkg->Scale(1.0 / bkg->Integral());

    //sig->GetXaxis()->SetRangeUser(0, 20);
    //sig->GetYaxis()->SetRangeUser(0, 0.3);
    //bkg->GetXaxis()->SetRangeUser(0, 20);
    //bkg->GetYaxis()->SetRangeUser(0, 0.3);

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
    //c1->Print("Next_to_Momentum.png");
    c1->Print("Momentum.png");
}
