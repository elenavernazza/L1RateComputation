#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TString.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <TPaletteAxis.h>
#include <TLatex.h>
#include <cstdlib>

// for the ZeroBias sample 362616 we can look at this page
// https://cmsoms.cern.ch/cms/runs/lumisection?cms_run=362616&cms_run_sequence=GLOBAL-RUN
// thisLumiRun = 20.5 → 2.050E34;

using namespace std;

void PlotRate_2D_ptj1_mjj(Int_t fixed_bin_y, Int_t fixed_bin_w, THnF* rates_4D, TString output, Int_t bin0, Int_t xmin0, Int_t xmax0, Int_t bin2, Int_t xmin2, Int_t xmax2, Int_t xmin3) 
{
    // Plot 2D_rate with a fixed muonPt value and jetPt2 value 
    TH2F* rate_2D_ptj1_mjj = new TH2F("rate_2D_ptj1_mjj","rate_2D_ptj1_mjj", bin0, xmin0, xmax0, bin2, xmin2, xmax2);
    Int_t fixed_cut_y = rates_4D->GetAxis(1)->GetBinLowEdge(fixed_bin_y);
    Int_t fixed_cut_w = xmin3 + fixed_bin_w - 1;

    Float_t rate_value = 0;
    Int_t vec_bin[4] = {0, 0, 0, 0};
    for (Int_t I_bin_x = 1 ; I_bin_x <= bin0 ; ++I_bin_x)
    {
        for (Int_t I_bin_z = 1 ; I_bin_z <= bin2 ; ++I_bin_z)
        {
            Int_t vec_bin[4] = {I_bin_x, fixed_bin_y, I_bin_z, fixed_bin_w};
            rate_value = rates_4D->GetBinContent(vec_bin);
            rate_2D_ptj1_mjj->SetBinContent(I_bin_x, I_bin_z, rate_value);
        }
    }

    TCanvas* c_2D_1 = new TCanvas("c_2D_1","c_2D_1",700.,550.);
    c_2D_1->cd();
    c_2D_1->SetRightMargin(0.16); // important for labels going outside the canvas!
    rate_2D_ptj1_mjj->GetXaxis()->SetTitle("p_{T}^{j1} > X [GeV]");
    rate_2D_ptj1_mjj->GetXaxis()->SetTitleOffset(1.3);
    rate_2D_ptj1_mjj->GetYaxis()->SetTitle("M_{jj} > Y [GeV]");
    rate_2D_ptj1_mjj->GetYaxis()->SetTitleOffset(1.3);
    rate_2D_ptj1_mjj->GetZaxis()->SetTitle("rate [kHz]");
    rate_2D_ptj1_mjj->GetZaxis()->SetTitleOffset(1.3);
    rate_2D_ptj1_mjj->SetStats(0);
    rate_2D_ptj1_mjj->Draw("COLZ");
    rate_2D_ptj1_mjj->SetMinimum(0);
    rate_2D_ptj1_mjj->SetMaximum(5);
    rate_2D_ptj1_mjj->SetTitle("");

    TLatex Tex11;
    Tex11.SetTextSize(0.03);
    Tex11.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex11.Draw("same");

    TLatex Tex22;
    Tex22.SetTextSize(0.035);
    Tex22.SetTextAlign(31);
    Tex22.DrawLatexNDC(0.88,0.91,"(14 TeV)");
    Tex22.Draw("same");

    TLatex Tex33;
    Tex33.SetTextSize(0.035);
    Tex33.SetTextAlign(12);
    TString printout = Form("Zero Bias L1 rate p_{T}^{jet2} > %i GeV && p_{T}^{#mu} > %i GeV", (int)fixed_cut_y, (int)fixed_cut_w);
    Tex33.DrawLatexNDC(0.25,0.96, printout);
    Tex33.Draw("same");

    rate_2D_ptj1_mjj->Write();

    TString OutFileName;
    c_2D_1->SaveAs(Form(output + "/total_rate_2D_ptj1_X_ptj2_%i_mjj_Y_ptmu_%i.png", (int)fixed_cut_y, (int)fixed_cut_w));
    c_2D_1->SaveAs(Form(output + "/total_rate_2D_ptj1_X_ptj2_%i_mjj_Y_ptmu_%i.pdf", (int)fixed_cut_y, (int)fixed_cut_w));

    TH2F *clone_rate_2D_ptj1_mjj = (TH2F*) rate_2D_ptj1_mjj->Clone();
    clone_rate_2D_ptj1_mjj->SetName("clone_rate_2D_ptj1_mjj");
    // draw contour line
    double contours[2];
    contours[0] = 0.95;
    contours[1] = 1.05;
    clone_rate_2D_ptj1_mjj->SetContour(2,contours);
    clone_rate_2D_ptj1_mjj->Draw("cont3 list same");
    clone_rate_2D_ptj1_mjj->SetLineColor(kRed);
    clone_rate_2D_ptj1_mjj->SetLineStyle(kRed);

    c_2D_1->SaveAs(Form(output + "/total_rate_2D_ptj1_X_ptj2_%i_mjj_Y_ptmu_%i_line.png", (int)fixed_cut_y, (int)fixed_cut_w));
    c_2D_1->SaveAs(Form(output + "/total_rate_2D_ptj1_X_ptj2_%i_mjj_Y_ptmu_%i_line.pdf", (int)fixed_cut_y, (int)fixed_cut_w));
    c_2D_1->Close();
}

void Plot_2D()
{

    TString Indir = "/afs/cern.ch/user/e/evernazz/L1Studies/RateComputation/CMSSW_13_0_0_pre2/src/L1Trigger/L1TNtuples/TotalRate/Run362616_L1_DoubleJet_X_Y_Mass_MinZ_30_30_MuW_OpenQual";
    TString FileName_rate = Indir + "/Rates_4D.root";
    cout << "\nReading rates from :\n" << FileName_rate << endl;

    TFile f (FileName_rate.Data(),"READ");
    THnF* rates_4D = (THnF*)f.Get("rates_4D");

    Int_t bins[4] = {rates_4D->GetAxis(0)->GetNbins(), rates_4D->GetAxis(1)->GetNbins(), rates_4D->GetAxis(2)->GetNbins(), rates_4D->GetAxis(3)->GetNbins()};
    Double_t xmin[4] = {rates_4D->GetAxis(0)->GetBinLowEdge(1), rates_4D->GetAxis(1)->GetBinLowEdge(1), rates_4D->GetAxis(2)->GetBinLowEdge(1), rates_4D->GetAxis(3)->GetBinLowEdge(1)};
    Double_t xmax[4] = {rates_4D->GetAxis(0)->GetBinUpEdge(rates_4D->GetAxis(0)->GetNbins()), \
                        rates_4D->GetAxis(1)->GetBinUpEdge(rates_4D->GetAxis(1)->GetNbins()), \
                        rates_4D->GetAxis(2)->GetBinUpEdge(rates_4D->GetAxis(2)->GetNbins()), \
                        rates_4D->GetAxis(3)->GetBinUpEdge(rates_4D->GetAxis(3)->GetNbins())};

    Float_t Total = 1.21205e+07;

    float nb = 2450.;
    float scale_rate = 0.001*(nb*11245.6);
    float thisLumiRun = 2.050E34;
    float scaleToLumi = 2.000E34;
    float scale_lumi = scaleToLumi/thisLumiRun;

    cout << "Total number of events = " << Total << endl;

    TString output2D = Indir + "/rate_2D_plots";
    cout << "Saving plots in :\n" << output2D << endl;
    int ret1 = system("mkdir -p "+output2D);

    // --------------- Plotting ---------------

    TString rootname2D = Indir+"/total_rate_2D.root";
    TFile *f2D=TFile::Open(rootname2D, "recreate");
    f2D->cd();

    PlotRate_2D_ptj1_mjj(1, 1, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(1, 2, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(1, 3, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(1, 4, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(1, 5, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(1, 6, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(1, 7, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(1, 8, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 1, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 2, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 3, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 4, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 5, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 6, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 7, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(2, 8, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 1, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 2, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 3, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 4, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 5, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 6, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 7, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
    PlotRate_2D_ptj1_mjj(3, 8, rates_4D, output2D, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);

    f2D->Close();

}
