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

// Header file for the classes stored in the TTree if any.
#include "interface/L1AnalysisL1UpgradeDataFormat.h"
// #include "interface/L1AnalysisL1CaloClusterDataFormat.h"
#include "interface/L1AnalysisEventDataFormat.h"

// for the ZeroBias sample 362616 we can look at this page
// https://cmsoms.cern.ch/cms/runs/lumisection?cms_run=362616&cms_run_sequence=GLOBAL-RUN
// thisLumiRun = 20.5 → 2.050E34;

using namespace std;

void PlotRate_2D_ptj1_mjj(Int_t fixed_bin_y, Int_t fixed_bin_w, THnF* rates_4D, TString output, Int_t bin0, Int_t xmin0, Int_t xmax0, Int_t bin2, Int_t xmin2, Int_t xmax2, Int_t xmin3) 
{
    // Compute 2D_rate with a fixed muonPt value and jetPt2 value 
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
    rate_2D_ptj1_mjj->Write();

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

void Rate()
{

    Int_t bins[4] = {14, 14, 30, 12};
    Double_t xmin[4] = {30., 30., 200., 3.};
    Double_t xmax[4] = {100., 100., 800., 15.};

    Float_t Max_Events = -1;
    Float_t Total = 0;
    Float_t Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 = 0.;

    cout << "Begin loop" << endl;

    TString path = "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/elfontan/condor/2022EphZB_run362616_126X";
    TString output = "/afs/cern.ch/user/e/evernazz/L1Studies/RateComputation/CMSSW_13_0_0_pre2/src/L1Trigger/L1TNtuples/TotalRate/Run362616_L1_DoubleJet_X_Y_Mass_MinZ_30_30_MuW_OpenQual";
    int ret = system("mkdir -p "+output);

    vector <TH3F*> ptjet1_ptjet2_jetmass_ptmu;
    for (Int_t m = 0; m < bins[3]; m++) ptjet1_ptjet2_jetmass_ptmu.push_back(nullptr); 
    for (Int_t m = 0; m < bins[3]; m++)
    {
        int num = xmin[3] + m;
        TString histname = Form("ptjet1_ptjet2_jetmass_ptmu_%d", num);
        ptjet1_ptjet2_jetmass_ptmu.at(m) = new TH3F(histname, histname, bins[0], xmin[0], xmax[0], bins[1], xmin[1], xmax[1], bins[2], xmin[2], xmax[2]);
    }

    float nb = 2450.;
    float scale_rate = 0.001*(nb*11245.6);
    float thisLumiRun = 2.050E34;
    float scaleToLumi = 2.000E34;
    float scale_lumi = scaleToLumi/thisLumiRun;

    TString NameL1EventTree = "l1EventTree/L1EventTree";
    TString NameL1UpgradeTree = "l1UpgradeTree/L1UpgradeTree";

    TChain dataL1Event(NameL1EventTree.Data());
    TChain dataL1Upgrade(NameL1UpgradeTree.Data());

    TString FileName_in = path + "/*.root";
    dataL1Event.Add(FileName_in.Data());
    dataL1Upgrade.Add(FileName_in.Data());

    dataL1Event.SetMakeClass(1);
    dataL1Upgrade.SetMakeClass(1);

    ULong_t                     event;
    UInt_t                      in_l1tnJets;
    vector<float>               in_l1tPtJet;
    vector<float>               in_l1tEtaJet;
    vector<float>               in_l1tPhiJet;
    vector<short int>           in_l1tBxJet;
    vector<float>               in_l1tMuPt;
    vector<float>               in_l1tMuEta;
    vector<float>               in_l1tMuPhi;
    vector<short int>           in_l1tMuBx;
    vector<unsigned short int>  in_l1tMuQual;

    TBranch *b_L1Event_event;
    TBranch *b_L1Upgrade_l1tnJets;
    TBranch *b_L1Upgrade_l1tPtJet;
    TBranch *b_L1Upgrade_l1tPhiJet;
    TBranch *b_L1Upgrade_l1tEtaJet;
    TBranch *b_L1Upgrade_l1tBxJet;
    TBranch *b_L1Upgrade_l1tMuPt;
    TBranch *b_L1Upgrade_l1tMuEta;
    TBranch *b_L1Upgrade_l1tMuPhi;
    TBranch *b_L1Upgrade_l1tMuBx;
    TBranch *b_L1Upgrade_l1tMuQual;

    dataL1Event.SetBranchAddress("event",      &event,         &b_L1Event_event);
    dataL1Upgrade.SetBranchAddress("nJets",    &in_l1tnJets,   &b_L1Upgrade_l1tnJets);
    dataL1Upgrade.SetBranchAddress("jetEt",    &in_l1tPtJet,   &b_L1Upgrade_l1tPtJet);
    dataL1Upgrade.SetBranchAddress("jetEta",   &in_l1tEtaJet,  &b_L1Upgrade_l1tEtaJet);
    dataL1Upgrade.SetBranchAddress("jetPhi",   &in_l1tPhiJet,  &b_L1Upgrade_l1tPhiJet);
    dataL1Upgrade.SetBranchAddress("jetBx",    &in_l1tBxJet,   &b_L1Upgrade_l1tBxJet);
    dataL1Upgrade.SetBranchAddress("muonEt",   &in_l1tMuPt,    &b_L1Upgrade_l1tMuPt);
    dataL1Upgrade.SetBranchAddress("muonEta",  &in_l1tMuEta,   &b_L1Upgrade_l1tMuEta);
    dataL1Upgrade.SetBranchAddress("muonPhi",  &in_l1tMuPhi,   &b_L1Upgrade_l1tMuPhi);
    dataL1Upgrade.SetBranchAddress("muonBx",   &in_l1tMuBx,    &b_L1Upgrade_l1tMuBx);
    dataL1Upgrade.SetBranchAddress("muonQual", &in_l1tMuQual,  &b_L1Upgrade_l1tMuQual);

    int nEntries = dataL1Upgrade.GetEntries();
    cout << "TChain has " << nEntries << " events" << endl;
    if (Max_Events == -1)
    {
        Max_Events = nEntries;
    }
    cout << "Analysing " << Max_Events << " events\n" << endl;

    for(int i = 0 ; i < Max_Events ; ++i)
    {
        if (i%100000 == 0) cout << "Analysing entry " << i << endl;
        scale_lumi = scaleToLumi/thisLumiRun;
        dataL1Event.GetEntry(i);
        dataL1Upgrade.GetEntry(i);
        ++ Total;

        bool L1_objects_VBF = in_l1tPtJet.size() > 1 && in_l1tMuPt.size() > 0 ;

        TLorentzVector myGoodOnlineMuon;
        TLorentzVector myGoodOnlineJet1;
        TLorentzVector myGoodOnlineJet2;
        TLorentzVector myGoodOnlineDiJet;
        float myGoodOnlineMjj = -1.;

        if (L1_objects_VBF)
        {

            // take the first (most enegrgetic) muon passing the open quality requirement
            for (UInt_t iMuon = 0; iMuon < in_l1tMuPt.size(); ++iMuon)
            {
                // IMPORTANT: only select objects in BX == 0
                if (in_l1tMuBx.at(iMuon) != 0) continue;
                if (in_l1tMuQual.at(iMuon) < 4) continue; // open quality: 4, 5, 6, 7, 8, 9, 10, 11, 12,13, 14, 15
                myGoodOnlineMuon.SetPtEtaPhiM(in_l1tMuPt.at(iMuon),in_l1tMuEta.at(iMuon),in_l1tMuPhi.at(iMuon),0.105);
                break;
            }

            // loop among the pairs of jets (80,30) giving the highest mjj
            for (UInt_t iL1Jet1 = 0 ; iL1Jet1 < in_l1tPtJet.size() ; ++iL1Jet1)
            {
                // IMPORTANT: only select objects in BX == 0
                if (in_l1tBxJet.at(iL1Jet1) != 0) continue;
                TLorentzVector myOnlineJet1;
                myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet.at(iL1Jet1),in_l1tEtaJet.at(iL1Jet1),in_l1tPhiJet.at(iL1Jet1),0.);
                if (myOnlineJet1.Pt() < 30) continue;

                for (UInt_t jL1Jet2 = iL1Jet1 + 1 ; jL1Jet2 < in_l1tPtJet.size() ; ++jL1Jet2)
                {
                    // IMPORTANT: only select objects in BX == 0
                    if (in_l1tBxJet.at(jL1Jet2) != 0) continue;
                    TLorentzVector myOnlineJet2;
                    myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet.at(jL1Jet2),in_l1tEtaJet.at(jL1Jet2),in_l1tPhiJet.at(jL1Jet2),0.);
                    if (myOnlineJet2.Pt() < 30) continue;
                    
                    TLorentzVector myOnlineDiJet = myOnlineJet1 + myOnlineJet2;
                    float myOnlineMjj = myOnlineDiJet.M();

                    if (myOnlineMjj > myGoodOnlineMjj)
                    {
                        myGoodOnlineMjj = myOnlineMjj;
                        myGoodOnlineJet1 = myOnlineJet1;
                        myGoodOnlineJet2 = myOnlineJet2;
                        myGoodOnlineDiJet = myOnlineJet1 + myOnlineJet2;
                    }
                }  
            }

            for (Int_t n_cut = 0; n_cut < bins[3]; n_cut++)
            {
                int muon_cut_min = n_cut + xmin[3];
                int muon_cut_max = n_cut + 1 + xmin[3];
                if (n_cut == 0)
                {
                    if (myGoodOnlineMuon.Pt() < muon_cut_max)
                    {
                        ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(myGoodOnlineJet1.Pt(), myGoodOnlineJet2.Pt(), myGoodOnlineDiJet.M(), scale_lumi);
                    }
                }
                else if (n_cut == bins[3]-1)
                {
                    if (myGoodOnlineMuon.Pt() >= muon_cut_min)
                    {
                        ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(myGoodOnlineJet1.Pt(), myGoodOnlineJet2.Pt(), myGoodOnlineDiJet.M(), scale_lumi);
                    }
                }
                else
                {
                    if (myGoodOnlineMuon.Pt() >= muon_cut_min && myGoodOnlineMuon.Pt() < muon_cut_max)
                    {
                        ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(myGoodOnlineJet1.Pt(), myGoodOnlineJet2.Pt(), myGoodOnlineDiJet.M(), scale_lumi);
                    }
                }
            }
        }

    }

    cout << "Total number of events = " << Total << endl;

    // Compute 4D histogram, where at each combination of jetPt1, jetPt2, Mjj and muonPt corresponds a rate
    THnF* rates_4D = new THnF("rates_4D","rates_4D", 4, bins, xmin, xmax);
    THnF* acceptances_4D = new THnF("acceptances_4D","acceptances_4D", 4, bins, xmin, xmax);
    float integral_xyzw = -1.;
    int combinations = 0;

    for (Int_t i_bin_x = 1 ; i_bin_x <= bins[0] ; ++i_bin_x)
    {
        for (Int_t i_bin_y = 1 ; i_bin_y <= bins[1] ; ++i_bin_y)
        {
            for (Int_t i_bin_z = 1 ; i_bin_z <= bins[2] ; ++i_bin_z)
            {

                // compute the 3D integral over ptj1, ptj2, mjj for each ptmu bin (from bin 1 to the overflow bin)
                vector<float> int_xyz ;
                for (Int_t w = 0; w < bins[3] ; w++)
                {
                    int_xyz.push_back(ptjet1_ptjet2_jetmass_ptmu.at(w)->Integral(i_bin_x,bins[0]+1,i_bin_y,bins[1]+1,i_bin_z,bins[2]+1)/Total*scale_rate);
                }

                for (Int_t i_bin_w = 1 ; i_bin_w <= bins[3] ; ++i_bin_w)
                {
                    // compute the 4D integral over ptmu 
                    integral_xyzw = 0;
                    for (Int_t i_vec_w = i_bin_w - 1; i_vec_w < bins[3]; ++i_vec_w) integral_xyzw += int_xyz.at(i_vec_w);

                    // fill in the rate_4d histogram to be saved
                    Int_t v_bin[4] = {i_bin_x, i_bin_y, i_bin_z, i_bin_w};
                    rates_4D->SetBinContent(v_bin, integral_xyzw);
                    acceptances_4D->SetBinContent(v_bin, integral_xyzw*Total/scale_rate);
                    if (integral_xyzw < 1.05 && integral_xyzw > 0.95 && i_bin_x >= i_bin_y) ++ combinations;

                }
            }
        }
    }

    // --------------- Plotting ---------------

    TString rootname2D = output+"/total_rate_2D.root";
    TFile *f2D=TFile::Open(rootname2D, "recreate");
    f2D->cd();

    TString output2D = output + "/rate_2D_plots";
    int ret1 = system("mkdir -p "+output2D);

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

    vector <TCanvas*> canvas_mu;
    for (Int_t m =0; m < bins[3]; m++)
    {
        canvas_mu.push_back(nullptr);
    }

    TString rootname = output+"/rate_3D_ptjet1_ptjet2_jetmass_ptmu.root";
    TFile *f=TFile::Open(rootname, "recreate");
    f->cd();

    for (Int_t m = 0; m < bins[3]; m++)
    {
        int num = xmin[3] + m;
        TString canvname = Form("/rate_3D_ptjet1_ptjet2_jetmass_ptmu_%d", num);
        TString pngname = Form(output+"/rate_3D_ptjet1_ptjet2_jetmass_ptmu_%d.png", num);
        TString pdfname = Form(output+"/rate_3D_ptjet1_ptjet2_jetmass_ptmu_%d.pdf", num);
        canvas_mu.at(m) = new TCanvas(canvname,canvname,900.,650.);
        ptjet1_ptjet2_jetmass_ptmu.at(m)->SetTitle("");
        ptjet1_ptjet2_jetmass_ptmu.at(m)->Draw("BOX2 Z");
        ptjet1_ptjet2_jetmass_ptmu.at(m)->GetXaxis()->SetTitle("p_{T}^{Jet1} > X [GeV]");
        ptjet1_ptjet2_jetmass_ptmu.at(m)->GetXaxis()->SetTitleOffset(1.8);
        ptjet1_ptjet2_jetmass_ptmu.at(m)->GetYaxis()->SetTitle("p_{T}^{Jet2} > Y [GeV]");
        ptjet1_ptjet2_jetmass_ptmu.at(m)->GetYaxis()->SetTitleOffset(1.8);
        ptjet1_ptjet2_jetmass_ptmu.at(m)->GetZaxis()->SetTitle("m_{jj} > Z [GeV]");
        ptjet1_ptjet2_jetmass_ptmu.at(m)->GetZaxis()->SetTitleOffset(1.4);
        ptjet1_ptjet2_jetmass_ptmu.at(m)->SetMinimum(0);
        ptjet1_ptjet2_jetmass_ptmu.at(m)->SetMaximum(500);
        ptjet1_ptjet2_jetmass_ptmu.at(m)->Write();
        gPad->Update();
        TPaletteAxis* pal = (TPaletteAxis*) ptjet1_ptjet2_jetmass_ptmu.at(m)->GetListOfFunctions()->FindObject("palette");
        pal->SetX1NDC(0.926); // It doesn't work (Why??)
        // ptjet1_ptjet2_jetmass_ptmu.at(m)->SaveAs(rootname);

        TLatex Tex4;
        Tex4.SetTextSize(0.035);
        Tex4.DrawLatexNDC(0.3,0.96,Form("Zero Bias L1 acceptance p_{T}^{#mu} > %d GeV", num));
        Tex4.Draw("same");

        canvas_mu.at(m)->SaveAs(pngname);
        canvas_mu.at(m)->SaveAs(pdfname);
        canvas_mu.at(m)->Close();
    }

    f->Close();

    rates_4D->SaveAs(output+"/Rates_4D.root");
    acceptances_4D->SaveAs(output+"/Acceptances_4D.root");

    cout << "\nNumber of bins = " << rates_4D->GetNbins() << endl;
    cout << "Number of good combinations = " << combinations << endl;

}
