// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <THStack.h>
#include <TFitResultPtr.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TEfficiency.h>
#include <TMath.h>
#include "TLorentzVector.h"
#include <TRandom3.h>
#include "TSystem.h"
#include "TROOT.h"
#include <TGraph2D.h>
#include <TRandom.h>
#include <TF2.h>

// C++ includes
#include <iostream>
#include <iomanip>
using namespace std;
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector> // Need this for generate dictionary for nested vectors
using namespace std;

// Include customized functions and constants
#include "Helpers.h"

float vSize = 30.;

// Pick six evenly spaced points in each non-dead region
vector<double> generatePoints(double start)
{
  vector<double> points;
  double step = 7.0; // Interval between points
  double current;

  if (start == -300) {
      current = start + 1.45; // Adjusted start for the edge case
  } else if (start == 256.55) {
      current = start + 7.0; // Start for the symmetric edge case
  } else {
      current = start + step; // Regular start
  }

  for (int i = 0; i < 6; ++i) {
      points.push_back(current + i * step);
  }
  return points;
}

// Draw ND hadronic hits plot
void NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA()
{
  gROOT->Reset();
  //gStyle->SetOptStat(0); // Remove Stat Box
  //File with coefficients histogram
  TFile* FileWithCoeffs = new TFile("FileWithCoeffsNuMu.root", "READ");
  FileWithCoeffs->cd();
  TH1D* CoefficientsHist = (TH1D*) FileWithCoeffs->Get("CoeffPRISMUpTo3mOA");
  CoefficientsHist->SetDirectory(0);
  FileWithCoeffs->Close();
  // Input FDroot file
  TString FileIn = "RootFiles/FDGeoEff_72040717_0.root";
  //
  // Read branch from input trees
  //
  // Read effValues
  TChain *t_effValues = new TChain("effValues");
  t_effValues->Add(FileIn.Data());
  Int_t iwritten;
  Double_t ND_LAr_dtctr_pos;
  Double_t ND_LAr_vtx_pos;
  Double_t ND_GeoEff;
  int validThrows;
  int NPassedThrows;

  vector<uint64_t> *ThrowResults = nullptr; //all throw results
  vector<float> *VetoEnergyEventsPass = nullptr; // veto energy of events that pass the throw
  vector<float> *TrimEnergyEventsPass = nullptr;
  vector<float> *TotalEnergyEventsPass = nullptr;
  vector<float> *VtxYEventsPass = nullptr;
  vector<float> *VtxZEventsPass = nullptr;
  t_effValues->SetBranchAddress("iwritten",         &iwritten);
  t_effValues->SetBranchAddress("ND_LAr_dtctr_pos",   &ND_LAr_dtctr_pos);
  t_effValues->SetBranchAddress("ND_LAr_vtx_pos",       &ND_LAr_vtx_pos);
  t_effValues->SetBranchAddress("ND_GeoEff",        &ND_GeoEff);
  t_effValues->SetBranchAddress("ThrowResults",      &ThrowResults);
  t_effValues->SetBranchAddress("validThrows",  &validThrows);
  t_effValues->SetBranchAddress("NPassedThrows", &NPassedThrows);
  t_effValues->SetBranchAddress("VetoEnergyEventsPass", &VetoEnergyEventsPass);
  t_effValues->SetBranchAddress("TrimEnergyEventsPass",  &TrimEnergyEventsPass);
  t_effValues->SetBranchAddress("TotalEnergyEventsPass",  &TotalEnergyEventsPass);
  t_effValues->SetBranchAddress("VtxYEventsPass",  &VtxYEventsPass);
  t_effValues->SetBranchAddress("VtxZEventsPass",  &VtxZEventsPass);
  //hardon info
  TChain *t_effTree = new TChain("effTreeND");
  Float_t totEnergyFDatND_f;
  t_effTree->Add(FileIn.Data());
  int ND_Sim_n_hadronic_Edep_b;
  vector<float> *HadronHitEdeps =0; // Hadron hit segment energy deposits [MeV]
  // vector<vector<float>> *CurrentThrowDepsX = 0; // Coordinates of hadron hits X after random throws
  // vector<vector<float>> *CurrentThrowDepsY =0; // Coordinates of hadron hits Y after random throws
  // vector<vector<float>> *CurrentThrowDepsZ = 0; // Coordinates of hadron hits Z after random throws
  vector<float> *CurrentThrowTotE = 0;
  vector<vector<float>> *ND_OffAxis_Sim_hadronic_hit_xyz=0; // coordinates of hadron hits before random throws

  double ND_Gen_numu_E; // Energy of generator level neutrino [GeV]
  double ND_E_vis_true; // True visible energy of neutrino [GeV]

  t_effTree->SetBranchAddress("ND_Sim_n_hadronic_Edep_b",         &ND_Sim_n_hadronic_Edep_b);
  // t_effTree->SetBranchAddress("CurrentThrowDepsX",         &CurrentThrowDepsX);
  // t_effTree->SetBranchAddress("CurrentThrowDepsY",         &CurrentThrowDepsY);
  // t_effTree->SetBranchAddress("CurrentThrowDepsZ",         &CurrentThrowDepsZ);
  t_effTree->SetBranchAddress("CurrentThrowTotE",          &CurrentThrowTotE);
  t_effTree->SetBranchAddress("HadronHitEdeps",            &HadronHitEdeps);
  t_effTree->SetBranchAddress("ND_E_vis_true",             &ND_E_vis_true);
  t_effTree->SetBranchAddress("ND_Gen_numu_E",             &ND_Gen_numu_E);
  t_effValues->SetBranchAddress("totEnergyFDatND_f",   &totEnergyFDatND_f);
  // t_effTree->SetBranchAddress("ND_OffAxis_Sim_hadronic_hit_xyz",            &ND_OffAxis_Sim_hadronic_hit_xyz);


  // Read PosVec
  TChain *t_PosVec = new TChain("effPosND");
  t_PosVec->Add(FileIn.Data());
  vector<Double_t> *ND_LAr_dtctr_pos_vec = 0; // unit: cm, ND off-axis choices for each FD evt
  vector<Double_t> *ND_vtx_vx_vec = 0;       // unit: cm, vtx x choices for each FD evt in ND volume
  vector<Int_t> *iwritten_vec = 0;
  TBranch *b_ND_LAr_dtctr_pos_vec = 0;
  TBranch *b_ND_vtx_vx_vec = 0;
  TBranch *b_iwritten_vec = 0;
  t_PosVec->SetBranchAddress("ND_LAr_dtctr_pos_vec",      &ND_LAr_dtctr_pos_vec , &b_ND_LAr_dtctr_pos_vec);
  t_PosVec->SetBranchAddress("ND_vtx_vx_vec",             &ND_vtx_vx_vec,         &b_ND_vtx_vx_vec);
  t_PosVec->SetBranchAddress("iwritten_vec",              &iwritten_vec,          &b_iwritten_vec);

  Long64_t tentry = t_PosVec->LoadTree(0);
  b_ND_LAr_dtctr_pos_vec->GetEntry(tentry);
  b_ND_vtx_vx_vec->GetEntry(tentry);
  b_iwritten_vec->GetEntry(tentry);

  Int_t ND_LAr_dtctr_pos_vec_size = ND_LAr_dtctr_pos_vec->size();
  Int_t ND_vtx_vx_vec_size = ND_vtx_vx_vec->size();
  Int_t iwritten_vec_size = iwritten_vec->size();
  Int_t tot_size = ND_LAr_dtctr_pos_vec_size*ND_vtx_vx_vec_size;
  Int_t hadronhit_n_plots = tot_size*N_throws;
  Int_t nentries = t_effValues->GetEntries();

  cout<<" tot size: "<<tot_size<<endl;
  Int_t nentrieseffTree = t_effTree->GetEntries();


  for(Int_t i = 0; i< nentrieseffTree; i++){
    t_effValues->GetEntry(i);

    if (i==0){
      int count = 0;
      cout<<" Entry: "<<i<<endl;
      for(float valueVeto : *VetoEnergyEventsPass){
        //cout<<" veto E: "<<valueVeto<<endl;
        if(valueVeto > 30.00)
          cout<<" more than 30 not correct!! "<<endl;
      }
      for(uint64_t value : *ThrowResults){
        count+=1;
      }
    }
  }
  // vector<Int_t> a_ND_off_axis_pos_vec = {-2800, -1600, 0};
  // vector<Double_t> a_ND_off_axis_pos_vec = {0, -10, -25};
  vector<Double_t> a_ND_off_axis_pos_vec = {0, -1.75, -2, -4, 5.75, -8, -9.75, -12, -13.75, -16, -17.75, -20, -21.75, -24, -25.75, -26.25, -28, -28.25};
 // vector<Double_t> a_ND_off_axis_pos_vec = {4,2, 0, -2, -4, -6, -8, -10, -12, -14, -16, -18, -20, -22, -24, -26, -28, -30, -32};
  vector<Double_t> a_ND_vtx_vx_vec;
  vector<Double_t> OAPosition;
  // Generate six evenly spaced point in each non-dead region
  vector<pair<double, double>> non_dead_regions = {
      {-300, -256.55}, {-253.95, -204.95}, {-203.45, -154.45}, {-151.85, -102.85},
      {-101.35, -52.35}, {-49.75, -0.75}, {0.75, 49.75}, {52.35, 101.35},
      {102.85, 151.85}, {154.45, 203.45}, {204.95, 253.95}, {256.55, 300}
  };

  // Iterate over each non-dead region and generate points
  for (const auto& region : non_dead_regions) {
      auto points = generatePoints(region.first);
      // std::cout << "Points between " << region.first << " and " << region.second << ": ";
      for (auto point : points) {
          std::cout << point << " ";
          a_ND_vtx_vx_vec.emplace_back(point);
      }
      std::cout << std::endl;
  }

  //mockup Coefficients (based on 1 histo..)
  vector<Double_t> Coefficients;

  // Define the non-dead regions
  // if (true)
  // for (auto x : a_ND_vtx_vx_vec){
  //     // std::cout <<" vtx: "<< x << ",  ";
  //     if (x<-250)
  //       Coefficients.emplace_back(69*1e-6);
  //     if(x>=-250 && x<-200)
  //       Coefficients.emplace_back(100*1e-6);
  //     if(x>=-200 && x<-150)
  //       Coefficients.emplace_back(20*1e-6);
  //     if(x>=-150 && x<-100)
  //       Coefficients.emplace_back(25*1e-6);
  //     if(x>=-100 && x<-50)
  //       Coefficients.emplace_back(28*1e-6);
  //     if(x>=-50 && x<0)
  //       Coefficients.emplace_back(29*1e-6);
  //     if(x>=0 && x< 50)
  //       Coefficients.emplace_back(32*1e-6);
  //     if(x>=50 && x<100)
  //       Coefficients.emplace_back(27*1e-6);
  //     if(x>=100 && x<150)
  //       Coefficients.emplace_back(19*1e-6);
  //     if(x>=150 && x<=300)
  //       Coefficients.emplace_back(10*1e-6);
  // }

  // TGraph* GraphCoeffVsVtxX = new TGraph(a_ND_vtx_vx_vec.size(), a_ND_vtx_vx_vec.data(),Coefficients.data());

  // Set Palette
  gStyle->SetPalette(55);
  //  gStyle->SetOptStat(0);


  int nFDEvents = iwritten_vec->size();
  int nvtxXpositions = a_ND_vtx_vx_vec.size();
  int nDetPos = a_ND_off_axis_pos_vec.size();

  double OAPos;
  double OAPos2;
  double CoefficientsAtOAPos;
  double WeightEventsAtOaPos;
  cout<<" nFDEvents = "<<nFDEvents<<" nXpos = "<<nvtxXpositions<<endl;

  TH1D* HistVetoE[nFDEvents][nvtxXpositions];
  TH1D* HistEtrim[nFDEvents][nvtxXpositions];
  TH1D* HistEtrimDetPos[nFDEvents][nvtxXpositions][nDetPos];
  TH1D* HistEtrimDetPosCoeff1[nFDEvents][nvtxXpositions][nDetPos];
  TH1D* HistEtrimAllVtxX[nFDEvents];
  TH1D* HistEtrimAllVtxXTimesCoeff[nFDEvents];
  TH1D* HistOAPosVtxX[nFDEvents][nvtxXpositions];




  // TCanvas** CanvasTrimE = new TCanvas*[nvtxXpositions];

  // save histo with total hadronic energy at FD for all FD events
  TH1D* hist_FDTotEnergy = new TH1D("hist_FDTotEnergy", "hist_FDTotEnergy", 25000, 0, 25000);
  //save histo with total neutrino energy at FD for all FD events
  TH1D* hist_EnuFDEnergy = new TH1D("hist_EnuFDEnergy", "hist_EnuFDEnergy", 25000, 0, 25);
  // save histo with visible neutrino energy at FD for all FD events
  TH1D* hist_visEnuFDEnergy = new TH1D("hist_visEnuFDEnergy", "hist_visEnuFDEnergy", 25000, 0, 25);

  TGraph* PlotEfficiencyVsVtxX[nFDEvents];
  TGraph* EfficiencyVsOAPos[nFDEvents];
  // TH1D* HistOAPos = new TH1D("HistOAPos", "HistOAPos", 67, -30.5, 3);
  TH1D* HistOAPos[nFDEvents];
  TH2D* CoefficientsAtOAPosHist = new TH2D("CoefficientsAtOAPosHist", "CoefficientsAtOAPosHist", 67, -30.5, 3, 60, -0.3, 0.3);


  TFile* FileWithHistoInfo = new TFile("FileWithHist_10EventsFromJobSub_NoTrimX_CoeffsAndOAPoswithSameEff.root", "RECREATE");

  //===want to look only at Event3 with Eff =1 at all Vtx_x for now ====

  //first get the weights/distribution of events at OA postions -> later on weight the histograms to these POT-like histo
  for (Int_t i_iwritten = 0; i_iwritten<10; i_iwritten++)
  { HistOAPos[i_iwritten] = new TH1D(Form("HistOAPos_FDEvt_%d", i_iwritten), Form("HistOAPos_FDEvt_%d", i_iwritten), 67, -30.5, 3);
      for (Double_t i_ND_LAr_vtx_pos: a_ND_vtx_vx_vec)
      {
        Int_t i_entry = tot_size * i_iwritten;
        for (i_entry ; i_entry < tot_size * (i_iwritten+1); i_entry++ )
        {
          t_effTree->GetEntry(i_entry);
          t_effValues->GetEntry(i_entry);

          //cout<<" i_iwritten "<<i_iwritten<<" totEnergyFDatND_f " <<totEnergyFDatND_f<<endl;

          if ( ND_LAr_vtx_pos == i_ND_LAr_vtx_pos ){
            if(i_ND_LAr_vtx_pos == -298.55) {//only want to fill the histos for 1 vtxX
              hist_FDTotEnergy->Fill(totEnergyFDatND_f);
              hist_EnuFDEnergy->Fill(ND_Gen_numu_E);
              hist_visEnuFDEnergy->Fill(ND_E_vis_true);
            }

            Int_t i_detposInCoefRange = 0;
            for (Double_t i_ND_LAr_dtctr_pos: a_ND_off_axis_pos_vec)
            {
              i_detposInCoefRange +=1;
              // //calculate OAPos=vtx_x+det_pos
              OAPos2 = ND_LAr_vtx_pos/100.0 + a_ND_off_axis_pos_vec[i_detposInCoefRange-1];
              HistOAPos[i_iwritten]->Fill(OAPos2);
            }
          }
        }
      }
      //cout<<" first loop iwritten: "<<i_iwritten<<endl;
  }

 // start the loop with efficiencies etc

  for (Int_t i_iwritten = 0; i_iwritten<nFDEvents; i_iwritten++)
  {
    cout<<" i_iwritten: "<<i_iwritten<<endl;
    Double_t x_ND_LAr_vtx_pos[ND_vtx_vx_vec_size];
    Double_t y_geoeff[ND_vtx_vx_vec_size];
    Double_t y_geoEffOAPos[ND_vtx_vx_vec_size * nDetPos];
    Double_t X_OAPos[ND_vtx_vx_vec_size * nDetPos];

    TString HistEtrimAllVtxX_name = Form("HistEtrimAllVtxX_FDEvt_%d", i_iwritten);
    HistEtrimAllVtxX[i_iwritten] = new TH1D(HistEtrimAllVtxX_name, HistEtrimAllVtxX_name, 25000, 0, 25000);
    TString HistEtrimAllVtxXTimesCoeff_name = Form("HistEtrimAllVtxXTimesCoeff_FDEvt_%d", i_iwritten);
    HistEtrimAllVtxXTimesCoeff[i_iwritten] = new TH1D(HistEtrimAllVtxXTimesCoeff_name, HistEtrimAllVtxXTimesCoeff_name, 25000, 0, 25000);



    Int_t n_plot = 0;
    Int_t i_n_plot = 0;
    Int_t i_vtxX_plot=0;
    Int_t iOAPostAtVtxX = 0;





      for (Double_t i_ND_LAr_vtx_pos: a_ND_vtx_vx_vec)
      {

        i_vtxX_plot +=1;

        TString HistVetoE_name = Form("HistVetoE_FDEvt_%d_vtxXpost_%f", i_iwritten, i_ND_LAr_vtx_pos);
        HistVetoE[i_iwritten][i_vtxX_plot-1] = new TH1D(HistVetoE_name, HistVetoE_name, 60, 0, 60);
        TString HistEtrim_name = Form("HistEtrim_FDEvt_%d_vtxXpost_%f", i_iwritten, i_ND_LAr_vtx_pos);
        HistEtrim[i_iwritten][i_vtxX_plot-1] = new TH1D(HistEtrim_name, HistEtrim_name, 25000, 0, 25000);

        TString HistOAPosVtxX_name = Form("HistOAPosVtxX_FDEvt_%d_vtxXpost_%f", i_iwritten, i_ND_LAr_vtx_pos);
        HistOAPosVtxX[i_iwritten][i_vtxX_plot-1] = new TH1D(HistOAPosVtxX_name, HistOAPosVtxX_name, 72, -300, 300);

        // cout<<" coeff: "<<Coefficients[i_vtxX_plot-1]<<endl;
        Int_t i_detpos_Hist=0;

        Int_t i_entry = tot_size * i_iwritten;
        //cout<<" i entry: "<<i_entry<<endl;
        for (i_entry ; i_entry < tot_size * (i_iwritten+1); i_entry++ )
        {
          t_effTree->GetEntry(i_entry);
          t_effValues->GetEntry(i_entry);

          //cout<<" i_iwritten "<<i_iwritten<<" totEnergyFDatND_f " <<totEnergyFDatND_f<<endl;

          if ( ND_LAr_vtx_pos == i_ND_LAr_vtx_pos ){

            x_ND_LAr_vtx_pos[i_vtxX_plot-1] = ND_LAr_vtx_pos/100;
            y_geoeff[i_vtxX_plot-1] = ND_GeoEff;
            PlotEfficiencyVsVtxX[i_iwritten] = new TGraph(a_ND_vtx_vx_vec.size(), x_ND_LAr_vtx_pos, y_geoeff);


            // cout<<"here??"<<endl;
            // Double_t offset_X = i_ND_LAr_dtctr_pos;
            // Store hadron hit info into vector
            vector<float> ThrowDepsX_hit;
            vector<vector<float>> ThrowDepsX;
            ThrowDepsX.clear();
            int it_throw_x_counter =0;
            double averageTotalE ; //totalE at vtx_x = 7.75

            int throw_couter = 0;

            int nthrowsToLoop = NPassedThrows; //this is going to be the validThrows


            for (Int_t ithrow = 0; ithrow < nthrowsToLoop; ithrow++ ){
                // cout << "ithrow: " << ithrow <<endl;
                // for all events
                // n_plot = i_n_plot*N_throws + ithrow;

                // calculate energy outside the ND active volume
                // Double_t CurrentThrowoutEnergyND = 0.;
                // Double_t CurrentThrowvetoEnergyND = 0.;
                //
                // CurrentThrowvetoEnergyND = VetoEnergyEventsPass->at(ithrow);
                //
                // HistVetoE[i_iwritten][i_vtxX_plot-1]->Fill(CurrentThrowvetoEnergyND);
                  // cout<<" tot energy: "<<CurrentThrowTotE->at(ithrow)<<endl;
                HistEtrim[i_iwritten][i_vtxX_plot-1]->Fill(TrimEnergyEventsPass->at(ithrow));
                // cout<<" fil hit : "<<TrimEnergyEventsPass->at(ithrow)<<endl;
                HistEtrimAllVtxX[i_iwritten]->Fill(TrimEnergyEventsPass->at(ithrow), ND_GeoEff/(nthrowsToLoop*a_ND_vtx_vx_vec.size()));


              } //end throw
              // cout<<" end of throw loop"<<endl;

              PlotEfficiencyVsVtxX[i_iwritten]->SetTitle(Form("TotalFD Energy = %.2f MeV", totEnergyFDatND_f));


              HistEtrim[i_iwritten][i_vtxX_plot-1]->Scale(ND_GeoEff/nthrowsToLoop);
              HistEtrim[i_iwritten][i_vtxX_plot-1]->GetYaxis()->SetTitle("norm");
              HistEtrim[i_iwritten][i_vtxX_plot-1]->GetXaxis()->SetTitle("FD E_{trim} (MeV)");
              // HistEtrim[i_iwritten][i_vtxX_plot-1]->Draw("hist");
              //here loop over DetPos here -> for now just assume same efficiency at every det position..
              Int_t i_detpos = 0;


              for (Double_t i_ND_LAr_dtctr_pos: a_ND_off_axis_pos_vec)
              {
                i_detpos+=1;
                iOAPostAtVtxX += 1;

                // //calculate OAPos=vtx_x+det_pos
                OAPos = ND_LAr_vtx_pos/100.0 + a_ND_off_axis_pos_vec[i_detpos-1];
                y_geoEffOAPos[iOAPostAtVtxX-1] = ND_GeoEff;
                X_OAPos[iOAPostAtVtxX-1] = OAPos;
                EfficiencyVsOAPos[i_iwritten] =  new TGraph(ND_vtx_vx_vec_size * nDetPos, X_OAPos, y_geoEffOAPos);



                TString HistEtrimDetPos_name = Form("HistEtrim_FDEvt_%d_vtxXpost_%f_DetPos_%f", i_iwritten, i_ND_LAr_vtx_pos, a_ND_off_axis_pos_vec[i_detpos-1] );
                //same efficiency at all detecto positions means same events passing the cuts so same HisEtrim[i_vtxX_plot-1]


                HistEtrimDetPos[i_iwritten][i_vtxX_plot-1][i_detpos-1] = (TH1D*) HistEtrim[i_iwritten][i_vtxX_plot-1]->Clone();

                HistEtrimDetPos[i_iwritten][i_vtxX_plot-1][i_detpos-1]->SetName(HistEtrimDetPos_name);

                // TString HistEtrimDetPosCoeff1_name = Form("HistEtrimCoeff1_FDEvt_%d_vtxXpost_%f_DetPos_%f", i_iwritten, i_ND_LAr_vtx_pos, a_ND_off_axis_pos_vec[i_detpos-1] );
                // HistEtrimDetPosCoeff1[i_iwritten][i_vtxX_plot-1][i_detpos-1] = (TH1D*) HistEtrim[i_iwritten][i_vtxX_plot-1]->Clone();
                // HistEtrimDetPosCoeff1[i_iwritten][i_vtxX_plot-1][i_detpos-1]->SetName(HistEtrimDetPosCoeff1_name);



                CoefficientsAtOAPos = CoefficientsHist->GetBinContent(CoefficientsHist->FindBin(OAPos));
                CoefficientsAtOAPosHist->Fill(OAPos, CoefficientsAtOAPos);
                WeightEventsAtOaPos = HistOAPos[i_iwritten]->GetBinContent(HistOAPos[i_iwritten]->FindBin(OAPos));

                //multiply histo by coefficients value and weight to the nr of entries (i.e how many Etrim histograms in a given OA pos with 0.5 cm width)
                HistEtrimDetPos[i_iwritten][i_vtxX_plot-1][i_detpos-1]->Scale(CoefficientsAtOAPos * 1.0/WeightEventsAtOaPos);
                HistEtrimDetPos[i_iwritten][i_vtxX_plot-1][i_detpos-1]->Write(HistEtrimDetPos[i_iwritten][i_vtxX_plot-1][i_detpos-1]->GetName());





              }//end LAr pos

              HistEtrim[i_iwritten][i_vtxX_plot-1]->Write(HistEtrim[i_iwritten][i_vtxX_plot-1]->GetName());
              delete HistEtrim[i_iwritten][i_vtxX_plot-1];

            //  cout<<" integral: "<<HistEtrim[i_iwritten][i_vtxX_plot-1]->Integral()<< " nd geo: "<<ND_GeoEff;
          }// end vtx selection
        }//end ientry

      }//end vtx pos inside LAr




    PlotEfficiencyVsVtxX[i_iwritten]->SetMarkerStyle(20);
    // PlotEfficiencyVsVtxX[i_iwritten]->Draw("AP");
    PlotEfficiencyVsVtxX[i_iwritten]->Write(Form("GraphEffficiency_FDEventNr_%d",i_iwritten));
    EfficiencyVsOAPos[i_iwritten]->Write(Form("GraphEffficiency_AllNDDetPos_FDEventNr_%d",i_iwritten));
    HistEtrimAllVtxX[i_iwritten]->Write(HistEtrimAllVtxX[i_iwritten]->GetName());

    HistEtrimAllVtxXTimesCoeff[i_iwritten] = (TH1D*)HistEtrimDetPos[i_iwritten][0][0]->Clone();
    HistEtrimAllVtxXTimesCoeff[i_iwritten]->Reset();
    HistEtrimAllVtxXTimesCoeff[i_iwritten]->SetName(HistEtrimAllVtxXTimesCoeff_name);
    HistEtrimAllVtxXTimesCoeff[i_iwritten]->SetTitle(Form("TotalFD Energy = %.2f MeV", totEnergyFDatND_f));

    // HistEtrimAllVtxXTimesCoeffEq1[i_iwritten] = (TH1D*)HistEtrimDetPosCoeff1[i_iwritten][0][0]->Clone();
    // HistEtrimAllVtxXTimesCoeffEq1[i_iwritten]->Reset();
    // HistEtrimAllVtxXTimesCoeffEq1[i_iwritten]->SetName("HistEtrimAllVtxXTimesCoeffEq1");
    //scale histogram with coeffs applied
    // HistEtrimAllVtxXTimesCoeff[i_iwritten]->Scale(ND_GeoEff/(nthrowsToLoop*a_ND_vtx_vx_vec.size());
    for(int ivtxX = 0; ivtxX < nvtxXpositions; ivtxX++){
      for(int iDetPos = 0; iDetPos < nDetPos; iDetPos++){
        if(HistEtrimDetPos[i_iwritten][ivtxX][iDetPos]->GetEntries() > 0){
          HistEtrimAllVtxXTimesCoeff[i_iwritten]->Add(HistEtrimDetPos[i_iwritten][ivtxX][iDetPos]);
        }
        // if(HistEtrimDetPosCoeff1[i_iwritten][ivtxX][iDetPos]->GetEntries() > 0)
        //   HistEtrimAllVtxXTimesCoeffEq1[i_iwritten]->Add(HistEtrimDetPosCoeff1[i_iwritten][ivtxX][iDetPos]);

          //cout<<" name "<<HistEtrimAllVtxXTimesCoeff[i_iwritten]->GetName()<<endl;
          // //add up the other histos at vtx_x and DetPos
          // if(ivtxX!=0 && iDetPos!=0)

            //cout<<"ivtxx: "<<ivtxX<< " iDetPos "<< iDetPos<<" mean: "<<HistEtrimDetPos[i_iwritten][ivtxX][iDetPos]->GetMean()<<endl;
            delete HistEtrimDetPos[i_iwritten][ivtxX][iDetPos];

      }
    }

    // HistEtrimAllVtxXTimesCoeff[i_iwritten]->Scale(1.0/(nDetPos));
    cout<<" total DetPos: "<<nDetPos<<endl;

    HistEtrimAllVtxXTimesCoeff[i_iwritten]->Write(HistEtrimAllVtxXTimesCoeff[i_iwritten]->GetName());

    // HistEtrimAllVtxXTimesCoeffEq1[i_iwritten]->Write("HistEtrimAllVtxXTimesCoeffEq1");



    // cout<<" ndet pos = "<<nDetPos<<endl;

    HistOAPos[i_iwritten]->Write(Form("HistOAPos_FDEvt_%d", i_iwritten));

  }//end iwritten



  CoefficientsAtOAPosHist->Write("CoefficientsAtOAPosHist");
  CoefficientsHist->Write("CoefficientsHist");
//  hist_FDTotEnergy->Scale(1.0/72); // we have 72 vtxX positions therefore instead of 1 event in the histo we have 72 entries
  hist_FDTotEnergy->Write("hist_FDTotEnergy");
  hist_EnuFDEnergy->Write("hist_EnuFDEnergy");
  hist_visEnuFDEnergy->Write("hist_visEnuFDEnergy");

  // HistOAPos->Draw("hist");

  //FileWithHistoInfo->Write();
  FileWithHistoInfo->Close();

  delete hist_FDTotEnergy;
  delete hist_EnuFDEnergy;
  delete hist_visEnuFDEnergy;

  cout<<" before delete "<<endl;


  // delete all canvas





}
