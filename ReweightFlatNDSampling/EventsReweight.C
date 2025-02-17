void EventsReweight(){

  TFile* file0m = TFile::Open("FileWithHistos_andTree_0m_280kA.root");
  TFile* file2m = TFile::Open("FileWithHistos_andTree_2m.root");

  TTree *tree0m = (TTree*)file0m->Get("treeVars");
  TTree *tree2m = (TTree*)file2m->Get("treeVars");

  double Ev;
  double vtx_x;
  double abspos_x;
  double det_x;
  bool SignalLikeOnlyNuMuCC;

  vector<double> EvAll;
  vector<double> vtx_xAll;
  vector<double> abspos_All;

  tree0m->SetBranchAddress("Ev", &Ev);
  tree0m->SetBranchAddress("vtx_x", &vtx_x);
  tree0m->SetBranchAddress("abspos_x", &abspos_x);
  tree0m->SetBranchAddress("SignalLikeOnlyNuMuCC", &SignalLikeOnlyNuMuCC);

  Long64_t nEntries0m = tree0m->GetEntries();
  for(int i = 0; i<nEntries0m; i++){
    tree0m->GetEntry(i);
    if(SignalLikeOnlyNuMuCC){
      EvAll.push_back(Ev);
      vtx_xAll.push_back(vtx_x);
      abspos_All.push_back(abspos_x);
    }
  }

  tree2m->SetBranchAddress("Ev", &Ev);
  tree2m->SetBranchAddress("vtx_x", &vtx_x);
  tree2m->SetBranchAddress("abspos_x", &abspos_x);
  tree2m->SetBranchAddress("SignalLikeOnlyNuMuCC", &SignalLikeOnlyNuMuCC);

  Long64_t nEntries2m = tree2m->GetEntries();
  for(int i = 0; i<nEntries2m; i++){
    tree2m->GetEntry(i);
    if(SignalLikeOnlyNuMuCC){
      EvAll.push_back(Ev);
      vtx_xAll.push_back(vtx_x);
      abspos_All.push_back(abspos_x);
    }
  }

  const int nDetSlices = 12;
  TH2D* EventsVsDetSliceVsOApos = new TH2D("EventsVsDetSliceVsOApos", "EventsVsDetSliceVsOApos",  64, -2, 30,12 , -3, 3);
  for(int i = 0; i < vtx_xAll.size(); i++){
    EventsVsDetSliceVsOApos->Fill( abspos_All[i]/100, vtx_xAll[i]/100);
    // cout<<vtx_xAll[i]<<endl;
  }

  EventsVsDetSliceVsOApos->Draw("COLZ");

  // flux bins (for now by hand but will make it from function TODO)
  vector<double> LowBinEdges = {-2, -1.5, -1, -0.5, 0, 0, 0, 0, 0};
  vector<double> HighBinEdges = {2, 2, 2, 2, 2.5, 3, 3.5, 4, 2};

  int nFluxBins = LowBinEdges.size();
  TH2D* EventsVsDetSliceVsOAposAtFluxBin[nFluxBins];
  //Tuee Enu energy binning
  std::vector<double> edges;
  edges.emplace_back(0.);
  // edges.emplace_back(0.15);
  // edges.emplace_back(0.25);
  // edges.emplace_back(0.35);
  // edges.emplace_back(0.45);
  edges.emplace_back(0.5);
  while (edges.back() < 2.) {
    edges.emplace_back(edges.back() + 0.04);
  }
  while (edges.back() < 3.) {
    edges.emplace_back(edges.back() + 0.08);
  }
  while (edges.back() < 4.) {
    edges.emplace_back(edges.back() + 0.1);
  }
  edges.emplace_back(4.5);
  edges.emplace_back(5.);
  edges.emplace_back(6.);
  edges.emplace_back(10.);
  edges.emplace_back(120.);

  int N = edges.size() -1;
  Double_t edges2[N+1];
  for(int i = 0; i<N+1; i++)
    edges2[i] = edges[i];

  TH1D* AverageEnu[nFluxBins];
  TH1D* HistTrueEnufluxBinDetSlice[nFluxBins][nDetSlices];

  for(int ifluxBin = 0; ifluxBin<nFluxBins; ifluxBin++){
    for(int binIndexDetSlice = 0; binIndexDetSlice<nDetSlices; binIndexDetSlice++){
      HistTrueEnufluxBinDetSlice[ifluxBin][binIndexDetSlice] = new TH1D(Form("EnuFlux_fluxBin_%d_detSlice_%d", ifluxBin, binIndexDetSlice),Form("EnuFlux_fluxBin_%d_detSlice_%d", ifluxBin, binIndexDetSlice), N, edges2 );
    }
  }

  for(int ifluxBin = 0; ifluxBin<nFluxBins; ifluxBin++){
    TString flux_name = Form("fluxhisto_%d", ifluxBin);
    EventsVsDetSliceVsOAposAtFluxBin[ifluxBin] = new TH2D(flux_name, flux_name, 64, -2, 30,12 , -3, 3);
    TString Enu_name = Form("Enuatfluxbin_%d", ifluxBin);
    AverageEnu[ifluxBin] = new TH1D(Enu_name, Enu_name, N, edges2 );
    for(int i = 0; i<EvAll.size(); i++){
      if(abspos_All[i]/100>=LowBinEdges[ifluxBin] && abspos_All[i]/100 <= HighBinEdges[ifluxBin]){
        EventsVsDetSliceVsOAposAtFluxBin[ifluxBin]->Fill( abspos_All[i]/100, vtx_xAll[i]/100);
        AverageEnu[ifluxBin]->Fill(EvAll[i]);
        int binIndexDetSlice = EventsVsDetSliceVsOApos->GetYaxis()->FindBin(vtx_xAll[i] / 100) - 1;
        // cout<<" vtx: "<<vtx_xAll[i] / 100<<" bin index: "<<binIndexDetSlice<<endl;
        HistTrueEnufluxBinDetSlice[ifluxBin][binIndexDetSlice]->Fill(EvAll[i]);
      }
    }

  }


  TCanvas* CanvasFluxbinsHisto = new TCanvas("CanvasFluxbinsHisto", "CanvasFluxbinsHisto");
  CanvasFluxbinsHisto->Divide(3,3);
  for(int ifluxBin = 0; ifluxBin<nFluxBins; ifluxBin++){
    CanvasFluxbinsHisto->cd(ifluxBin+1);
    EventsVsDetSliceVsOAposAtFluxBin[ifluxBin]->SetTitle(Form("Flux bin %d from %.2f to %.2f m ", ifluxBin,LowBinEdges[ifluxBin], HighBinEdges[ifluxBin] ));
    EventsVsDetSliceVsOAposAtFluxBin[ifluxBin]->GetXaxis()->SetRangeUser(-2,5);
    EventsVsDetSliceVsOAposAtFluxBin[ifluxBin]->Draw("COLZ");
  }

    TCanvas* CanvasAverageEAtFluxbinsHisto = new TCanvas("CanvasAverageEAtFluxbinsHisto", "CanvasAverageEAtFluxbinsHisto");
    CanvasAverageEAtFluxbinsHisto->Divide(3,3);

    for(int ifluxBin = 0; ifluxBin<nFluxBins; ifluxBin++){
      CanvasAverageEAtFluxbinsHisto->cd(ifluxBin+1);
      AverageEnu[ifluxBin]->SetTitle(Form("Average Enu at Flux bin %d ", ifluxBin ));
      AverageEnu[ifluxBin]->GetXaxis()->SetRangeUser(0,10);

      AverageEnu[ifluxBin]->Scale(1, "width");
      AverageEnu[ifluxBin]->Scale(1.0/(nDetSlices-4));
      AverageEnu[ifluxBin]->GetYaxis()->SetRangeUser(0, 240E3);
      AverageEnu[ifluxBin]->SetLineWidth(2);
      AverageEnu[ifluxBin]->SetLineColor(kOrange+1);
      AverageEnu[ifluxBin]->Draw("hist");
    }

    TCanvas* CanvasEventsAtFluxBinDiffDetSlices = new TCanvas("CanvasEventsAtFluxBinDiffDetSlices", "CanvasEventsAtFluxBinDiffDetSlices");
    CanvasEventsAtFluxBinDiffDetSlices->Divide(4,2);
    // will plot them 1 by 1 in terms of flux bins for now -- cross check everything is ok
    int ifluxBinPlot = 0;
    for(int i = 2; i<nDetSlices-2; i++){ //first and last 2 detector bins are always empty
      CanvasEventsAtFluxBinDiffDetSlices->cd(i+1-2);
      HistTrueEnufluxBinDetSlice[ifluxBinPlot][i]->SetLineWidth(2);
      HistTrueEnufluxBinDetSlice[ifluxBinPlot][i]->GetXaxis()->SetRangeUser(0,10);
      HistTrueEnufluxBinDetSlice[ifluxBinPlot][i]->Scale(1, "width");
      HistTrueEnufluxBinDetSlice[ifluxBinPlot][i]->GetYaxis()->SetRangeUser(0,300E3);
      HistTrueEnufluxBinDetSlice[ifluxBinPlot][i]->Draw("hist");
      AverageEnu[ifluxBinPlot]->SetLineColor(kOrange+1);
      AverageEnu[ifluxBinPlot]->Draw("histsames");
    }

    // calculate the weights
    TH1D* HistWeights[nFluxBins][nDetSlices];
    for(int iFluxBin = 0; iFluxBin < nFluxBins; iFluxBin++){
      for(int idetSlice = 2; idetSlice < nDetSlices-2; idetSlice++){
        HistWeights[iFluxBin][idetSlice] = (TH1D*) AverageEnu[iFluxBin]->Clone();
        HistWeights[iFluxBin][idetSlice]->Divide(HistTrueEnufluxBinDetSlice[iFluxBin][idetSlice]);
        HistWeights[iFluxBin][idetSlice]->SetName(Form("EnuFlux_fluxBin_%d_detSlice_%d", iFluxBin, idetSlice));
        HistWeights[iFluxBin][idetSlice]->SetTitle(Form("EnuFlux_fluxBin_%d_detSlice_%d", iFluxBin, idetSlice));
      }
    }

    TCanvas* CanvasWeightsAtFluxBinDiffDetSlices = new TCanvas("CanvasWeightsAtFluxBinDiffDetSlices", "CanvasWeightsAtFluxBinDiffDetSlices");
    CanvasWeightsAtFluxBinDiffDetSlices->Divide(4,2);
    // will plot them 1 by 1 in terms of flux bins for now -- cross check everything is ok
    for(int i = 2; i<nDetSlices-2; i++){ //first and last 2 detector bins are always empty
      CanvasWeightsAtFluxBinDiffDetSlices->cd(i+1-2);
      HistWeights[ifluxBinPlot][i]->SetLineWidth(2);
      HistWeights[ifluxBinPlot][i]->SetLineColor(2);
      HistWeights[ifluxBinPlot][i]->GetXaxis()->SetRangeUser(0,10);
      HistWeights[ifluxBinPlot][i]->GetYaxis()->SetRangeUser(0, 1.5);
      HistWeights[ifluxBinPlot][i]->Draw("hist");
    }







  cout<<" size Ev: "<<EvAll.size()<<endl;

}
