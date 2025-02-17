void PlotAndReweightOnAxis(){

  TFile* file0m = TFile::Open("FileWithHistos_0m_280kA_trueEhistos.root");
  TH1D* EventsInVtxX_0m = (TH1D*) file0m->Get("EventsInVtxX");
  TH1D* EventsInVtxXWithVetoCut_0m = (TH1D*) file0m->Get("EventsInVtxXWithVetoCut");
  TH2D* EventsInVtxXvsOffAxisslice_0m = (TH2D*) file0m->Get("EventsInVtxXvsOffAxisslice");
  TH1D* HistTrueEnuOnAxis = (TH1D*) file0m->Get("HistTrueEnuOnAxis");
  int nDetSlices = 12;
  TH1D* HistTrueEnuDetSlice[nDetSlices];

  TCanvas* CanvasTrueEnuNoWeight = new TCanvas("CanvasTrueEnuNoWeight", "CanvasTrueEnuNoWeight");
  CanvasTrueEnuNoWeight->Divide(3,4);

  int nEntries = 0;

  HistTrueEnuOnAxis->Scale(1, "width");
  HistTrueEnuOnAxis->Scale(1.0 / (nDetSlices-4)); //this will be the average flux

  for(int i = 0; i<nDetSlices; i++){
    TString histName = (Form("EnuFlux_detSlice_%d", i));
    HistTrueEnuDetSlice[i] = (TH1D*) file0m->Get(histName);
    HistTrueEnuDetSlice[i]->SetDirectory(0);
    CanvasTrueEnuNoWeight->cd(i+1);
    HistTrueEnuDetSlice[i]->Scale(1, "width");
    HistTrueEnuDetSlice[i]->Draw("hist");
    nEntries += HistTrueEnuDetSlice[i]->GetEntries();


  }

  cout<<" entries all det slice hist: "<<nEntries<<" entries Etru = "<<HistTrueEnuOnAxis->GetEntries()<<endl;


  HistTrueEnuOnAxis->SetLineColor(kOrange+2);
  HistTrueEnuOnAxis->SetLineWidth(2);
  HistTrueEnuOnAxis->GetXaxis()->SetRangeUser(0,10);
  TCanvas* CanvasTrueEnuAllBinsVsPerBin = new TCanvas("CanvasTrueEnuAllBinsVsPerBin", "CanvasTrueEnuAllBinsVsPerBin");
  CanvasTrueEnuAllBinsVsPerBin->Divide(4,2);
  HistTrueEnuOnAxis->GetYaxis()->SetRangeUser(0, 145E3);
  for(int i = 2; i<nDetSlices-2; i++){
    CanvasTrueEnuAllBinsVsPerBin->cd(i+1-2);

    HistTrueEnuOnAxis->Draw("hist");
    HistTrueEnuDetSlice[i]->SetLineWidth(2);
    HistTrueEnuDetSlice[i]->Draw("histsames");

  }

  TH1D* HistWeights[nDetSlices];
  TCanvas* CanvasWights = new TCanvas("CanvasWights", "CanvasWights");
  CanvasWights->Divide(4,2);
  for(int i = 2; i<nDetSlices-2; i++){
    HistWeights[i] = (TH1D*) HistTrueEnuOnAxis->Clone();
    HistWeights[i]->Divide(HistTrueEnuDetSlice[i]);
    HistWeights[i]->SetLineColor(2);
    HistWeights[i]->SetLineWidth(2);

    CanvasWights->cd(i+1-2);
    HistWeights[i]->GetYaxis()->SetRangeUser(0.8, 1.22);
    HistWeights[i]->Draw("hist");
  }

  TCanvas* EtrueSpectrumCanvas = new TCanvas("EtrueSpectrumCanvas", "EtrueSpectrumCanvas");
  EtrueSpectrumCanvas->cd();

  HistTrueEnuOnAxis->Draw("hist");



}
