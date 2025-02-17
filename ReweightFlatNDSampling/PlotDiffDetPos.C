void PlotDiffDetPos(){

  TFile* file0m = TFile::Open("FileWithHistos_0m_280kA.root");
  TFile* file1p75m = TFile::Open("FileWithHistos_1p75m.root");
  TFile* file2m = TFile::Open("FileWithHistos_2m.root");

  TH1D* EventsInVtxX_0m = (TH1D*) file0m->Get("EventsInVtxX");
  TH1D* EventsInVtxX_2m = (TH1D*) file2m->Get("EventsInVtxX");
  TH1D* EventsInVtxX_1p75m = (TH1D*) file1p75m->Get("EventsInVtxX");

  TH1D* EventsInVtxXWithVetoCut_0m = (TH1D*) file0m->Get("EventsInVtxXWithVetoCut");
  TH1D* EventsInVtxXWithVetoCut_2m = (TH1D*) file2m->Get("EventsInVtxXWithVetoCut");
  TH1D* EventsInVtxXWithVetoCut_1p75m = (TH1D*) file1p75m->Get("EventsInVtxXWithVetoCut");

  TH2D* EventsInVtxXvsOffAxisslice_0m = (TH2D*) file0m->Get("EventsInVtxXvsOffAxisslice");
  TH2D* EventsInVtxXvsOffAxisslice_2m = (TH2D*) file2m->Get("EventsInVtxXvsOffAxisslice");
  TH2D* EventsInVtxXvsOffAxisslice_1p75m = (TH2D*) file1p75m->Get("EventsInVtxXvsOffAxisslice");

  TH2D* EventsInVtxXvsOffAxissliceWithVetoCut_0m = (TH2D*) file0m->Get("EventsInVtxXvsOffAxissliceWithVetoCut");
  TH2D* EventsInVtxXvsOffAxissliceWithVetoCut_2m = (TH2D*) file2m->Get("EventsInVtxXvsOffAxissliceWithVetoCut");
  TH2D* EventsInVtxXvsOffAxissliceWithVetoCut_1p75m = (TH2D*) file1p75m->Get("EventsInVtxXvsOffAxissliceWithVetoCut");

  //add all together

  TH1D* EventsInVtxX_AllDetPos = (TH1D*) EventsInVtxX_0m->Clone();
  EventsInVtxX_AllDetPos->Add(EventsInVtxX_2m);
  EventsInVtxX_AllDetPos->Add(EventsInVtxX_1p75m);
  TH1D* EventsInVtxXWithVetoCut_AllDetPos = (TH1D*) EventsInVtxXWithVetoCut_0m->Clone();
  EventsInVtxXWithVetoCut_AllDetPos->Add(EventsInVtxXWithVetoCut_2m);
  EventsInVtxXWithVetoCut_AllDetPos->Add(EventsInVtxXWithVetoCut_1p75m);

  TCanvas* CanvasEventsInVtxX = new TCanvas("CanvasEventsInVtxX", "CanvasEventsInVtxX");
  EventsInVtxX_AllDetPos->Draw("hist");
  EventsInVtxXWithVetoCut_AllDetPos->Draw("histsames");
  TLegend* leg = new TLegend();
  leg->AddEntry(EventsInVtxX_AllDetPos, "no cut");
  leg->AddEntry(EventsInVtxX_AllDetPos, "had veto cut");

  TH2D* EventsInVtxXvsOffAxisslice_AllDetPos = (TH2D*) EventsInVtxXvsOffAxisslice_0m->Clone();
  EventsInVtxXvsOffAxisslice_AllDetPos->Add(EventsInVtxXvsOffAxisslice_2m);
  // EventsInVtxXvsOffAxisslice_AllDetPos->Add(EventsInVtxXvsOffAxisslice_1p75m);
  TCanvas* CanvasEvVsVtxXvsDetSlice = new TCanvas("CanvasEvVsVtxXvsDetSlice", "CanvasEvVsVtxXvsDetSlice");
  CanvasEvVsVtxXvsDetSlice->cd();
  EventsInVtxXvsOffAxisslice_AllDetPos->Draw("COLZ");

  TH2D* EventsInVtxXvsOffAxissliceWithVetoCut_AllDetPos = (TH2D*) EventsInVtxXvsOffAxissliceWithVetoCut_0m->Clone();
  EventsInVtxXvsOffAxissliceWithVetoCut_AllDetPos->Add(EventsInVtxXvsOffAxissliceWithVetoCut_2m);
    EventsInVtxXvsOffAxissliceWithVetoCut_AllDetPos->Add(EventsInVtxXvsOffAxissliceWithVetoCut_1p75m);
  TCanvas* CanvasEvVsVtxXvsDetSliceVetoCut = new TCanvas("CanvasEvVsVtxXvsDetSliceVetoCut", "CanvasEvVsVtxXvsDetSliceVetoCut");
  CanvasEvVsVtxXvsDetSliceVetoCut->cd();
  EventsInVtxXvsOffAxissliceWithVetoCut_AllDetPos->Draw("COLZ");

  TCanvas *CanvasProjOASlice = new TCanvas("CanvasProjOASlice", "CanvasProjOASlice");
  CanvasProjOASlice->cd();
  EventsInVtxXvsOffAxissliceWithVetoCut_AllDetPos->ProjectionX()->Draw("hist");

}
