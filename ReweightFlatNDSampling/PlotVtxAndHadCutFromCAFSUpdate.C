inline bool IsInNDFV(double pos_x_cm, double pos_y_cm, double pos_z_cm) {
  bool inDeadRegion = false;
  for (int i = -3; i <= 3; ++i) {
    // 0.5cm cathode in the middle of each module, plus 0.5cm buffer
    double cathode_center = i * 102.1;
    if (pos_x_cm > cathode_center - 0.75 && pos_x_cm < cathode_center + 0.75)
      inDeadRegion = true;

    // 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel
    // plane, x2) don't worry about outer boundary because events are only
    // generated in active Ar + insides
    double module_boundary = i * 102.1 + 51.05;
    if (i <= 2 && pos_x_cm > module_boundary - 1.3 &&
        pos_x_cm < module_boundary + 1.3)
          inDeadRegion = true;
  }
  for (int i = 1; i <= 4; ++i) {
    // module boundaries in z are 1.8cm (0.4cm ArCLight plane + 0.5cm module
    // wall, x2) module is 102.1cm wide, but only 101.8cm long due to cathode
    // (0.5cm) being absent in length but ArCLight is 0.1cm thicker than pixel
    // plane so it's 0.3cm difference positions are off-set by 0.6 because I
    // defined 0 to be the upstream edge based on the active volume by
    // inspecting a plot, and aparently missed by 3 mm, but whatever add 8mm =
    // 2 pad buffer due to worse position resolution in spatial dimension z
    // compared to timing direction x so total FV gap will be 1.8 + 2*0.8
    // = 3.4cm
    double module_boundary = i * 101.8 - 0.6;
    if (pos_z_cm > module_boundary - 1.7 && pos_z_cm < module_boundary + 1.7)
      inDeadRegion = true;
  }
  //cout<<" inDeadRegion "<<inDeadRegion<<endl;

  return (abs(pos_x_cm) < 200 && abs(pos_y_cm) < 100 && pos_z_cm > 50 &&
          pos_z_cm < 350 && !inDeadRegion);
}
void PlotVtxAndHadCutFromCAFSUpdate(){

  TFile* CAFFile = new TFile("/localscratch/icaracas/ImprovePRISMPred/RootFiles/CAFv7_280kA_0m_Chained.root", "READ");
  CAFFile->cd();
  TTree* caf = (TTree*) CAFFile->Get("cafTree");
  assert(caf);

  double vtx_x;
  double vtx_y;
  double vtx_z;
  double det_x;

  double abspos_x;

  int nuPDG;
  int nuPDGunosc;
  int isCC;

  int reco_numu;
  int muon_contained;
  int muon_tracker;
  int reco_q;
  double Ehad_veto;

  double Ev;
  double Ev_reco;

  //vars for VisETrue definition
  double VisEtrue;
  double eP;
  double ePip;
  double ePim;
  double ePi0;
  double eOther;
  double LepE;
  double HadE;

  //vars for EvusRecoND
  double VisERecoND;
  double eRecoP;
  double eRecoPip;
  double eRecoPim;
  double eRecoPi0;
  double eRecoOther;
  double HadEvisReco_ND;
  double Elep_reco;


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

  cout<<edges.size()<<endl;



  caf->SetBranchAddress("vtx_x", &vtx_x);
  caf->SetBranchAddress("vtx_y", &vtx_y);
  caf->SetBranchAddress("vtx_z", &vtx_z);
  caf->SetBranchAddress("det_x", &det_x);

  caf->SetBranchAddress("nuPDG", &nuPDG);
  caf->SetBranchAddress("nuPDGunosc", &nuPDGunosc);
  caf->SetBranchAddress("isCC", &isCC);

  caf->SetBranchAddress("reco_numu", &reco_numu);
  caf->SetBranchAddress("muon_contained", &muon_contained);
  caf->SetBranchAddress("muon_tracker", &muon_tracker);
  caf->SetBranchAddress("reco_q", &reco_q);
  caf->SetBranchAddress("Ehad_veto", &Ehad_veto);

  caf->SetBranchAddress("Ev", &Ev);
  caf->SetBranchAddress("Ev_reco", &Ev_reco);
  //vars for VisETrue
  caf->SetBranchAddress("eP", &eP);
  caf->SetBranchAddress("ePip", &ePip);
  caf->SetBranchAddress("ePim", &ePim);
  caf->SetBranchAddress("ePi0", &ePi0);
  caf->SetBranchAddress("eOther", &eOther);
  caf->SetBranchAddress("LepE", &LepE);
  //vars for VisERecoND
  caf->SetBranchAddress("eRecoP", &eRecoP);
  caf->SetBranchAddress("eRecoPip", &eRecoPip);
  caf->SetBranchAddress("eRecoPim", &eRecoPim);
  caf->SetBranchAddress("eRecoPi0", &eRecoPi0);
  caf->SetBranchAddress("eRecoOther", &eRecoOther);
  caf->SetBranchAddress("Elep_reco", &Elep_reco);


  size_t nevs = caf->GetEntries();

  bool kIsNumuCC;
  bool kIsTrueFV;
  bool kPRISMNDSignal_True_numu; // events user for the Coefficients calculation (True numu CC events)
  bool kIsOutOfTheDesert;
  bool kIsAntiNu;

  bool SignalLikeOnlyNuMuCC;
  bool OnlyNuMuRecoAndQ;

  bool CompleteRecoSelectionNu;

  // TH1D* EvDistribALLSimEvents = new TH1D("EvDistribALLSimEvents", "EvDistribALLSimEvents", N, edges2);
  // TH1D* EvDistribMatch = new TH1D("EvDistribMatch", "EvDistribMatch", N, edges2);
  // TH1D* EvDistribNoRecoCut = new TH1D("EvDistribNoRecoCut", "EvDistribNoRecoCut", 100, 0, 100);
  // TH1D* EvDistribNoRecoCutNC = new TH1D("EvDistribNoRecoCutNC", "EvDistribNoRecoCutNC",100, 0, 100);
  // TH1D* EvDistribNoRecoCutWrongLepton = new TH1D("EvDistribNoRecoCutWrongLepton", "EvDistribNoRecoCutWrongLepton",100, 0, 100);
  // TH1D* EvDistribNoRecoCutMinusBkg = new TH1D("EvDistribNoRecoCutMinusBkg", "EvDistribNoRecoCutMinusBkg", 100, 0, 100);
  // TH1D* EvDistribNoRecoCutNoKantiNuCut = new TH1D("EvDistribNoRecoCutNoKantiNuCut","EvDistribNoRecoCutNoKantiNuCut",100, 0, 100);
  // TH1D* EvDistribALLSelectionCuts = new TH1D("EvDistribALLSelectionCuts","EvDistribALLSelectionCuts",N, edges2);
  // TH1D* EvDistribOnlyNumRecoAndQ = new TH1D("EvDistribOnlyNumRecoAndQ", "EvDistribOnlyNumRecoAndQ", N, edges2);
  //
  //
  // TH1D* EvRecoDistribALLSimEvents = new TH1D("EvRecoDistribALLSimEvents", "EvRecoDistribALLSimEvents", N, edges2);
  // TH1D* EvRecoDistribOnlyNumRecoAndQ = new TH1D("EvRecoDistribOnlyNumRecoAndQ", "EvRecoDistribOnlyNumRecoAndQ", N, edges2);
  // TH1D* EvRecoDistribMatch = new TH1D("EvRecoDistribMatch", "EvRecoDistribMatch", N, edges2);
  // TH1D* EvRecoDistribALLSelectionCuts = new TH1D("EvRecoDistribALLSelectionCuts","EvRecoDistribALLSelectionCuts",N, edges2);
  //
  // TH1D* TrueEvisALLSimEvents = new TH1D("TrueEvisALLSimEvents", "TrueEvisALLSimEvents", N, edges2);
  // TH1D* TrueEvisOnlyNumRecoAndQ = new TH1D("TrueEvisOnlyNumRecoAndQ", "TrueEvisOnlyNumRecoAndQ", N, edges2);
  // TH1D* TrueEvisMatch = new TH1D("TrueEvisMatch", "TrueEvisMatch", N, edges2);
  // TH1D* TrueEvisALLSelectionCuts = new TH1D("TrueEvisALLSelectionCuts","TrueEvisALLSelectionCuts",N, edges2);
  //
  // TH1D* RecoEvisALLSimEvents = new TH1D("RecoEvisALLSimEvents", "RecoEvisALLSimEvents", N, edges2);
  // TH1D* RecoEvisOnlyNumRecoAndQ = new TH1D("RecoEvisOnlyNumRecoAndQ", "RecoEvisOnlyNumRecoAndQ", N, edges2);
  // TH1D* RecoEvisMatch = new TH1D("RecoEvisMatch", "RecoEvisMatch", N, edges2);
  // TH1D* RecoEvisALLSelectionCuts = new TH1D("RecoEvisALLSelectionCuts","RecoEvisALLSelectionCuts",N, edges2);
  //
  // TH2D* EvVsEvRecoDistribALLSimEvents = new TH2D("EvVsEvRecoDistribALLSimEvents", "EvVsEvRecoDistribALLSimEvents", N, edges2, N, edges2);
  // TH2D* EvVsEvRecoDistribMatch = new TH2D("EvVsEvRecoDistribMatch", "EvVsEvRecoDistribMatch", N, edges2, N, edges2);
  // TH2D* EvVsEvRecoDistribOnlyNumRecoAndQ = new TH2D("EvVsEvRecoDistribOnlyNumRecoAndQ", "EvVsEvRecoDistribOnlyNumRecoAndQ", N, edges2,N, edges2);
  // TH2D* EvVsEvRecoDistribALLSelectionCuts = new TH2D("EvVsEvRecoDistribALLSelectionCuts", "EvVsEvRecoDistribALLSelectionCuts", N, edges2, N, edges2);
  //
  // TH1D* nuPDGALLSimEvents = new TH1D("nuPDGALLSimEvents", "nuPDGALLSimEvents", 40, -20, 20);
  // TH1D* nuPDGMatch = new TH1D("nuPDGMatch", "nuPDGMatch", 40, -20, 20);
  // TH1D* nuPDGNoRecoCut = new TH1D("nuPDGNoRecoCut", "nuPDGNoRecoCut", 40, -20, 20);
  //
  // TH1D* EfficiencyEvDistribMatch = new TH1D("EfficiencyEvDistribMatch", "EfficiencyEvDistribMatch", N, edges2);
  // TH1D* EfficiencyEvDistribOnlyNumRecoAndQ = new TH1D("EfficiencyEvDistribOnlyNumRecoAndQ", "EfficiencyEvDistribOnlyNumRecoAndQ", N, edges2);
  // TH1D* EfficiencyEvDistribALLSelectionCuts = new TH1D("EfficiencyEvDistribALLSelectionCuts", "EfficiencyEvDistribALLSelectionCuts", N, edges2);

  TH1D* EventsInVtxX = new TH1D("EventsInVtxX", "EventsInVtxX", 600, -300, 300);
  TH1D* EventsInVtxXWithVetoCut = new TH1D("EventsInVtxXWithVetoCut", "EventsInVtxXWithVetoCut", 600, -300, 300);
  TH1D* EventsInVtxY = new TH1D("EventsInVtxY", "EventsInVtxY", 200, -100, 100);
  TH1D* EventsInVtxYWithVetoCut = new TH1D("EventsInVtxYWithVetoCut", "EventsInVtxYWithVetoCut", 200, -100, 100);
  TH1D* EventsInVtxZ = new TH1D("EventsInVtxZ", "EventsInVtxZ", 300, 50, 350);
  TH1D* EventsInVtxZWithVetoCut = new TH1D("EventsInVtxZWithVetoCut", "EventsInVtxZWithVetoCut", 300, 50, 350);
  TH2D* EventsInVtxXvsOffAxisslice = new TH2D("EventsInVtxXvsOffAxisslice", "EventsInVtxXvsOffAxisslice", 68, -200, 3200, 12, -300, 300);
  TH2D* EventsInVtxXvsOffAxissliceWithVetoCut = new TH2D("EventsInVtxXvsOffAxissliceWithVetoCut", "EventsInVtxXvsOffAxissliceWithVetoCut", 68, -200, 3200, 12, -300, 300);

  int nDetSlices = 12;
  TH1D* HistTrueEnuDetPosOffAxis = new TH1D("HistTrueEnuDetPosOffAxis", "HistTrueEnuDetPosOffAxis", N, edges2);
  TH1D* HistTrueEnuDetSlice[nDetSlices];

  for(int binIndexDetSlice = 0; binIndexDetSlice<nDetSlices; binIndexDetSlice++){
    HistTrueEnuDetSlice[binIndexDetSlice] = new TH1D(Form("EnuFlux_detSlice_%d", binIndexDetSlice),Form("EnuFlux_detSlice_%d", binIndexDetSlice), N, edges2 );
  }

  TFile* FileWithHistos = new TFile("FileWithHistos_andTree_0m_280kA.root", "RECREATE");
  FileWithHistos->cd();

  //Ttree with vars of interest
  TTree* treeVars = new TTree("treeVars", "tree for storing bin reweighting vars");
  treeVars->Branch("Ev", &Ev, "Ev/D");
  treeVars->Branch("vtx_x", &vtx_x, "vtx_x/D");
  treeVars->Branch("abspos_x", &abspos_x, "abspos_x/D");
  treeVars->Branch("det_x", &det_x, "det_x/D");
  treeVars->Branch("SignalLikeOnlyNuMuCC", &SignalLikeOnlyNuMuCC, "SignalLikeOnlyNuMuCC/O");
  treeVars->Branch("Ehad_veto", &Ehad_veto, "Ehad_veto/O");


  for (Long64_t n = 0; n < nevs; ++n) {
    caf->GetEntry(n);

    cout<<" ev nr: "<<n<<" vtx_x " << vtx_x<<" det_x "<< det_x<<endl;

    abspos_x = vtx_x + abs(det_x);
    //define VisETrue
    // HadE = eP + ePip + ePim + ePi0 + eOther;
    // VisEtrue = LepE + HadE;
    //
    // //define VisERecoND
    // HadEvisReco_ND = eRecoP + eRecoPip + eRecoPim + eRecoPi0 + eRecoOther;
    // VisERecoND = HadEvisReco_ND + Elep_reco;

    if(isCC == 1 && abs(nuPDGunosc) == 14 && abs(nuPDG)==14)
      kIsNumuCC=true;
    else
      kIsNumuCC=false;

    kIsTrueFV = IsInNDFV(vtx_x,vtx_y,vtx_z); //seems like it's always true for ND
    if(vtx_x<200)
      kIsOutOfTheDesert=true;
    else
      kIsOutOfTheDesert=false;

    if(nuPDG<0)
      kIsAntiNu=true;
    else
      kIsAntiNu=false;

    if(kIsNumuCC && !kIsAntiNu && kIsTrueFV && kIsOutOfTheDesert)
      SignalLikeOnlyNuMuCC = true;
    else
      SignalLikeOnlyNuMuCC = false;

    if(SignalLikeOnlyNuMuCC && reco_numu && reco_q == -1)
      OnlyNuMuRecoAndQ = true;
    else
      OnlyNuMuRecoAndQ = false;

    if(SignalLikeOnlyNuMuCC && reco_numu && (muon_contained || muon_tracker) && reco_q == -1 && Ehad_veto<30) //isonlynumucc has to be here because we care for after subtracting bkg
      CompleteRecoSelectionNu = true;
    else
      CompleteRecoSelectionNu = false;


    if(SignalLikeOnlyNuMuCC){
      EventsInVtxX->Fill(vtx_x);
      EventsInVtxXvsOffAxisslice->Fill(abspos_x, vtx_x);

    // EventsInVtxY->Fill(vtx_y);
    // EventsInVtxZ->Fill(vtx_z);

      if (Ehad_veto<30){
        EventsInVtxXWithVetoCut->Fill(vtx_x);
        EventsInVtxXvsOffAxissliceWithVetoCut->Fill(abspos_x, vtx_x);
        // EventsInVtxYWithVetoCut->Fill(vtx_y);
        // EventsInVtxZWithVetoCut->Fill(vtx_z);
      }
      int binIndexDetSlice = EventsInVtxXvsOffAxisslice->GetYaxis()->FindBin(vtx_x) - 1;
      HistTrueEnuDetPosOffAxis->Fill(Ev);
      HistTrueEnuDetSlice[binIndexDetSlice]->Fill(Ev);
    }


    //
    // EvDistribALLSimEvents->Fill(Ev);
    // EvVsEvRecoDistribALLSimEvents->Fill(Ev,Ev_reco);
    // nuPDGALLSimEvents->Fill(nuPDG);
    // EvRecoDistribALLSimEvents->Fill(Ev_reco);
    // TrueEvisALLSimEvents->Fill(VisEtrue);
    // RecoEvisALLSimEvents->Fill(VisERecoND);
    //
    //
    // if(SignalLikeOnlyNuMuCC){
    //   EvDistribMatch->Fill(Ev);
    //   EvRecoDistribMatch->Fill(Ev_reco);
    //   EvVsEvRecoDistribMatch->Fill(Ev, Ev_reco);
    //   nuPDGMatch->Fill(nuPDG);
    //   TrueEvisMatch->Fill(VisEtrue);
    //   RecoEvisMatch->Fill(VisERecoND);
    // }
    // if(!kIsAntiNu && kIsTrueFV && kIsOutOfTheDesert){
    //   EvDistribNoRecoCut->Fill(Ev);
    //   nuPDGNoRecoCut->Fill(nuPDG);
    //   if(!isCC)
    //     EvDistribNoRecoCutNC->Fill(Ev);
    //   if(isCC && nuPDG!=14)
    //     EvDistribNoRecoCutWrongLepton->Fill(Ev);
    // }
    //
    // if( kIsTrueFV && kIsOutOfTheDesert){
    //   EvDistribNoRecoCutNoKantiNuCut->Fill(Ev);
    // }
    // if(CompleteRecoSelectionNu){
    //   EvDistribALLSelectionCuts->Fill(Ev);
    //   EvRecoDistribALLSelectionCuts->Fill(Ev_reco);
    //   EvVsEvRecoDistribALLSelectionCuts->Fill(Ev, Ev_reco);
    //   TrueEvisALLSelectionCuts->Fill(VisEtrue);
    //   RecoEvisALLSelectionCuts->Fill(VisERecoND);
    // }
    //
    // if(OnlyNuMuRecoAndQ){
    //   EvDistribOnlyNumRecoAndQ->Fill(Ev);
    //   EvRecoDistribOnlyNumRecoAndQ->Fill(Ev_reco);
    //   EvVsEvRecoDistribOnlyNumRecoAndQ->Fill(Ev, Ev_reco);
    //   TrueEvisOnlyNumRecoAndQ->Fill(VisEtrue);
    //   RecoEvisOnlyNumRecoAndQ->Fill(VisERecoND);
    // }






    // if(abs(nuPDGunosc)==14)
    //   cout<<"nuPDGunosc = "<<nuPDGunosc<<" nupdg "<<nuPDG<<endl;
    // if(!kIsOutOfTheDesert)
  //  cout<<" ev nr: "<<n<<" is NuMuCC? "<<kIsNumuCC<<" kIsTrueFV" <<kIsTrueFV<<" kisantinu: "<<kIsAntiNu<<endl;

  treeVars->Fill();
  }


  TCanvas* CanvasVtxX = new TCanvas("CanvasVtxX", "CanvasVtxX");
  CanvasVtxX->cd();
  EventsInVtxX->SetLineWidth(2);
  EventsInVtxX->SetLineColor(1);
  EventsInVtxX->Draw("hist");
  EventsInVtxXWithVetoCut->SetLineWidth(2);
  EventsInVtxXWithVetoCut->SetLineColor(2);
  EventsInVtxXWithVetoCut->Draw("histsames");
  TLegend* legVtxX = new TLegend();
  legVtxX->AddEntry(EventsInVtxX, "All Events");
  legVtxX->AddEntry(EventsInVtxXWithVetoCut, "Events With Ehad < 30");
  legVtxX->Draw("same");

  TH1D* SelectedOverGeneratedInVtxX = (TH1D*) EventsInVtxXWithVetoCut->Clone();
  SelectedOverGeneratedInVtxX->Divide(EventsInVtxX);
  SelectedOverGeneratedInVtxX->SetLineWidth(2);
  TCanvas* CanvasSelOverGen = new TCanvas("CanvasSelOverGen", "CanvasSelOverGen");
  SelectedOverGeneratedInVtxX->Draw("hist");

  TCanvas* CanvasOffAxisSliceVsDetSliceAllEvents = new TCanvas("CanvasOffAxisSliceVsDetSliceAllEvents", "CanvasOffAxisSliceVsDetSliceAllEvents");
  CanvasOffAxisSliceVsDetSliceAllEvents->cd();
  // EventsInVtxXvsOffAxisslice->SetMinimum(508000);
  // EventsInVtxXvsOffAxisslice->SetMaximum(540000);
  EventsInVtxXvsOffAxisslice->Draw("COLZ");

  TCanvas* CanvasOffAxisSliceVsDetSliceHadCutEvents = new TCanvas("CanvasOffAxisSliceVsDetSliceHadCutEvents", "CanvasOffAxisSliceVsDetSliceHadCutEvents");
  CanvasOffAxisSliceVsDetSliceHadCutEvents->cd();
  // EventsInVtxXvsOffAxissliceWithVetoCut->SetMinimum(330000);
  // EventsInVtxXvsOffAxissliceWithVetoCut->SetMaximum(364000);
  EventsInVtxXvsOffAxissliceWithVetoCut->Draw("COLZ");

  EventsInVtxXvsOffAxissliceWithVetoCut->Write();
  EventsInVtxXvsOffAxisslice->Write();
  SelectedOverGeneratedInVtxX->Write();
  EventsInVtxXWithVetoCut->Write();
  EventsInVtxX->Write();
  HistTrueEnuDetPosOffAxis->Write();
  for(int i = 0; i< nDetSlices; i++ ){
    HistTrueEnuDetSlice[i]->Write();
    delete HistTrueEnuDetSlice[i];
  }
  treeVars->Write();
  FileWithHistos->Write();
  FileWithHistos->Close();

  //
  // TCanvas* CanvasVtxY = new TCanvas("CanvasVtxY", "CanvasVtxY");
  // CanvasVtxY->cd();
  // EventsInVtxY->SetLineWidth(2);
  // EventsInVtxY->SetLineColor(1);
  // EventsInVtxY->Draw("hist");
  // EventsInVtxYWithVetoCut->SetLineWidth(2);
  // EventsInVtxYWithVetoCut->SetLineColor(2);
  // EventsInVtxYWithVetoCut->Draw("histsames");
  // TLegend* legVtxY = new TLegend();
  // legVtxY->AddEntry(EventsInVtxY, "All Events");
  // legVtxY->AddEntry(EventsInVtxYWithVetoCut, "Events With Ehad < 30");
  // legVtxY->Draw("same");
  //
  // TCanvas* CanvasVtxZ = new TCanvas("CanvasVtxZ", "CanvasVtxZ");
  // CanvasVtxZ->cd();
  // EventsInVtxZ->SetLineWidth(2);
  // EventsInVtxZ->SetLineColor(1);
  // EventsInVtxZ->Draw("hist");
  // EventsInVtxZWithVetoCut->SetLineWidth(2);
  // EventsInVtxZWithVetoCut->SetLineColor(2);
  // EventsInVtxZWithVetoCut->Draw("histsames");
  // TLegend* legVtxZ = new TLegend();
  // legVtxZ->AddEntry(EventsInVtxZ, "All Events");
  // legVtxZ->AddEntry(EventsInVtxZWithVetoCut, "Events With Ehad < 30");
  // legVtxZ->Draw("same");
  //
  // kPRISMNDSignal_True_numu =
  //   kIsNumuCC && (!kIsAntiNu) && ana::kIsTrueFV && kIsOutOfTheDesert;


}
