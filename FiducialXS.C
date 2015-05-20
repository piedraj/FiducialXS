#define FiducialXS_cxx
#include "FiducialXS.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>


TFile* output_file;

TH1F* h_events;
TH1F* h_ttbar;
TH1F* h_ttbar_selected;
TH1F* h_ttbar_emu;
TH1F* h_ttbar_emu_fiducial;
TH1F* h_ttbar_emu_fiducial_selected;
TH1F* h_ttbar_emu_nonfiducial_selected;


void FiducialXS::Loop(Int_t index)
{
  gSystem->mkdir("rootfiles", kTRUE);

  TString suffix = (index < 0) ? "" : Form("_%d", index);

  output_file = new TFile("rootfiles/test" + suffix + ".root", "recreate");

  h_events                         = new TH1F("h_events",                         "", 3, 0, 3);
  h_ttbar                          = new TH1F("h_ttbar",                          "", 3, 0, 3);
  h_ttbar_selected                 = new TH1F("h_ttbar_selected",                 "", 3, 0, 3);
  h_ttbar_emu                      = new TH1F("h_ttbar_emu",                      "", 3, 0, 3);
  h_ttbar_emu_fiducial             = new TH1F("h_ttbar_emu_fiducial",             "", 3, 0, 3);
  h_ttbar_emu_fiducial_selected    = new TH1F("h_ttbar_emu_fiducial_selected",    "", 3, 0, 3);
  h_ttbar_emu_nonfiducial_selected = new TH1F("h_ttbar_emu_nonfiducial_selected", "", 3, 0, 3);


  // Loop
  //----------------------------------------------------------------------------
  if (fChain == 0) return;

  Long64_t nentries = 100;//fChain->GetEntriesFast();

  for (Long64_t jentry=0; jentry<nentries; jentry++) {

    h_events->Fill(1);

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;

    h_ttbar->Fill(1);
  }


  // Write the output
  //----------------------------------------------------------------------------
  output_file->Write();
  output_file->Close();
}
