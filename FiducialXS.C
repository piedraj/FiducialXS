#define FiducialXS_cxx
#include "FiducialXS.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
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

  output_file = new TFile("rootfiles/fiducial" + suffix + ".root", "recreate");

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

  Long64_t nentries = 10000;//fChain->GetEntries();

  for (Long64_t jentry=0; jentry<nentries; jentry++) {

    h_events->Fill(1);

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    h_ttbar->Fill(1);


    // Count emu events
    //--------------------------------------------------------------------------
    int n_gen_electron = T_Gen_PromptElec_pdgId->size();
    int n_gen_muon     = T_Gen_PromptMuon_pdgId->size();

    int n_gen_tau_electron = 0;
    int n_gen_tau_muon     = 0;

    int gen_tau_electron_index = -1;
    int gen_tau_muon_index     = -1;

    for (UInt_t i=0; i<T_Gen_PromptTau_pdgId->size(); i++) {
      if (T_Gen_PromptTau_IsLepDec->at(i)) {
	if (abs(T_Gen_PromptTau_LepDec_pdgId->at(i)) == 11)
	  {
	    n_gen_tau_electron++;
	    gen_tau_electron_index = i;
	  }
	if (abs(T_Gen_PromptTau_LepDec_pdgId->at(i)) == 13)
	  {
	    n_gen_tau_muon++;
	    gen_tau_muon_index = i;
	  }
      }
    }


    // Build the gen TLorentzVectors
    //--------------------------------------------------------------------------
    TLorentzVector gen_electron;
    TLorentzVector gen_muon;
	
    int gen_electron_pdgId = 0;
    int gen_muon_pdgId     = 0;

    if (n_gen_electron + n_gen_tau_electron == 1)
      {
	if (n_gen_electron == 1)
	  {
	    gen_electron.SetPxPyPzE(T_Gen_PromptElec_Px->at(0),
				    T_Gen_PromptElec_Py->at(0),
				    T_Gen_PromptElec_Pz->at(0),
				    T_Gen_PromptElec_Energy->at(0));

	    gen_electron_pdgId = T_Gen_PromptElec_pdgId->at(0);
	  }
	else if (n_gen_tau_electron == 1)
	  {
	    gen_electron.SetPxPyPzE(T_Gen_PromptTau_LepDec_Px->at(gen_tau_electron_index),
				    T_Gen_PromptTau_LepDec_Py->at(gen_tau_electron_index),
				    T_Gen_PromptTau_LepDec_Pz->at(gen_tau_electron_index),
				    T_Gen_PromptTau_LepDec_Energy->at(gen_tau_electron_index));

	    gen_electron_pdgId = T_Gen_PromptTau_LepDec_pdgId->at(gen_tau_electron_index);
	  }
      }

    if (n_gen_muon + n_gen_tau_muon == 1)
      {
	if (n_gen_muon == 1)
	  {
	    gen_muon.SetPxPyPzE(T_Gen_PromptMuon_Px->at(0),
				T_Gen_PromptMuon_Py->at(0),
				T_Gen_PromptMuon_Pz->at(0),
				T_Gen_PromptMuon_Energy->at(0));
	    
	    gen_muon_pdgId = T_Gen_PromptMuon_pdgId->at(0);
	  }
	else if (n_gen_tau_muon == 1)
	  {
	    gen_muon.SetPxPyPzE(T_Gen_PromptTau_LepDec_Px->at(gen_tau_muon_index),
				T_Gen_PromptTau_LepDec_Py->at(gen_tau_muon_index),
				T_Gen_PromptTau_LepDec_Pz->at(gen_tau_muon_index),
				T_Gen_PromptTau_LepDec_Energy->at(gen_tau_muon_index));
	    
	    gen_muon_pdgId = T_Gen_PromptTau_LepDec_pdgId->at(gen_tau_muon_index);
	  }
      }


    // Decide if the event is emu, fiducial
    //--------------------------------------------------------------------------
    bool is_ttbar_emu = (gen_electron_pdgId * gen_muon_pdgId < 0);

    bool is_ttbar_emu_fiducial = is_ttbar_emu;
    
    if (is_ttbar_emu)
      {
	is_ttbar_emu_fiducial &= (gen_electron.Pt() > 20);
	is_ttbar_emu_fiducial &= (fabs(gen_electron.Eta()) < 2.4);
	is_ttbar_emu_fiducial &= (gen_muon.Pt() > 20);
	is_ttbar_emu_fiducial &= (fabs(gen_muon.Eta()) < 2.4);
      }


    if (is_ttbar_emu)          h_ttbar_emu->Fill(1);
    if (is_ttbar_emu_fiducial) h_ttbar_emu_fiducial->Fill(1);
  }


  // Write the output
  //----------------------------------------------------------------------------
  printf("\n");
  printf(" events:             %.0f\n", h_events->Integral());
  printf(" ttbar:              %.0f\n", h_ttbar->Integral());
  printf(" ttbar_emu:          %.0f\n", h_ttbar_emu->Integral());
  printf(" ttbar_emu_fiducial: %.0f\n", h_ttbar_emu_fiducial->Integral());
  printf("\n");

  output_file->Write();
  output_file->Close();
}
