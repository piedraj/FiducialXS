#define FiducialXS_cxx
#include "FiducialXS.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>


enum {NONE=-1, EE=0, MM=1, EM=2};

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

  Long64_t nentries = 20000;//fChain->GetEntries();

  printf("\n Will run on %lld events\n\n", nentries);

  for (Long64_t jentry=0; jentry<nentries; jentry++) {

    h_events->Fill(1);

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    if (jentry%10000 == 0) std::cout << "." << std::flush;

    h_ttbar->Fill(1);


    // Count emu GEN events
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


    // Build the GEN TLorentzVectors
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


    // Get the good reconstructed leptons
    //--------------------------------------------------------------------------
    std::vector<TLorentzVector> vElec;
    std::vector<TLorentzVector> vMuon;

    std::vector<int> cElec;
    std::vector<int> cMuon;

    int n_loose_lepton = 0;

    for (UInt_t i=0; i<T_Elec_Pt->size(); i++) {

      if (isTightElec(i))
	{
	  TLorentzVector gElec(T_Elec_Px->at(i), T_Elec_Py->at(i),
			       T_Elec_Pz->at(i), T_Elec_Energy->at(i));

	  vElec.push_back(gElec);
	  cElec.push_back(T_Elec_Charge->at(i));
	}
      else if (isLooseElec(i)) n_loose_lepton++;
    }

    for (UInt_t i=0; i<T_Muon_Pt->size(); i++) {

      if (isTightMuon(i))
	{
	  TLorentzVector gMuon(T_Muon_Px->at(i), T_Muon_Py->at(i),
			       T_Muon_Pz->at(i), T_Muon_Energy->at(i));

	  vMuon.push_back(gMuon);
	  cMuon.push_back(T_Muon_Charge->at(i));
	}
      else if (isLooseMuon(i)) n_loose_lepton++;
    }


    // Select the emu channel
    //--------------------------------------------------------------------------
    int n_elec = vElec.size();
    int n_muon = vMuon.size();
    
    bool  opposite_sign = false;
    float mll           = -999;
    int   njet          = 0;
    int   nbjet         = 0;

    if (n_elec == 1 && n_muon == 1)
      {
	opposite_sign = (cElec[0] * cMuon[0] < 0);

	mll = (vElec[0] + vMuon[0]).M();


	// Jet selection
	//----------------------------------------------------------------------
	for (UInt_t i=0 ; i<T_JetAKCHS_Et->size(); i++) {

	  TLorentzVector jet(T_JetAKCHS_Px->at(i), T_JetAKCHS_Py->at(i),
                             T_JetAKCHS_Pz->at(i), T_JetAKCHS_Energy->at(i));
	  
	  if (T_JetAKCHS_Et->at(i) > 30. &&
	      fabs(T_JetAKCHS_Eta->at(i)) < 2.4 &&
	      (jet.DeltaR(vElec[0]) >= 0.4) &&
	      (jet.DeltaR(vMuon[0]) >= 0.4) &&
	      passJetID(i)) {

            njet++;

            bool is_bjet = T_JetAKCHS_Tag_CombInclusiveSVtxV2->at(i) > 0.423;
            if (is_bjet) nbjet++;
	  }
	}
      }

    
    // Get the selected yields
    //--------------------------------------------------------------------------
    if (n_elec == 1 &&
	n_muon == 1 &&
	n_loose_lepton == 0 &&
	opposite_sign &&
	mll > 20 &&
	njet > 1 &&
	nbjet > 0)
      {
	h_ttbar_selected->Fill(1);

	if (is_ttbar_emu_fiducial) h_ttbar_emu_fiducial_selected->Fill(1);
	else if (is_ttbar_emu) h_ttbar_emu_nonfiducial_selected->Fill(1);
      }
  }


  // Write the output
  //----------------------------------------------------------------------------
  float n_events                         = h_events                        ->Integral();
  float n_ttbar                          = h_ttbar                         ->Integral();
  float n_ttbar_selected                 = h_ttbar_selected                ->Integral();
  float n_ttbar_emu                      = h_ttbar_emu                     ->Integral();
  float n_ttbar_emu_fiducial             = h_ttbar_emu_fiducial            ->Integral();
  float n_ttbar_emu_fiducial_selected    = h_ttbar_emu_fiducial_selected   ->Integral();
  float n_ttbar_emu_nonfiducial_selected = h_ttbar_emu_nonfiducial_selected->Integral();

  float eff      = (n_ttbar > 0) ? n_ttbar_selected / n_ttbar : -999;
  float eff_fid  = (n_ttbar_emu_fiducial > 0) ? n_ttbar_emu_fiducial_selected / n_ttbar_emu_fiducial : -999;
  float f        = (n_ttbar_emu_fiducial_selected > 0) ? n_ttbar_emu_nonfiducial_selected / n_ttbar_emu_fiducial_selected : -999;
  float xs_ratio = eff / eff_fid / (1. + f);

  printf("\n\n");
  printf("--------------------------------------------------\n");
  printf(" N(events)                             = %.0f\n", n_events);
  printf("--------------------------------------------------\n");
  printf(" N(ttbar inclusive)                    = %.0f\n", n_ttbar);
  printf(" N(ttbar inclusive selected)           = %.0f\n", n_ttbar_selected);
  printf(" N(ttbar -> emu selected)              = %.0f (not used)\n", n_ttbar_emu_fiducial_selected + n_ttbar_emu_nonfiducial_selected);
  printf(" total efficiency eff                  = %5.2f%s\n", 1e2 * eff, "%");
  printf("--------------------------------------------------\n");
  printf(" N(ttbar -> emu)                       = %.0f\n", n_ttbar_emu);
  printf(" N(ttbar -> emu fiducial)              = %.0f\n", n_ttbar_emu_fiducial);
  printf(" N(ttbar -> emu fiducial selected)     = %.0f\n", n_ttbar_emu_fiducial_selected);
  printf(" N(ttbar -> emu non-fiducial selected) = %.0f\n", n_ttbar_emu_nonfiducial_selected);
  printf(" fiducial efficiency eff_fid           = %5.2f%s\n", 1e2 * eff_fid, "%");
  printf(" contamination fraction f              = %5.2f%s\n", 1e2 * f, "%");
  printf("--------------------------------------------------\n");
  printf(" xs_fid/xs = eff/eff_fid/(1+f)         = %.4f\n", xs_ratio);
  printf("--------------------------------------------------\n");
  printf("\n");
  
  output_file->Write();
  output_file->Close();
}


bool FiducialXS::isTightMuon( unsigned iMuon, float minPt) 
{
  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < minPt)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  
  // POG Tight Muons definition              
  if (!T_Muon_IsGlobalMuon->at(iMuon)                            ) return false;
  if (!T_Muon_IsPFMuon->at(iMuon)                                ) return false;
  if (T_Muon_NormChi2GTrk->at(iMuon)                       >= 10.) return false;
  if (T_Muon_NValidHitsGTrk->at(iMuon)                     <= 0  ) return false;
  if (T_Muon_NumOfMatchedStations->at(iMuon)               <= 1  ) return false;                     
  //if (TMath::Abs(T_Muon_IPwrtAveBSInTrack->at(iMuon))      >= 0.2) return false; 
  //if (TMath::Abs(T_Muon_vz->at(iMuon) - T_Vertex_z->at(0)) >= 0.5) return false;
  if (TMath::Abs(T_Muon_dxyInTrack->at(iMuon))             >= 0.2) return false; 
  if (TMath::Abs(T_Muon_dzInTrack ->at(iMuon))             >= 0.5) return false;
  if (T_Muon_NLayers->at(iMuon)                            <= 5  ) return false;
  if (T_Muon_NValidPixelHitsInTrk->at(iMuon)               <= 0  ) return false;

  float relIso = muonIsolation(iMuon);
  
  if (relIso > 0.12) return false;
  
  return true;
}

bool FiducialXS::isLooseMuon( unsigned iMuon, float minPt) 
{
  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < minPt)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  float relIso = muonIsolation(iMuon);
  if (relIso > 0.20)               return false; 
  if (T_Muon_IsPFMuon->at(iMuon) == 0) return false;
  if (T_Muon_IsGlobalMuon->at(iMuon) == 0 && T_Muon_IsTrackerMuonArbitrated->at(iMuon) == 0) return false; 
  return true;
}


float FiducialXS::muonIsolation(unsigned iMuon)
{
  float isolation = T_Muon_chargedHadronIsoR04->at(iMuon) +
    std::max(0.0,
	     T_Muon_neutralHadronIsoR04->at(iMuon) +
	     T_Muon_photonIsoR04->at(iMuon) -
	     0.5*T_Muon_sumPUPtR04->at(iMuon));

  float relative_isolation = isolation / T_Muon_Pt->at(iMuon);

  return relative_isolation;
}


bool FiducialXS::isTightElec(unsigned int iElec, float minPt){
  
  TLorentzVector lep( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                      T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  if (lep.Pt() < minPt)                return false;
  if (TMath::Abs(lep.Eta()) > 2.5)     return false;

  float relIso =  elecIsolation(iElec);

  // Medium ID requirements
  bool passMediumID = false;
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.010399 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.007641 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.032643 &&    
       T_Elec_HtoE->at(iElec)                               < 0.060662 && 
       relIso                                               < 0.097213 &&
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) - 
                  T_Elec_eSuperClusterOverP->at(iElec) / 
                  T_Elec_ecalEnergy->at(iElec))             < 0.153897 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.011811 && 
       //TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.070775 &&
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.070775 &&
       T_Elec_nLost->at(iElec)                              <= 1       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
      passMediumID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.029524 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.009285 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.042447 &&    
       T_Elec_HtoE->at(iElec)                               < 0.104263 && 
       relIso                                               < 0.116708 &&
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) - 
       T_Elec_eSuperClusterOverP->at(iElec) /
       T_Elec_ecalEnergy->at(iElec))                        < 0.137468 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.051682 && 
       //TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.180720 &&
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.180720 &&
       T_Elec_nLost->at(iElec)                              <= 1       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
       passMediumID = true;
  }
  
  if (!passMediumID) return false;
   
  return true;
}

bool FiducialXS::isLooseElec(unsigned int iElec, float minPt){
  
  TLorentzVector lep( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                      T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  if (lep.Pt() < minPt)                return false;
  if (TMath::Abs(lep.Eta()) > 2.5)     return false;

  float relIso =  elecIsolation(iElec);

  bool passVetoID = false;
  // veto ID requirements
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.011100 &&
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.016315 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.252044 &&       
       T_Elec_HtoE->at(iElec)                               < 0.345843 && 
       relIso                                               < 0.164369 &&
       ( 1.0/T_Elec_ecalEnergy->at(iElec) -
         T_Elec_eSuperClusterOverP->at(iElec) /
         T_Elec_ecalEnergy->at(iElec))                      < 0.248070 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.060279 && 
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.800538 &&
       T_Elec_nLost->at(iElec)                              <= 2       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
      passVetoID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.033987 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.010671 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.245263 &&    
       T_Elec_HtoE->at(iElec)                               < 0.134691 && 
       relIso                                               < 0.212604 &&
       ( 1.0/T_Elec_ecalEnergy->at(iElec) -
         T_Elec_eSuperClusterOverP->at(iElec) /
         T_Elec_ecalEnergy->at(iElec))                      < 0.157160 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.273097 && 
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.885860 &&
       T_Elec_nLost->at(iElec)                              <= 3       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
       passVetoID = true;
  }  
  
  if (!passVetoID) return false;
  
  return true;
}


float FiducialXS::elecIsolation(unsigned iElec)
{
  float isolation = T_Elec_sumChargedHadronPt->at(iElec) +
    std::max(0.0,
	     T_Elec_sumNeutralHadronEt->at(iElec) +
	     T_Elec_sumPhotonEt->at(iElec) -
	     0.5*T_Elec_sumPUPt->at(iElec));

  float relative_isolation = isolation / T_Elec_Pt->at(iElec);

  return relative_isolation;
}


bool FiducialXS::passJetID(unsigned iJet) {

  if ( !(T_JetAKCHS_nDaughters->at(iJet)        > 1   ) ) return false;
  if ( !(T_JetAKCHS_NeutHadEnergyFrac->at(iJet) < 0.99) ) return false;
  if ( !(T_JetAKCHS_NeutEmEnergyFrac ->at(iJet) < 0.99) ) return false;
  if (TMath::Abs(T_JetAKCHS_Eta->at(iJet)) < 2.5){
  if ( !(T_JetAKCHS_CharEmEnergyFrac->at(iJet)  < 0.99) ) return false;
  if ( !(T_JetAKCHS_CharHadEnergyFrac->at(iJet) > 0.00) ) return false;
  if ( !(T_JetAKCHS_ChargedMultiplicity->at(iJet) > 0 ) ) return false;
  return true;
  }
  return false;
}
