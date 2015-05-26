#define FiducialXS_cxx
#include "FiducialXS.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>


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

  Long64_t nentries = fChain->GetEntries();

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

    for (UInt_t i=0; i<T_Gen_PromptTau_IsLepDec->size(); i++) {
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

    // Electron loop
    for (UInt_t i=0; i<T_Elec_Pt->size(); i++) {

      if (!isFiducialElec(i)) continue;

      if (isMediumElec(i))
	{
	  TLorentzVector gElec(T_Elec_Px->at(i), T_Elec_Py->at(i),
			       T_Elec_Pz->at(i), T_Elec_Energy->at(i));
	  
	  vElec.push_back(gElec);
	  cElec.push_back(T_Elec_Charge->at(i));
	}
      else if (isVetoElec(i)) n_loose_lepton++;
    }

    // Muon loop
    for (UInt_t i=0; i<T_Muon_Pt->size(); i++) {

      if (!isFiducialMuon(i)) continue;

      if (isTightMuon(i) && muonIsolation(i) < 0.12)
	{
	  TLorentzVector gMuon(T_Muon_Px->at(i), T_Muon_Py->at(i),
			       T_Muon_Pz->at(i), T_Muon_Energy->at(i));

	  vMuon.push_back(gMuon);
	  cMuon.push_back(T_Muon_Charge->at(i));
	}
      else if (isLooseMuon(i) && muonIsolation(i) < 0.20) n_loose_lepton++;
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
	      (jet.DeltaR(vElec[0]) > 0.4) &&
	      (jet.DeltaR(vMuon[0]) > 0.4) &&
	      passJetID(i)) {

            njet++;

	    //            bool is_bjet = T_JetAKCHS_Tag_CombInclusiveSVtxV2->at(i) > 0.423;
            bool is_bjet = T_JetAKCHS_Tag_CombSVtx->at(i) > 0.423;
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


  printf("\n");


  // Write the output
  //----------------------------------------------------------------------------
  output_file->Write();
  output_file->Close();
}


//------------------------------------------------------------------------------
// isFiducialMuon
//------------------------------------------------------------------------------
bool FiducialXS::isFiducialMuon(unsigned int iMuon) 
{
  bool isFiducial = true;

  isFiducial &= (T_Muon_Pt->at(iMuon)        > 20.);
  isFiducial &= (fabs(T_Muon_Eta->at(iMuon)) < 2.4);

  return isFiducial;
}


//------------------------------------------------------------------------------
// isTightMuon
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015#Tight_Muon
//------------------------------------------------------------------------------
bool FiducialXS::isTightMuon(unsigned int iMuon) 
{
  bool isTight = true;

  isTight &= (T_Muon_IsGlobalMuon->at(iMuon)              );
  isTight &= (T_Muon_IsPFMuon->at(iMuon)                  );
  isTight &= (T_Muon_NormChi2GTrk->at(iMuon)         < 10.);
  isTight &= (T_Muon_NValidHitsGTrk->at(iMuon)       > 0  );
  isTight &= (T_Muon_NumOfMatchedStations->at(iMuon) > 1  );                     
  isTight &= (fabs(T_Muon_dxyInTrack->at(iMuon))     < 0.2); 
  isTight &= (fabs(T_Muon_dzInTrack ->at(iMuon))     < 0.5);
  isTight &= (T_Muon_NValidPixelHitsInTrk->at(iMuon) > 0  );
  isTight &= (T_Muon_NLayers->at(iMuon)              > 5  );

  return isTight;
}


//------------------------------------------------------------------------------
// isLooseMuon
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015#Loose_Muon
//------------------------------------------------------------------------------
bool FiducialXS::isLooseMuon(unsigned int iMuon) 
{
  bool isLoose = true;

  isLoose &= (T_Muon_IsPFMuon->at(iMuon));
  isLoose &= (T_Muon_IsGlobalMuon->at(iMuon) || T_Muon_IsTrackerMuonArbitrated->at(iMuon));

  return isLoose;
}


//------------------------------------------------------------------------------
// muonIsolation
//------------------------------------------------------------------------------
float FiducialXS::muonIsolation(unsigned int iMuon)
{
  float isolation =
    T_Muon_chargedHadronIsoR04->at(iMuon) +
    std::max(0.0,
	     T_Muon_neutralHadronIsoR04->at(iMuon) +
	     T_Muon_photonIsoR04->at(iMuon) -
	     0.5*T_Muon_sumPUPtR04->at(iMuon));

  float relative_isolation = isolation / T_Muon_Pt->at(iMuon);

  return relative_isolation;
}


//------------------------------------------------------------------------------
// isFiducialElec
//------------------------------------------------------------------------------
bool FiducialXS::isFiducialElec(unsigned int iElec)
{
  bool isFiducial = true;

  isFiducial &= (T_Elec_Pt->at(iElec)        > 20.);
  isFiducial &= (fabs(T_Elec_Eta->at(iElec)) < 2.5);

  return isFiducial;
}


//------------------------------------------------------------------------------
// isMediumElec
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=14#PHYS14_selection_all_conditions
//------------------------------------------------------------------------------
bool FiducialXS::isMediumElec(unsigned int iElec)
{
  float sceta = fabs(T_Elec_SC_Eta->at(iElec));

  if (sceta > 1.4442 && sceta < 1.566) return false;

  float relIso =  elecIsolation(iElec);

  float ooEmooP = (1. - T_Elec_eSuperClusterOverP->at(iElec)) / T_Elec_ecalEnergy->at(iElec);

  bool isMedium = false;
  
  if (sceta < 1.479) {
    if (T_Elec_sigmaIetaIetaFull5by5->at(iElec) < 0.010399 && 
	fabs(T_Elec_deltaEtaIn->at(iElec))      < 0.007641 && 
	fabs(T_Elec_deltaPhiIn->at(iElec))      < 0.032643 &&    
	T_Elec_HtoE->at(iElec)                  < 0.060662 && 
	relIso                                  < 0.097213 &&
	fabs(ooEmooP)                           < 0.153897 &&
	fabs(T_Elec_IPwrtPV->at(iElec))         < 0.011811 && 
	fabs(T_Elec_dzwrtPV->at(iElec))         < 0.070775 &&
	T_Elec_nLost->at(iElec)                 < 2        && 
	T_Elec_passConversionVeto->at(iElec)    > 0)
      isMedium = true;
  }
  else {
    if (T_Elec_sigmaIetaIetaFull5by5->at(iElec) < 0.029524 && 
	fabs(T_Elec_deltaEtaIn->at(iElec))      < 0.009285 && 
	fabs(T_Elec_deltaPhiIn->at(iElec))      < 0.042447 &&    
	T_Elec_HtoE->at(iElec)                  < 0.104263 && 
	relIso                                  < 0.116708 &&
	fabs(ooEmooP)                           < 0.137468 &&
	fabs(T_Elec_IPwrtPV->at(iElec))         < 0.051682 && 
	fabs(T_Elec_dzwrtPV->at(iElec))         < 0.180720 &&
	T_Elec_nLost->at(iElec)                 < 2        && 
	T_Elec_passConversionVeto->at(iElec)    > 0)
      isMedium = true;
  }
  
  return isMedium;
}


//------------------------------------------------------------------------------
// isVetoElec
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=14#PHYS14_selection_all_conditions
//------------------------------------------------------------------------------
bool FiducialXS::isVetoElec(unsigned int iElec)
{
  float sceta = fabs(T_Elec_SC_Eta->at(iElec));

  if (sceta > 1.4442 && sceta < 1.566) return false;

  float relIso = elecIsolation(iElec);

  float ooEmooP = (1. - T_Elec_eSuperClusterOverP->at(iElec)) / T_Elec_ecalEnergy->at(iElec);

  bool isVeto = false;

  if (sceta < 1.479) {
    if (T_Elec_sigmaIetaIetaFull5by5->at(iElec) < 0.011100 &&
	fabs(T_Elec_deltaEtaIn->at(iElec))      < 0.016315 && 
	fabs(T_Elec_deltaPhiIn->at(iElec))      < 0.252044 &&       
	T_Elec_HtoE->at(iElec)                  < 0.345843 && 
	relIso                                  < 0.164369 &&
	fabs(ooEmooP)                           < 0.248070 &&
	fabs(T_Elec_IPwrtPV->at(iElec))         < 0.060279 && 
	fabs(T_Elec_dzwrtPV->at(iElec))         < 0.800538 &&
	T_Elec_nLost->at(iElec)                 < 3        && 
	T_Elec_passConversionVeto->at(iElec)    > 0)
      isVeto = true;
  }
  else {
    if (T_Elec_sigmaIetaIetaFull5by5->at(iElec) < 0.033987 && 
	fabs(T_Elec_deltaEtaIn->at(iElec))      < 0.010671 && 
	fabs(T_Elec_deltaPhiIn->at(iElec))      < 0.245263 &&    
	T_Elec_HtoE->at(iElec)                  < 0.134691 && 
	relIso                                  < 0.212604 &&
	fabs(ooEmooP)                           < 0.157160 &&
	fabs(T_Elec_IPwrtPV->at(iElec))         < 0.273097 && 
	fabs(T_Elec_dzwrtPV->at(iElec))         < 0.885860 &&
	T_Elec_nLost->at(iElec)                 < 4        && 
	T_Elec_passConversionVeto->at(iElec)    > 0)
      isVeto = true;
  }
  
  return isVeto;
}


//------------------------------------------------------------------------------
// elecIsolation
//------------------------------------------------------------------------------
float FiducialXS::elecIsolation(unsigned int iElec)
{
  float isolation =
    T_Elec_sumChargedHadronPt->at(iElec) +
    std::max(0.0,
	     T_Elec_sumNeutralHadronEt->at(iElec) +
	     T_Elec_sumPhotonEt->at(iElec) -
	     0.5*T_Elec_sumPUPt->at(iElec));

  float relative_isolation = isolation / T_Elec_Pt->at(iElec);

  return relative_isolation;
}


//------------------------------------------------------------------------------
// passJetID
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
//------------------------------------------------------------------------------
bool FiducialXS::passJetID(unsigned int iJet)
{
  Bool_t jetid = true;

  jetid &= (T_JetAKCHS_nDaughters->at(iJet)        > 1   );
  jetid &= (T_JetAKCHS_NeutHadEnergyFrac->at(iJet) < 0.99);
  jetid &= (T_JetAKCHS_NeutEmEnergyFrac ->at(iJet) < 0.99);

  if (fabs(T_JetAKCHS_Eta->at(iJet)) < 2.5)
    {
      jetid &= (T_JetAKCHS_CharEmEnergyFrac->at(iJet)    < 0.99);
      jetid &= (T_JetAKCHS_CharHadEnergyFrac->at(iJet)   > 0.  );
      jetid &= (T_JetAKCHS_ChargedMultiplicity->at(iJet) > 0   );
    }

  return jetid;
}
