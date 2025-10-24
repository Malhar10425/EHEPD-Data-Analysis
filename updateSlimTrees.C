//===========================================================
// updateSlimTrees.C - Update existing slim trees with new branches
//===========================================================

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include <vector>
#include <utility>

// --- Detect Delphes path automatically
const TString DELPHES_DIR = "/home/malhar/MG5_aMC_v3_6_4/Delphes";

// --- Add include paths for ROOT
R__ADD_INCLUDE_PATH(/home/malhar/MG5_aMC_v3_6_4/Delphes)
R__ADD_INCLUDE_PATH(/home/malhar/MG5_aMC_v3_6_4/Delphes/external)

// --- Include Delphes headers
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

void updateSlimTrees() {

  // --- Load Delphes shared library
  gSystem->Load(DELPHES_DIR + "/libDelphes.so");

  // --- Define all file pairs (input Delphes -> output SlimTree)
  std::vector<std::pair<TString, TString>> filePairs = {
    {"/media/malhar/MG5_storage/pp_tb_s_channel_single_top_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root", 
     "/media/malhar/MG5_storage/Silimed_Trees_1/pp_tb_s.root"},
    
    {"/media/malhar/MG5_storage/pp_tj_t_channel_single_top_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Silimed_Trees_1/pp_tj_t.root"},
    
    {"/media/malhar/MG5_storage/pp_tth_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Silimed_Trees_1/tth.root"},
    
    {"/media/malhar/MG5_storage/pp_ttw_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Silimed_Trees_1/pp_ttw.root"},
    
    {"/media/malhar/MG5_storage/pp_ww_NLO_4F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Silimed_Trees_1/pp_ww.root"},
    
    {"/media/malhar/MG5_storage/pp_wz_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Silimed_Trees_1/pp_wz.root"},
    
    {"/media/malhar/MG5_storage/pp_zz_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Silimed_Trees_1/pp_zz.root"}
  };

  // --- Process each file pair
  for (const auto& filePair : filePairs) {
    TString inputFile = filePair.first;
    TString outputFile = filePair.second;
    
    std::cout << "🔁 Processing: " << inputFile << std::endl;
    std::cout << "   Output: " << outputFile << std::endl;

    // --- Open input ROOT file
    TFile *f_in = TFile::Open(inputFile);
    if (!f_in || f_in->IsZombie()) {
      std::cerr << "❌ Error: cannot open " << inputFile << "!" << std::endl;
      continue;
    }

    // --- Get main tree
    TTree *t_in = (TTree*) f_in->Get("Delphes");
    if (!t_in) {
      std::cerr << "❌ Error: tree 'Delphes' not found in " << inputFile << "!" << std::endl;
      f_in->Close();
      continue;
    }

    // --- Set up input branches
    TClonesArray *branchPhoton = nullptr;
    TClonesArray *branchElectron = nullptr;
    TClonesArray *branchMuon = nullptr;
    TClonesArray *branchJet = nullptr;
    TClonesArray *branchMissingET = nullptr;

    t_in->SetBranchAddress("Photon", &branchPhoton);
    t_in->SetBranchAddress("Electron", &branchElectron);
    t_in->SetBranchAddress("Muon", &branchMuon);
    t_in->SetBranchAddress("Jet", &branchJet);
    t_in->SetBranchAddress("MissingET", &branchMissingET);

    // --- Output file and tree
    TFile *f_out = new TFile(outputFile, "RECREATE");
    TTree *t_out = new TTree("SlimTree", "Updated Reduced Delphes Tree");

    // --- EXISTING branches (from your original macro)
    Float_t photon_pt = 0;
    Float_t electron_pt = 0;
    Int_t n_jets = 0;

    // --- NEW branches to add
    Float_t electron_eta = 0;
    Float_t electron_phi = 0;
    Float_t photon_eta = 0;
    Float_t photon_phi = 0;
    Float_t muon_pt = 0;
    Float_t muon_eta = 0;
    Float_t muon_phi = 0;
    Float_t jet_pt = 0;
    Float_t jet_eta = 0;
    Float_t jet_phi = 0;
    Float_t bjet_pt = 0;
    Float_t bjet_eta = 0;
    Float_t bjet_phi = 0;
    Float_t met_pt = 0;
    Float_t met_phi = 0;

    // --- Create branches in output tree
    // Existing branches
    t_out->Branch("photon_pt", &photon_pt, "photon_pt/F");
    t_out->Branch("electron_pt", &electron_pt, "electron_pt/F");
    t_out->Branch("n_jets", &n_jets, "n_jets/I");
    
    // New branches
    t_out->Branch("electron_eta", &electron_eta, "electron_eta/F");
    t_out->Branch("electron_phi", &electron_phi, "electron_phi/F");
    t_out->Branch("photon_eta", &photon_eta, "photon_eta/F");
    t_out->Branch("photon_phi", &photon_phi, "photon_phi/F");
    t_out->Branch("muon_pt", &muon_pt, "muon_pt/F");
    t_out->Branch("muon_eta", &muon_eta, "muon_eta/F");
    t_out->Branch("muon_phi", &muon_phi, "muon_phi/F");
    t_out->Branch("jet_pt", &jet_pt, "jet_pt/F");
    t_out->Branch("jet_eta", &jet_eta, "jet_eta/F");
    t_out->Branch("jet_phi", &jet_phi, "jet_phi/F");
    t_out->Branch("bjet_pt", &bjet_pt, "bjet_pt/F");
    t_out->Branch("bjet_eta", &bjet_eta, "bjet_eta/F");
    t_out->Branch("bjet_phi", &bjet_phi, "bjet_phi/F");
    t_out->Branch("met_pt", &met_pt, "met_pt/F");
    t_out->Branch("met_phi", &met_phi, "met_phi/F");

    // --- Loop over events
    Long64_t nEntries = t_in->GetEntries();
    std::cout << "   Processing " << nEntries << " events..." << std::endl;
    
    for (Long64_t i = 0; i < nEntries; i++) {
      t_in->GetEntry(i);

      // --- Reset all variables for this event
      photon_pt = 0; photon_eta = 0; photon_phi = 0;
      electron_pt = 0; electron_eta = 0; electron_phi = 0;
      muon_pt = 0; muon_eta = 0; muon_phi = 0;
      jet_pt = 0; jet_eta = 0; jet_phi = 0;
      bjet_pt = 0; bjet_eta = 0; bjet_phi = 0;
      met_pt = 0; met_phi = 0;
      n_jets = branchJet->GetEntries();

      // --- Photon (leading)
      if (branchPhoton->GetEntries() > 0) {
        Photon *photon = (Photon*) branchPhoton->At(0);
        photon_pt = photon->PT;
        photon_eta = photon->Eta;
        photon_phi = photon->Phi;
      }

      // --- Electron (leading)
      if (branchElectron->GetEntries() > 0) {
        Electron *electron = (Electron*) branchElectron->At(0);
        electron_pt = electron->PT;
        electron_eta = electron->Eta;
        electron_phi = electron->Phi;
      }

      // --- Muon (leading)
      if (branchMuon->GetEntries() > 0) {
        Muon *muon = (Muon*) branchMuon->At(0);
        muon_pt = muon->PT;
        muon_eta = muon->Eta;
        muon_phi = muon->Phi;
      }

      // --- Jet (leading)
      if (branchJet->GetEntries() > 0) {
        Jet *jet = (Jet*) branchJet->At(0);
        jet_pt = jet->PT;
        jet_eta = jet->Eta;
        jet_phi = jet->Phi;
      }

      // --- B-jet (leading b-tagged jet)
      for (Int_t j = 0; j < branchJet->GetEntries(); j++) {
        Jet *jet = (Jet*) branchJet->At(j);
        if (jet->BTag == 1) { // Assuming BTag == 1 means b-jet
          bjet_pt = jet->PT;
          bjet_eta = jet->Eta;
          bjet_phi = jet->Phi;
          break; // Take first b-jet
        }
      }

      // --- Missing ET
      if (branchMissingET->GetEntries() > 0) {
        MissingET *met = (MissingET*) branchMissingET->At(0);
        met_pt = met->MET;
        met_phi = met->Phi;
      }

      t_out->Fill();
    }

    // --- Save and close
    f_out->cd();
    t_out->Write();
    f_out->Close();
    f_in->Close();

    std::cout << "✅ Updated slim tree saved to: " << outputFile << std::endl;
    std::cout << std::endl;
  }

  std::cout << "🎉 All 7 processes completed successfully!" << std::endl;
}
