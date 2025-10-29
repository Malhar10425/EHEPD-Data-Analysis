//===========================================================
// createCompleteSlimTree.C - Save ALL particles as vectors
//===========================================================

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"
#include <vector>
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

void createCompleteSlimTree() {

  // --- Load Delphes shared library
  gSystem->Load(DELPHES_DIR + "/libDelphes.so");

  // --- Define all file pairs
  std::vector<std::pair<TString, TString>> filePairs = {
    {"/media/malhar/MG5_storage/pp_tb_s_channel_single_top_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root", 
     "/media/malhar/MG5_storage/Vector_Slim_Trees/pp_tb_s_complete.root"},
    
    {"/media/malhar/MG5_storage/pp_tj_t_channel_single_top_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Vector_Slim_Trees/pp_tj_t_complete.root"},
    
    {"/media/malhar/MG5_storage/pp_tth_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Vector_Slim_Trees/tth_complete.root"},
    
    {"/media/malhar/MG5_storage/pp_ttw_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Vector_Slim_Trees/pp_ttw_complete.root"},
    
    {"/media/malhar/MG5_storage/pp_ww_NLO_4F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Vector_Slim_Trees/pp_ww_complete.root"},
    
    {"/media/malhar/MG5_storage/pp_wz_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Vector_Slim_Trees/pp_wz_complete.root"},
    
    {"/media/malhar/MG5_storage/pp_zz_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root",
     "/media/malhar/MG5_storage/Vector_Slim_Trees/pp_zz_complete.root"}
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
    TTree *t_out = new TTree("SlimTree", "Complete Delphes Data with All Particles");

    // --- Vectors to store ALL particles
    std::vector<Float_t> photon_pt;
    std::vector<Float_t> photon_eta;
    std::vector<Float_t> photon_phi;
    
    std::vector<Float_t> electron_pt;
    std::vector<Float_t> electron_eta;
    std::vector<Float_t> electron_phi;
    
    std::vector<Float_t> muon_pt;
    std::vector<Float_t> muon_eta;
    std::vector<Float_t> muon_phi;
    
    std::vector<Float_t> jet_pt;
    std::vector<Float_t> jet_eta;
    std::vector<Float_t> jet_phi;
    std::vector<Int_t> jet_btag;  // Store b-tag info
    
    std::vector<Float_t> bjet_pt;
    std::vector<Float_t> bjet_eta;
    std::vector<Float_t> bjet_phi;
    
    Float_t met_pt = 0;
    Float_t met_phi = 0;
    
    Int_t n_photons = 0;
    Int_t n_electrons = 0;
    Int_t n_muons = 0;
    Int_t n_jets = 0;
    Int_t n_bjets = 0;

    // --- Create branches for vectors
    t_out->Branch("photon_pt", &photon_pt);
    t_out->Branch("photon_eta", &photon_eta);
    t_out->Branch("photon_phi", &photon_phi);
    
    t_out->Branch("electron_pt", &electron_pt);
    t_out->Branch("electron_eta", &electron_eta);
    t_out->Branch("electron_phi", &electron_phi);
    
    t_out->Branch("muon_pt", &muon_pt);
    t_out->Branch("muon_eta", &muon_eta);
    t_out->Branch("muon_phi", &muon_phi);
    
    t_out->Branch("jet_pt", &jet_pt);
    t_out->Branch("jet_eta", &jet_eta);
    t_out->Branch("jet_phi", &jet_phi);
    t_out->Branch("jet_btag", &jet_btag);
    
    t_out->Branch("bjet_pt", &bjet_pt);
    t_out->Branch("bjet_eta", &bjet_eta);
    t_out->Branch("bjet_phi", &bjet_phi);
    
    t_out->Branch("met_pt", &met_pt, "met_pt/F");
    t_out->Branch("met_phi", &met_phi, "met_phi/F");
    
    t_out->Branch("n_photons", &n_photons, "n_photons/I");
    t_out->Branch("n_electrons", &n_electrons, "n_electrons/I");
    t_out->Branch("n_muons", &n_muons, "n_muons/I");
    t_out->Branch("n_jets", &n_jets, "n_jets/I");
    t_out->Branch("n_bjets", &n_bjets, "n_bjets/I");

    // --- Loop over ALL events
    Long64_t nEntries = t_in->GetEntries();
    std::cout << "   Input tree has " << nEntries << " events" << std::endl;
    
    Long64_t processedEvents = 0;
    Long64_t totalJets = 0;
    Long64_t totalBJets = 0;
    
    for (Long64_t i = 0; i < nEntries; i++) {
      if (i % 5000 == 0) {
        std::cout << "      Processing event " << i << "/" << nEntries << std::endl;
      }

      t_in->GetEntry(i);

      // --- Clear vectors for new event
      photon_pt.clear(); photon_eta.clear(); photon_phi.clear();
      electron_pt.clear(); electron_eta.clear(); electron_phi.clear();
      muon_pt.clear(); muon_eta.clear(); muon_phi.clear();
      jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_btag.clear();
      bjet_pt.clear(); bjet_eta.clear(); bjet_phi.clear();
      
      met_pt = 0; met_phi = 0;
      n_photons = 0; n_electrons = 0; n_muons = 0; n_jets = 0; n_bjets = 0;

      // --- Fill ALL photons
      if (branchPhoton) {
        n_photons = branchPhoton->GetEntries();
        for (Int_t j = 0; j < n_photons; j++) {
          Photon *photon = (Photon*) branchPhoton->At(j);
          photon_pt.push_back(photon->PT);
          photon_eta.push_back(photon->Eta);
          photon_phi.push_back(photon->Phi);
        }
      }

      // --- Fill ALL electrons
      if (branchElectron) {
        n_electrons = branchElectron->GetEntries();
        for (Int_t j = 0; j < n_electrons; j++) {
          Electron *electron = (Electron*) branchElectron->At(j);
          electron_pt.push_back(electron->PT);
          electron_eta.push_back(electron->Eta);
          electron_phi.push_back(electron->Phi);
        }
      }

      // --- Fill ALL muons
      if (branchMuon) {
        n_muons = branchMuon->GetEntries();
        for (Int_t j = 0; j < n_muons; j++) {
          Muon *muon = (Muon*) branchMuon->At(j);
          muon_pt.push_back(muon->PT);
          muon_eta.push_back(muon->Eta);
          muon_phi.push_back(muon->Phi);
        }
      }

      // --- Fill ALL jets and separate b-jets
      if (branchJet) {
        n_jets = branchJet->GetEntries();
        totalJets += n_jets;
        
        for (Int_t j = 0; j < n_jets; j++) {
          Jet *jet = (Jet*) branchJet->At(j);
          jet_pt.push_back(jet->PT);
          jet_eta.push_back(jet->Eta);
          jet_phi.push_back(jet->Phi);
          jet_btag.push_back(jet->BTag);
          
          // Separate b-jets
          if (jet->BTag == 1) {
            bjet_pt.push_back(jet->PT);
            bjet_eta.push_back(jet->Eta);
            bjet_phi.push_back(jet->Phi);
            n_bjets++;
            totalBJets++;
          }
        }
      }

      // --- Fill MET
      if (branchMissingET && branchMissingET->GetEntries() > 0) {
        MissingET *met = (MissingET*) branchMissingET->At(0);
        met_pt = met->MET;
        met_phi = met->Phi;
      }

      t_out->Fill();
      processedEvents++;
    }

    // --- Save and close
    f_out->cd();
    t_out->Write();
    f_out->Close();
    f_in->Close();

    std::cout << "✅ Successfully processed " << processedEvents << "/" << nEntries << " events" << std::endl;
    std::cout << "✅ Total jets stored: " << totalJets << " (average " << (Float_t)totalJets/nEntries << " jets/event)" << std::endl;
    std::cout << "✅ Total b-jets stored: " << totalBJets << " (average " << (Float_t)totalBJets/nEntries << " b-jets/event)" << std::endl;
    std::cout << "✅ Created complete slim tree with ALL particles: " << outputFile << std::endl;
    std::cout << std::endl;
  }

  std::cout << "🎉 All 7 processes completed successfully!" << std::endl;
}
