//===========================================================
// makeSlimTree.C  (Auto Delphes path version)
//===========================================================

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>

// --- Detect Delphes path automatically (change if you rename folder)
const TString DELPHES_DIR = "/home/malhar/MG5_aMC_v3_6_4/Delphes";

// --- Add include paths for ROOT
R__ADD_INCLUDE_PATH(/home/malhar/MG5_aMC_v3_6_4/Delphes)
R__ADD_INCLUDE_PATH(/home/malhar/MG5_aMC_v3_6_4/Delphes/external)

// --- Include Delphes headers
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

void makeSlimTree() {

  // --- Load Delphes shared library
  gSystem->Load(DELPHES_DIR + "/libDelphes.so");

  // --- Open input ROOT file
  TFile *f_in = TFile::Open("/media/malhar/MG5_storage/pp_zz_NLO_5F_TuneCP5_pythia8/Events/run_01_decayed_1/delphes_CP5.root");
  if (!f_in || f_in->IsZombie()) {
    std::cerr << "❌ Error: cannot open delphes.root!" << std::endl;
    return;
  }

  // --- Get main tree
  TTree *t_in = (TTree*) f_in->Get("Delphes");
  if (!t_in) {
    std::cerr << "❌ Error: tree 'Delphes' not found!" << std::endl;
    return;
  }

  // --- Set up branches
  TClonesArray *branchPhoton = nullptr;
  TClonesArray *branchElectron = nullptr;
  TClonesArray *branchJet = nullptr;

  t_in->SetBranchAddress("Photon", &branchPhoton);
  t_in->SetBranchAddress("Electron", &branchElectron);
  t_in->SetBranchAddress("Jet", &branchJet);

  // --- Output file and tree
  TFile *f_out = new TFile("/media/malhar/MG5_storage/Silimed_Trees_1/pp_zz.root", "RECREATE");
  TTree *t_out = new TTree("SlimTree", "Reduced Delphes Tree");

  Float_t photon_pt = 0;
  Float_t electron_pt = 0;
  Int_t n_jets = 0;

  t_out->Branch("photon_pt", &photon_pt, "photon_pt/F");
  t_out->Branch("electron_pt", &electron_pt, "electron_pt/F");
  t_out->Branch("n_jets", &n_jets, "n_jets/I");

  // --- Loop over events
  Long64_t nEntries = t_in->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    t_in->GetEntry(i);

    photon_pt = 0;
    if (branchPhoton->GetEntries() > 0) {
      Photon *photon = (Photon*) branchPhoton->At(0);
      photon_pt = photon->PT;
    }

    electron_pt = 0;
    if (branchElectron->GetEntries() > 0) {
      Electron *electron = (Electron*) branchElectron->At(0);
      electron_pt = electron->PT;
    }

    n_jets = branchJet->GetEntries();

    t_out->Fill();
  }

  // --- Save and close
  f_out->cd();
  t_out->Write();
  f_out->Close();
  f_in->Close();

  std::cout << "✅ Slimmed tree saved to slimmed.root" << std::endl;
}
