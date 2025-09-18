#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <vector>
#include <cmath>

void AnaKFPtree(
    // const char* filename="/mnt/e/sphenix/HF-QQ/output/K0s_reco_4M_0918.root",
    const char* filename="/mnt/e/sphenix/HF-QQ/output/K0s_bkgreco_4M_0918.root",
    const char* treename="DecayTree",
    // const char* outputname="/mnt/e/sphenix/HF-QQ/output/K0s_reco_4M_0918_ana.root"
    const char* outputname="/mnt/e/sphenix/HF-QQ/output/K0s_bkgreco_4M_0918_ana.root"
) 
{
    TFile outfile(outputname, "RECREATE");
    TH1D *h1_mass = new TH1D("h1_mass", "h1_mass; mass [GeV]; count", 100, 0, 1);
    
    TFile *file_in = TFile::Open(filename);
    if (!file_in || file_in->IsZombie()) {
        std::cerr<<"Error: cannot open "<<filename<<std::endl;
        return;
    }

    TTree *tree_data = (TTree*)file_in->Get(treename);
    if (!tree_data) {
        std::cerr<<"Error: cannot find tree "<<treename<<std::endl;
        return;
    }
               
    // branches read variables
    // std::vector <int> *clus_system = nullptr;
    // std::vector <double> *clus_X = nullptr;
    float K_S0_mass;

    tree_data->SetBranchAddress("K_S0_mass", & K_S0_mass);

    Long64_t nentries = tree_data->GetEntries();
    for (Long64_t i=0; i<nentries; ++i) 
    {
        tree_data->GetEntry(i);

        h1_mass->Fill(K_S0_mass);
    }

    outfile.cd();
    h1_mass->Write();
}


