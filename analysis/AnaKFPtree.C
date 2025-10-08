#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <vector>
#include <cmath>

void AnaKFPtree(
    const char* filename="/mnt/e/sphenix/HF-QQ/output/Conv_reco_4M.root",
    // const char* filename="/mnt/e/sphenix/HF-QQ/output/Conv_bkgreco_4M.root",
    // const char* filename="/sphenix/user/jzhang1/PhysicsAna/HF-QQ/output/PhotonConv/PhotonConv_reco_likesign/Reconstructed/53046/Conv_bkgreco_4M.root",
    const char* treename="DecayTree",
    const char* outputname="/mnt/e/sphenix/HF-QQ/output/Conv_reco_4M_ana.root"
    // const char* outputname="/mnt/e/sphenix/HF-QQ/output/Conv_bkgreco_4M_ana.root"
    // const char* outputname="/sphenix/user/jzhang1/PhysicsAna/HF-QQ/analysis/output/Conv_bkgreco_4M_ana.root"
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
    float mother_mass, mother_x, mother_y, mother_r;

    tree_data->SetBranchAddress("gamma_mass", & mother_mass);
    tree_data->SetBranchAddress("gamma_x", & mother_x);
    tree_data->SetBranchAddress("gamma_y", & mother_y);

    mother_r = std::sqrt(mother_x*mother_x + mother_y*mother_y);

    Long64_t nentries = tree_data->GetEntries();
    for (Long64_t i=0; i<nentries; ++i) 
    {
        tree_data->GetEntry(i);

        if(mother_r > 30) continue; // fiducial cut
        if(mother_r < 1)

        h1_mass->Fill(mother_mass);
    }

    outfile.cd();
    h1_mass->Write();

    std::cout<<"Ana Finished!"<<std::endl;
}


