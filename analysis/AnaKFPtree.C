#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <vector>
#include <cmath>

void AnaKFPtree(
    // const char* filename="/mnt/e/sphenix/HF-QQ/output/Conv_reco_4M_match.root",
    // const char* filename="/mnt/e/sphenix/HF-QQ/output/Conv_bkgreco_4M_match.root",
    const char* filename="/sphenix/user/jzhang1/PhysicsAna/HF-QQ/output/PhotonConv/PhotonConv_reco/Reconstructed/53046/Conv_reco_4M.root",
    const char* treename="DecayTree",
    // const char* outputname="/mnt/e/sphenix/HF-QQ/output/Conv_reco_4M_match_ana.root"
    // const char* outputname="/mnt/e/sphenix/HF-QQ/output/Conv_bkgreco_4M_match_ana.root"
    const char* outputname="/sphenix/user/jzhang1/PhysicsAna/HF-QQ/analysis/output/Conv_reco_4M_match_ana.root"
) 
{
    TFile outfile(outputname, "RECREATE");
    TH1D *h1_mass = new TH1D("h1_mass", "h1_mass; mass [GeV]; count", 500, 0, 0.5);
    TH1D *h1_Rxy = new TH1D("h1_Rxy", "h1_Rxy; Rxy [cm]; count", 300, 0, 30);
    
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
    float track_1_px, track_1_py, track_1_pz;
    float track_2_px, track_2_py, track_2_pz;

    tree_data->SetBranchAddress("gamma_mass", & mother_mass);
    tree_data->SetBranchAddress("gamma_x", & mother_x);
    tree_data->SetBranchAddress("gamma_y", & mother_y);

    tree_data->SetBranchAddress("track_1_px", & track_1_px);
    tree_data->SetBranchAddress("track_1_py", & track_1_py);
    tree_data->SetBranchAddress("track_1_pz", & track_1_pz);

    tree_data->SetBranchAddress("track_2_px", & track_2_px);
    tree_data->SetBranchAddress("track_2_py", & track_2_py);
    tree_data->SetBranchAddress("track_2_pz", & track_2_pz);

    Long64_t nentries = tree_data->GetEntries();
    for (Long64_t i=0; i<nentries; ++i) 
    {
        tree_data->GetEntry(i);

        mother_r = std::sqrt(mother_x*mother_x + mother_y*mother_y);

        double mpion = 0.13957; // GeV
        double p1 = sqrt(track_1_px*track_1_px + track_1_py*track_1_py + track_1_pz*track_1_pz);
        double p2 = sqrt(track_2_px*track_2_px + track_2_py*track_2_py + track_2_pz*track_2_pz);
        double Ee1 = sqrt(p1*p1 + mpion*mpion);
        double Ee2 = sqrt(p2*p2 + mpion*mpion);
        double m_recalc = sqrt(pow(Ee1+Ee2,2) - pow(track_1_px+track_2_px,2)
                               - pow(track_1_py+track_2_py,2) - pow(track_1_pz+track_2_pz,2));

        // h1_Rxy->Fill(mother_r);

        // if(mother_r > 15) continue; 
        if( ((mother_r > 2.5)&&(mother_r<4.5)) || ((mother_r>7)&&(mother_r<12)) )
        // if( (mother_r < 0.5)||(mother_r>20) )
        {
            h1_Rxy->Fill(mother_r);
            h1_mass->Fill(mother_mass);
        }
    }

    outfile.cd();
    h1_mass->Write();
    h1_Rxy->Write();

    std::cout<<"Ana Finished!"<<std::endl;
}


