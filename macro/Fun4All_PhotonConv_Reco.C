#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Calo_Calib.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>
#include <phool/getClass.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/TpcClusterQA.h>
#include <tpcqa/TpcRawHitQA.h>
#include <trackingqa/TpcSeedsQA.h>

#include <trackreco/PHActsTrackProjection.h>

#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>
#include <trackingdiagnostics/KshortReconstruction.h>

#include <trackbase_historic/SvtxTrackMap.h>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterBuilderTopo.h>

#include <calotrigger/TriggerRunInfoReco.h>

#include <stdio.h>

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wundefined-internal"

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#pragma GCC diagnostic pop

void KFPReco(std::string module_name = "KFPReco", std::string decaydescriptor = "K_S0 -> pi^+ pi^-", std::string outfile = "KFP.root", std::string trackmapName = "SvtxTrackMap", std::string containerName = "KFParticle");

float new_cemc_rad = 93.5;

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libepd.so)
R__LOAD_LIBRARY(libzdcinfo.so)

// trkr506:calo468 - 5:1
// trkr506:calo509 - 10:1
// 预筛 KFP：用 pt / chi2/ndf / ndf / crossing 控制规模（不依赖子探测器命中接口）
class PreKFPFilter : public SubsysReco {
 public:
    PreKFPFilter(double min_pt=0.5, double max_chi2ndf=30.0, unsigned min_ndf=30,
                 bool bunch0_only=true, unsigned long long max_pairs=200000ULL)
    : SubsysReco("PreKFPFilter"),
        m_minPt(min_pt), m_maxChi2NDF(max_chi2ndf), m_minNDF(min_ndf),
        m_bunch0(bunch0_only), m_maxPairs(max_pairs) {}

    int process_event(PHCompositeNode* topNode) override {
        auto trkmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
        if (!trkmap) return Fun4AllReturnCodes::EVENT_OK;

        size_t n_pos=0, n_neg=0, n_good=0;
        for (auto it = trkmap->begin(); it != trkmap->end(); ++it) {
            const SvtxTrack* t = it->second;
            if (!t) continue;

            // 1) 基本质量：pt / chi2/ndf / ndf(与命中数相关)
            const unsigned ndf = std::max(1u, t->get_ndf());
            const double chi2ndf = t->get_chisq() / ndf;
            if (t->get_pt() < m_minPt) continue;
            if (chi2ndf > m_maxChi2NDF) continue;
            if (ndf < m_minNDF) continue;

            // 2) 束团选择（如果 API 不支持 get_crossing，可把这几行注释掉）
            if (m_bunch0) {
            // 新接口通常叫 get_crossing()
            if (t->get_crossing() != 0) continue;
            }

            ++n_good;
            (t->get_charge() > 0 ? ++n_pos : ++n_neg);
        }

        // K_S^0 异号对上界；保守一点可改用 n_good*(n_good-1)/2
        const unsigned long long pairs = 1ULL * n_pos * n_neg;

        if (Verbosity() > 0) {
            auto rc = recoConsts::instance();
            std::cout << "[PreKFPFilter] run=" << rc->get_IntFlag("RUNNUMBER")
                    << " nTracks=" << trkmap->size()
                    << " n_good=" << n_good
                    << " (pos=" << n_pos << ", neg=" << n_neg << ")"
                    << " pairs_upper=" << pairs << std::endl;
        }

        if (pairs > m_maxPairs) {
            if (Verbosity() > 0)
            std::cout << "[PreKFPFilter] discard: pairs " << pairs
                        << " > " << m_maxPairs << std::endl;
            return Fun4AllReturnCodes::DISCARDEVENT;
        }
        return Fun4AllReturnCodes::EVENT_OK;
    }

    private:
    double m_minPt, m_maxChi2NDF;
    unsigned m_minNDF;
    bool m_bunch0;
    unsigned long long m_maxPairs;
};

void Fun4All_PhotonConv_Reco(
    const int nEvents = 10,
    const std::string trkr_clusterfilename = "DST_TRKR_CLUSTER_run2pp_ana505_2024p023_v001-00053046-00002.root",
    const std::string trkr_clusterdir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana505_2024p023_v001/DST_TRKR_CLUSTER/run_00053000_00053100/dst/",
    const std::string trkr_trackfilename = "DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053046-00002.root",
    const std::string trkr_trackdir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053000_00053100/dst/",
    const std::string calofilename = "DST_CALO_run2pp_ana509_2024p022_v001-00053046-00000.root",
    const std::string calodir = "/sphenix/lustre01/sphnxpro/production2/run2pp/physics/calocalib/ana509_2024p022_v001/run_00053000_00053100/",
    const std::string kfp_header = "outputKFParticle_",
    const std::string outdir = "./root",
    const int index = 0,
    const int stepsize = 10)
{
    // Input files names
    std::string inputTrkrTrackFile = trkr_trackdir + trkr_trackfilename;
    std::string inputCaloFile = calodir + calofilename;
    std::string inputTrkrClusterFile = trkr_clusterdir + trkr_clusterfilename;

    // get and set run number and segment from input file name
    std::pair<int, int>
        runseg = Fun4AllUtils::GetRunSegment(trkr_trackfilename);
    int runnumber = runseg.first;
    int segment = runseg.second;

    std::stringstream nice_runnumber;
    nice_runnumber << std::setw(8) << std::setfill('0') << to_string(runnumber);

    std::stringstream nice_segment;
    nice_segment << std::setw(5) << std::setfill('0') << to_string(segment);

    std::stringstream nice_index;
    nice_index << std::setw(5) << std::setfill('0') << to_string(index);

    std::stringstream nice_stepsize;
    nice_stepsize << std::setw(5) << std::setfill('0') << to_string(stepsize);

    std::cout << " run: " << runnumber
              << " samples: " << TRACKING::reco_tpc_maxtime_sample
              << " pre: " << TRACKING::reco_tpc_time_presample
              << " vdrift: " << G4TPC::tpc_drift_velocity_reco
              << std::endl;
    
    // distortion calibration mode
    /*
     * set to true to enable residuals in the TPC with
     * TPC clusters not participating to the ACTS track fit
     */
    G4TRACKING::SC_CALIBMODE = false;
    Enable::MVTX_APPLYMISALIGNMENT = true;
    ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
    TRACKING::pp_mode = true;

    auto se = Fun4AllServer::instance();
    se->Verbosity(1);

    auto rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", runnumber);
    rc->set_IntFlag("RUNSEGMENT", segment);

    // CDB and some fundamental setting need ask
    Enable::CDB = true;
    // rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag"); // whats this setting? "ProdA_2024"? "newcdbtag"?
    // rc->set_StringFlag("CDB_GLOBALTAG", CDB::global_tag);
    rc->set_StringFlag("CDB_GLOBALTAG", "2024p023");
    rc->set_uint64Flag("TIMESTAMP", runnumber);
    std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

    Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
    ingeo->AddFile(geofile);
    se->registerInputManager(ingeo);

    TpcReadoutInit( runnumber );
    // these lines show how to override the drift velocity and time offset values set in TpcReadoutInit
    G4TPC::tpc_drift_velocity_reco = 0.0073844; // cm/ns
    TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
    G4TPC::tpc_tzero_reco = -5*50;  // ns

    G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;

    // to turn on the default static corrections, enable the two lines below
    G4TPC::ENABLE_STATIC_CORRECTIONS = true;
    G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

    //to turn on the average corrections, enable the three lines below
    //note: these are designed to be used only if static corrections are also applied
    G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
    G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
    // to use a custom file instead of the database file:
    G4TPC::average_correction_filename = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
    
    G4MAGNET::magfield_rescale = 1;
    TrackingInit();

    auto hitsin_track = new Fun4AllDstInputManager("DSTin_track");
    hitsin_track->fileopen(inputTrkrTrackFile);
    se->registerInputManager(hitsin_track);

    auto hitsin_calo = new Fun4AllDstInputManager("DSTin_calo");
    hitsin_calo->fileopen(inputCaloFile);
    se->registerInputManager(hitsin_calo);

    auto hitsin_cluster = new Fun4AllDstInputManager("DSTin_cluster");
    hitsin_cluster->fileopen(inputTrkrClusterFile);
    se->registerInputManager(hitsin_cluster);

    /////////////////////////////////////////////////////
    // Set status of CALO towers, Calibrate towers, Cluster
    // Process_Calo_Calib(); // ???

    Global_Reco();

    // project track to calorimeter radius
    bool doEMcalRadiusCorr = true;
    auto projection = new PHActsTrackProjection("CaloProjection");
    new_cemc_rad = 102.9; // Virgile recommendation according to DetailedCalorimeterGeometry
    if (doEMcalRadiusCorr)
    {
        projection->setLayerRadius(SvtxTrack::CEMC, new_cemc_rad);
    }
    se->registerSubsystem(projection);

    // my calo reco
    // Calo Cluster reco with exact tower geom. setting
    std::cout<<"Begin my calo reco"<<std::endl;
    // Load the modified geometry
    CaloGeomMapping *cgm = new CaloGeomMapping();
    cgm->set_detector_name("CEMC");
    cgm->set_UseDetailedGeometry(true);
    se->registerSubsystem(cgm);

    //////////////////
    // EMC Clusters reco with exact tower geom. setting
    std::cout << "Building clusters" << std::endl;
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->set_UseDetailedGeometry(true);
    ClusterBuilder->set_threshold_energy(0.070);  // for when using basic calibration
    std::string emc_prof = getenv("CALIBRATIONROOT");
    emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
    //ClusterBuilder->set_UseAltZVertex(3); //0: GlobalVertexMap, 1: MbdVertexMap, 2: Nothing, 3: G4TruthInfo
    ClusterBuilder->LoadProfile(emc_prof);
    ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
    se->registerSubsystem(ClusterBuilder);

    //For particle flow studies
    RawClusterBuilderTopo* ClusterBuilder2 = new RawClusterBuilderTopo("EMcalRawClusterBuilderTopo2");
    ClusterBuilder2->Verbosity(0);
    ClusterBuilder2->set_nodename("TOPOCLUSTER_HCAL");
    ClusterBuilder2->set_enable_HCal(true);
    ClusterBuilder2->set_enable_EMCal(false);
    //ClusterBuilder2->set_noise(0.0025, 0.006, 0.03);
    ClusterBuilder2->set_noise(0.01, 0.03, 0.03);
    ClusterBuilder2->set_significance(4.0, 2.0, 1.0);
    ClusterBuilder2->allow_corner_neighbor(true);
    ClusterBuilder2->set_do_split(true);
    ClusterBuilder2->set_minE_local_max(1.0, 2.0, 0.5);
    ClusterBuilder2->set_R_shower(0.025);
    //se->registerSubsystem(ClusterBuilder2);

    // TriggerRunInfoReco ???
    TriggerRunInfoReco *triggerruninforeco = new TriggerRunInfoReco();
    se->registerSubsystem(triggerruninforeco);

    // // check track multiplicity, skip event if too large
    auto mon = new PreKFPFilter(
        /*min_pt=*/0.0,
        /*max_chi2ndf=*/1e9,
        /*min_ndf=*/0,
        /*bunch0_only=*/false,     // 先不要限定 crossing
        /*max_pairs=*/ULLONG_MAX   // 永不触发丢事件
    );
    mon->Verbosity(1);           // 设为 1 或 2，多点日志
    // se->registerSubsystem(mon);

    // output directory and file name setting
    string trailer = "_" + nice_runnumber.str() + "_" + nice_segment.str() + "_" + nice_index.str() + ".root";

    std::string PhotonConv_reco_dir = outdir + "/PhotonConv_reco/inReconstruction/" + to_string(runnumber) + "/";
    std::string PhotonConv_reco_filename = kfp_header + "PhotonConv_reco" + trailer;
    TString PhotonConv_outfile = PhotonConv_reco_dir + PhotonConv_reco_filename;
    std::cout<<"outfile "<<PhotonConv_outfile<<std::endl;
    {
        std::string makeDirectory = "mkdir -p " + PhotonConv_reco_dir;
        system(makeDirectory.c_str());
    }
    std::string PhotonConv_outfile_string(PhotonConv_outfile.Data());
    KFPReco("PhotonConvReco", "[gamma -> e^+ e^-]cc", PhotonConv_outfile_string);
    // KFPReco("PhotonConvReco", "[K_S0 -> pi^+ pi^-]cc", PhotonConv_outfile_string);

    std::string PhotonConv_reco_likesign_dir = outdir + "/PhotonConv_reco_likesign/inReconstruction/" + to_string  (runnumber) + "/";
    std::string PhotonConv_reco_likesign_filename = kfp_header + "PhotonConv_reco_likesign" + trailer;
    TString PhotonConv_likesign_outfile = PhotonConv_reco_likesign_dir + PhotonConv_reco_likesign_filename;
    std::cout<<"outfile "<<PhotonConv_likesign_outfile<<std::endl;
    {
        std::string makeDirectory = "mkdir -p " + PhotonConv_reco_likesign_dir;
        system(makeDirectory.c_str());
    }
    std::string PhotonConv_outfile_likesign_string(PhotonConv_likesign_outfile.Data());
    KFPReco("PhotonConvReco_likesign", "[gamma -> e^+ e^+]cc", PhotonConv_outfile_likesign_string);
    // KFPReco("PhotonConvReco_likesign", "[K_S0 -> pi^+ pi^+]cc", PhotonConv_outfile_likesign_string);

    se->skip(stepsize*index);
    se->run(nEvents);
    se->End();
    se->PrintTimer();
    CDBInterface::instance()->Print();

    ifstream file_PhotonConv_outfile_string(PhotonConv_outfile_string.c_str(), ios::binary | ios::ate);
    if (file_PhotonConv_outfile_string.good() && (file_PhotonConv_outfile_string.tellg() > 100))
    {
        string outputDirMove = outdir + "/PhotonConv_reco/Reconstructed/" + to_string(runnumber) + "/";
        string makeDirectoryMove = "mkdir -p " + outputDirMove;
        system(makeDirectoryMove.c_str());
        string moveOutput = "mv " + PhotonConv_outfile_string + " " + outputDirMove;
        std::cout << "moveOutput: " << moveOutput << std::endl;
        system(moveOutput.c_str());
    }

    ifstream file_PhotonConv_outfile_likesign_string(PhotonConv_outfile_likesign_string.c_str(), ios::binary | ios::ate);
    if (file_PhotonConv_outfile_likesign_string.good() && (file_PhotonConv_outfile_likesign_string.tellg() > 100))
    {
        string outputDirMove = outdir + "/PhotonConv_reco_likesign/Reconstructed/" + to_string(runnumber) + "/";
        string makeDirectoryMove = "mkdir -p " + outputDirMove;
        system(makeDirectoryMove.c_str());
        string moveOutput = "mv " + PhotonConv_outfile_likesign_string + " " + outputDirMove;
        std::cout << "moveOutput: " << moveOutput << std::endl;
        system(moveOutput.c_str());
    }

    delete se;
    std::cout << "Finished" << std::endl;
    gSystem->Exit(0);
}

// Photon conversion reconstruction with KFParticle
void KFPReco(std::string module_name = "KFPReco", std::string decaydescriptor = "K_S0 -> pi^+ pi^-", std::string outfile = "KFP.root", std::string trackmapName = "SvtxTrackMap", std::string containerName = "KFParticle")
{
    auto se = Fun4AllServer::instance();
    //KFParticle setup
    KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(module_name);
    kfparticle->Verbosity(0);
    kfparticle->setDecayDescriptor(decaydescriptor);

    kfparticle->setTrackMapNodeName(trackmapName);
    kfparticle->setContainerName(containerName);

    //Basic node selection and configuration
    kfparticle->dontUseGlobalVertex(false);
    kfparticle->requireTrackVertexBunchCrossingMatch(false);
    kfparticle->getAllPVInfo(false);
    kfparticle->allowZeroMassTracks(true);
    kfparticle->use2Dmatching(false);
    kfparticle->getTriggerInfo(true);
    kfparticle->getDetectorInfo(true);
    kfparticle->GetDetailedTracking(false);
    kfparticle->getCaloInfo(true);
    kfparticle->useFakePrimaryVertex(false);
    kfparticle->saveDST(false);
    kfparticle->saveParticleContainer(false);
    kfparticle->saveTrackContainer(false);
    kfparticle->magFieldFile("FIELDMAP_TRACKING");

    //PV to SV cuts
    kfparticle->constrainToPrimaryVertex(true);
    kfparticle->setMotherIPchi2(FLT_MAX);
    kfparticle->setFlightDistancechi2(-1.);
    kfparticle->setMinDIRA(-1.1);
    kfparticle->setDecayLengthRange(0., FLT_MAX);
    kfparticle->setDecayLengthRange_XY(-10000, FLT_MAX);
    kfparticle->setDecayTimeRange_XY(-10000, FLT_MAX);
    kfparticle->setDecayTimeRange(-10000, FLT_MAX);
    kfparticle->setMinDecayTimeSignificance(-1e5);
    kfparticle->setMinDecayLengthSignificance(-1e5);
    kfparticle->setMinDecayLengthSignificance_XY(-1e5);
    kfparticle->setMaximumDaughterDCA_XY(FLT_MAX);
    kfparticle->setMaximumDaughterDCA(FLT_MAX);

    //Track parameters
    kfparticle->bunchCrossingZeroOnly(true);
    kfparticle->setMinMVTXhits(0);
    kfparticle->setMinINTThits(0);
    kfparticle->setMinTPChits(35);
    kfparticle->setMinimumTrackPT(0.5);
    kfparticle->setMaximumTrackPTchi2(FLT_MAX);
    kfparticle->setMinimumTrackIPchi2(-1.);
    kfparticle->setMinimumTrackIP(-1.);
    kfparticle->setMaximumTrackchi2nDOF(100.);

    // //Track-Calo matching
    // kfparticle->set_emcal_radius_user(new_cemc_rad);
    // //narrow window
    // kfparticle->set_dphi_cut_low(-0.02); //rad
    // kfparticle->set_dphi_cut_high(0.09); //rad
    // kfparticle->set_dz_cut_low(-4); //cm
    // kfparticle->set_dz_cut_high(4); //cm
    // //loose window
    // /*
    // kfparticle->set_dphi_cut_low(-0.2); //rad
    // kfparticle->set_dphi_cut_high(0.2); //rad
    // kfparticle->set_dz_cut_low(-10); //cm
    // kfparticle->set_dz_cut_high(10); //cm
    // */
    // kfparticle->set_emcal_e_low_cut(0.2); //GeV
    // kfparticle->requireTrackEMCalMatch(true);

    //Vertex parameters
    kfparticle->setMaximumVertexchi2nDOF(FLT_MAX);

    //Parent parameters
    kfparticle->setMotherPT(0);
    kfparticle->setMinimumMass(-1);
    kfparticle->setMaximumMass(15);
    kfparticle->setMaximumMotherVertexVolume(FLT_MAX);

    kfparticle->setOutputName(outfile);

    se->registerSubsystem(kfparticle);
}

// // KFP K reco to ensure the code is working
// void KFPReco(std::string module_name = "KFPReco", std::string decaydescriptor = "K_S0 -> pi^+ pi^-", std::string outfile = "KFP.root", std::string trackmapName = "SvtxTrackMap", std::string containerName = "KFParticle")
// {   
//     bool save_tracks_to_DST = false;
//     bool dont_use_global_vertex = true;
//     bool require_track_and_vertex_match = true;
//     bool save_all_vtx_info = true;
//     bool extrapolate_track_to_SV = false;
//     bool constrain_phi_mass = false;
//     bool constrain_lambda_mass = false;
//     bool constrain_D_mass = false;
//     bool use_2D_matching = false;
//     bool get_trigger_info = true;
//     bool get_detector_info = true;
//     bool get_detailed_tracking = false;
//     bool get_dEdx_info = true;
//     bool constrain_to_primary_vertex = true;
//     bool use_pid = true;
//     bool use_calo = false;
//     float pid_frac = 0.4;

//     Fun4AllServer *se = Fun4AllServer::instance();

//     KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(module_name);
//     //KFParticle setup
//     kfparticle->Verbosity(0);
//     kfparticle->setDecayDescriptor(decaydescriptor);

//     kfparticle->usePID(use_pid);
//     kfparticle->setPIDacceptFraction(pid_frac);
//     kfparticle->dontUseGlobalVertex(dont_use_global_vertex);
//     kfparticle->requireTrackVertexBunchCrossingMatch(require_track_and_vertex_match);
//     kfparticle->extraolateTracksToSV(extrapolate_track_to_SV);
//     kfparticle->getAllPVInfo(save_all_vtx_info);
//     kfparticle->allowZeroMassTracks();
//     kfparticle->use2Dmatching(use_2D_matching);
//     kfparticle->getTriggerInfo(get_trigger_info);
//     kfparticle->getDetectorInfo(get_detector_info);
//     kfparticle->GetDetailedTracking(get_detailed_tracking);
//     kfparticle->saveDST(save_tracks_to_DST);
//     kfparticle->saveParticleContainer(false);
//     kfparticle->magFieldFile("FIELDMAP_TRACKING");

//     //PV to SV cuts
//     kfparticle->constrainToPrimaryVertex(constrain_to_primary_vertex);
//     //kfparticle->setMotherIPchi2(100);//for default navigator
//     kfparticle->setMotherIPchi2(FLT_MAX);//for direct navigator
//     kfparticle->setFlightDistancechi2(-1.);
//     kfparticle->setMinDIRA(0.88); //was .95
//     kfparticle->setDecayLengthRange(.2, FLT_MAX); //was 0.1 min

//     kfparticle->setDecayLengthRange_XY(-10000, FLT_MAX);
//     kfparticle->setDecayTimeRange_XY(-10000, FLT_MAX);
//     kfparticle->setDecayTimeRange(-10000, FLT_MAX);
//     kfparticle->setMinDecayTimeSignificance(-1e5);
//     kfparticle->setMinDecayLengthSignificance(-1e5);
//     kfparticle->setMinDecayLengthSignificance_XY(-1e5);
//     kfparticle->setMaximumDaughterDCA_XY(100);

//     //Track parameters
//     kfparticle->setMinimumTrackPT(0.0);
//     kfparticle->setMinimumTrackIPchi2(-1.);
//     kfparticle->setMinimumTrackIP(-1.);
//     kfparticle->setMaximumTrackchi2nDOF(100.);
//     kfparticle->setMinINTThits(1); //was 0
//     kfparticle->setMinMVTXhits(1); //was 0
//     kfparticle->setMinTPChits(20); //was 20

//     //Vertex parameters
//     kfparticle->setMaximumVertexchi2nDOF(20);
//     kfparticle->setMaximumDaughterDCA(0.1); //5 mm //was 0.5

//     //Parent parameters
//     kfparticle->setMotherPT(0);
//     kfparticle->setMinimumMass(0.40); //Check mass ranges
//     kfparticle->setMaximumMass(0.60);
//     kfparticle->setMaximumMotherVertexVolume(0.1);

//     kfparticle->setOutputName(outfile);

//     se->registerSubsystem(kfparticle);
// }
