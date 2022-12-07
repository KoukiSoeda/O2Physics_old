/// \author Koki Soeda
/// \since 25/11/2022

#include <iostream>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using MFTTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

struct pca_mc {
   //Number of bin
   Configurable<int> nBins{"nBins", 2000, "N bins in all histos"};
   //Histogram defined w/ HistogramRegistry
   HistogramRegistry registry{
      "registry",
      {
         {"nParticle", "nParticle", {HistType::kTH1F, {{nBins, -1000, 1000}}}},
         {"muon_vx", "; x (cm); counts", {HistType::kTH1F, {{200, -10, 10}}}},
         {"muon_vy", "; y (cm); counts", {HistType::kTH1F, {{200, -10, 10}}}},
         {"muon_vz", "; z (cm); counts", {HistType::kTH1F, {{500, -100, 25}}}},
         {"muon_pT", "; p_{T} (GeV/c); ", {HistType::kTH1F, {{100, 0.0, 5}}}},
         {"muon_eta", "; #eta; ", {HistType::kTH1F, {{35, -4.5, -1}}}},
         {"MuonMotherParticle", "; pdgCode; counts", {HistType::kTH1F, {{nBins, -1000, 1000}}}},
         {"MuonMotherParticle_vx", "; x (cm); counts", {HistType::kTH1F, {{200, -10, 10}}}},
         {"MuonMotherParticle_vy", "; y (cm); counts", {HistType::kTH1F, {{200, -10, 10}}}},
         {"MuonMotherParticle_vz", "; z (cm); counts", {HistType::kTH1F, {{500, -25, 25}}}}
      }
   };

   void process(MFTTracksLabeled const& mfttracks, aod::McParticles const& particles, aod::McTrackLabels const& labels){
      for(auto& mfttrack : mfttracks){
         //Load MFT tracks information
         auto MFTTrackId = mfttrack.mcParticle();
         auto CollisionId = mfttrack.collisionId();
         auto PDGCode = MFTTrackId.pdgCode();
         registry.get<TH1>(HIST("nParticle"))->Fill(PDGCode);

         //Get muon MC information
         if(abs(PDGCode)==13){
            auto mu_px = MFTTrackId.px();
            auto mu_py = MFTTrackId.py();
            auto mu_pz = MFTTrackId.pz();
            auto mu_vx = MFTTrackId.vx();
            auto mu_vy = MFTTrackId.vy();
            auto mu_vz = MFTTrackId.vz();
            auto mu_phi = MFTTrackId.phi();
            auto mu_eta = MFTTrackId.eta();
            auto mu_pt = MFTTrackId.pt();
            registry.get<TH1>(HIST("muon_vx"))->Fill(mu_vx);
            registry.get<TH1>(HIST("muon_vy"))->Fill(mu_vy);
            registry.get<TH1>(HIST("muon_vz"))->Fill(mu_vz);
            registry.get<TH1>(HIST("muon_pT"))->Fill(mu_pt);
            registry.get<TH1>(HIST("muon_eta"))->Fill(mu_eta);
         }
         
         //Get muon mother particle information
         if(MFTTrackId.has_mothers() && abs(PDGCode)==13){
            auto MotherParticles = MFTTrackId.mothersIds();
            auto mu_motherId = particles.rawIteratorAt(MotherParticles[0]);
            auto mu_MotherPDGCode = mu_motherId.pdgCode();
            auto mu_mother_vx = mu_motherId.vx();
            auto mu_mother_vy = mu_motherId.vy();
            auto mu_mother_vz = mu_motherId.vz();
            registry.get<TH1>(HIST("MuonMotherParticle"))->Fill(mu_MotherPDGCode);
            registry.get<TH1>(HIST("MuonMotherParticle_vx"))->Fill(mu_mother_vx);
            registry.get<TH1>(HIST("MuonMotherParticle_vy"))->Fill(mu_mother_vy);
            registry.get<TH1>(HIST("MuonMotherParticle_vz"))->Fill(mu_mother_vz);

            auto DaughterParticles = mu_motherId.daughtersIds();
            auto daughter_1 = particles.rawIteratorAt(DaughterParticles[0]);
            auto daughter_2 = particles.rawIteratorAt(DaughterParticles[1]);
            auto daughter_3 = particles.rawIteratorAt(DaughterParticles[2]);
            
            auto daughter_1_PDGCode = daughter_1.pdgCode();
            auto daughter_1_px = daughter_1.px();
            auto daughter_1_py = daughter_1.py();
            auto daughter_1_pz = daughter_1.pz();
            auto daughter_1_vx = daughter_1.vx();
            auto daughter_1_vy = daughter_1.vy();
            auto daughter_1_vz = daughter_1.vz();
            auto daughter_2_PDGCode = daughter_2.pdgCode();
            auto daughter_2_px = daughter_2.px();
            auto daughter_2_py = daughter_2.py();
            auto daughter_2_pz = daughter_2.pz();
            auto daughter_2_vx = daughter_2.vx();
            auto daughter_2_vy = daughter_2.vy();
            auto daughter_2_vz = daughter_2.vz();
            auto daughter_3_PDGCode = daughter_3.pdgCode();
            auto daughter_3_px = daughter_3.px();
            auto daughter_3_py = daughter_3.py();
            auto daughter_3_pz = daughter_3.pz();
            auto daughter_3_vx = daughter_3.vx();
            auto daughter_3_vy = daughter_3.vy();
            auto daughter_3_vz = daughter_3.vz();

            
            if(abs(mu_MotherPDGCode)!=211) cout << mu_MotherPDGCode << ",  " << daughter_1_PDGCode << ",  "<< daughter_2_PDGCode << ",  " << daughter_3_PDGCode << endl;
            //cout << daughter_1_PDGCode << ": " << daughter_1_vz << ",  "<< daughter_2_PDGCode << ": " << daughter_2_vz << ",  " << daughter_3_PDGCode << ": " << daughter_3_vz << ",  " << daughter_4_PDGCode << ":" << daughter_4_vz << endl;
         

         }

      }
   }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
   return WorkflowSpec{
      adaptAnalysisTask<pca_mc>(cfgc)
   };
}