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
         if(PDGCode==13){
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
         
         //Get mother particle information
         if(MFTTrackId.has_mothers() && PDGCode==13){
            auto MotherParticles = MFTTrackId.mothersIds();
            auto motherId = particles.rawIteratorAt(MotherParticles[0]);
            auto MotherPDGCode = motherId.pdgCode();
            auto mother_vx = motherId.vx();
            auto mother_vy = motherId.vy();
            auto mother_vz = motherId.vz();
            registry.get<TH1>(HIST("MuonMotherParticle"))->Fill(MotherPDGCode);
            registry.get<TH1>(HIST("MuonMotherParticle_vx"))->Fill(mother_vx);
            registry.get<TH1>(HIST("MuonMotherParticle_vy"))->Fill(mother_vy);
            registry.get<TH1>(HIST("MuonMotherParticle_vz"))->Fill(mother_vz);
         }

      }
   }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
   return WorkflowSpec{
      adaptAnalysisTask<pca_mc>(cfgc)
   };
}