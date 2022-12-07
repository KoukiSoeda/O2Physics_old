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
         {"MuonMotherParticle_vz", "; z (cm); counts", {HistType::kTH1F, {{500, -25, 25}}}},
         {"RecoMass", "; (GeV/c2); counts", {HistType::kTH1F, {{1000, 1.7, 2.0}}}}
      }
   };

   void process(MFTTracksLabeled const& mfttracks, aod::McParticles const& particles, aod::McTrackLabels const& labels){
      int i = 0;
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

            if(abs(mu_MotherPDGCode)==421){
               i++;
               float Daughter_E = 0;
               float Daughter_px = 0;
               float Daughter_py = 0;
               float Daughter_pz = 0;
               auto Daughters = mu_motherId.daughters_as<aod::McParticles>();
               //int ddN = Daughters.size();
               for(auto& Daughter : Daughters){
                  Daughter_E += Daughter.e();
                  Daughter_px += Daughter.px();
                  Daughter_py += Daughter.py();
                  Daughter_pz += Daughter.pz();
               }
               auto iMass = sqrt(pow(Daughter_E,2)-pow(Daughter_px,2)-pow(Daughter_py,2)-pow(Daughter_pz,2));
               cout << iMass << ": " << i << endl;
               registry.get<TH1>(HIST("RecoMass"))->Fill(iMass);
            }
         }


      }
   }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
   return WorkflowSpec{
      adaptAnalysisTask<pca_mc>(cfgc)
   };
}