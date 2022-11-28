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
   Configurable<int> nBins{"nBins", 1000, "N bins in all histos"};
   //Histogram defined w/ HistogramRegistry
   HistogramRegistry registry{
      "registry",
      {
         {"nParticle", "nParticle", {HistType::kTH1F, {{nBins, 0.0, 1000}}}},
         {"MuonPtEta", "; p_{T} (GeV/c); #eta; tracks", {HistType::kTH2F, {{600, 0 ,2*M_PI}, {35, -4.5, -1}}}},
         {"MuonMotherParticle", "; pdgCode; counts", {HistType::kTH1F, {{nBins, 0.0, 1000}}}}
      }
   };

   void process(MFTTracksLabeled const& mfttracks, aod::McParticles const& particles, aod::McTrackLabels const& labels){
      for(auto& mfttrack : mfttracks){
         auto MFTTrackId = mfttrack.mcParticle();
         auto CollisionId = mfttrack.collisionId();
         auto PDGCode = MFTTrackId.pdgCode();
         registry.get<TH1>(HIST("nParticle"))->Fill(abs(PDGCode));
         if(PDGCode==13){
            registry.get<TH2>(HIST("MuonPtEta"))->Fill(MFTTrackId.pt(), MFTTrackId.eta());
         }
         
         if(MFTTrackId.has_mothers() && PDGCode==13){
            auto MotherParticles = MFTTrackId.mothersIds();
            auto motherId = particles.rawIteratorAt(MotherParticles[0]);
            auto MotherPDGCode = motherId.pdgCode();
            registry.get<TH1>(HIST("MuonMotherParticle"))->Fill(MotherPDGCode);
         }

      }
   }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
   return WorkflowSpec{
      adaptAnalysisTask<pca_mc>(cfgc)
   };
}