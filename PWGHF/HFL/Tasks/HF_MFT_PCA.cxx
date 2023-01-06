/// \author Koki Soeda
/// \since 16/12/2022

#include <iostream>
#include <cmath>
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Common/Core/RecoDecay.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"


using namespace std;
using namespace o2; 
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct vertexdistribution {
  OutputObj<TH1F> vertex{TH1F("vertex", "vertex", 400, -20, 20)};

  Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};

  // loop over MC truth McCollisions
  void process(aod::McCollision const& mcCollision)
  {
    vertex->Fill(mcCollision.posZ());
  }
};

struct MftTrackInfo {
   Service<TDatabasePDG> pdg;
   Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};

   HistogramRegistry registry{"registry"};

   void init(InitContext&)
   {
      ConfigurableAxis PtAxis{"PtAxis", {1001, -0.05, 10.05}, "pt axis for histograms"};
      AxisSpec ptRecoAxis = {1500, 0, 15, "#it{p}_{T}_{Reco}"};
      AxisSpec trackTypeAxis = {6, -0.5, 5.5, "Track Type"};
      AxisSpec dcaxAxis = {1000, -5.0, 5.0, "DCA {x or y} (cm)"};
      AxisSpec dcaAxis = {1000, 0.0, 100.0, "DCA {xy} (cm)"};
      AxisSpec zvtxAxis = {400, -20.0, 20.0, "zvtx (cm)"};
      AxisSpec etaRecoAxis = {150, -5, -2, "#eta_{Reco}"};
      AxisSpec rAbsAxis = {100, 0, 100, "R_{abs}"};
      AxisSpec pdcaAxis = {450, 0, 450, "p_{DCA}"};
      AxisSpec chi2GlobalAxis = {170, -1.5, 150.5, "#chi^{2} global"};
      AxisSpec chi2MCHMFTAxis = {170, -1.5, 150.5, "#chi^{2} MCH-MFT"};
      AxisSpec chi2MCHMIDAxis = {170, -1.5, 150.5, "#chi^{2} MCH-MID"};
      HistogramConfigSpec HistVariable({HistType::kTHnSparseF, {ptRecoAxis, dcaxAxis, dcaxAxis, dcaAxis, zvtxAxis}});

      registry.add("hBasicDist", "", HistVariable);
      registry.add("hpTEta", "; pT (GeV/c); #eta", {HistType::kTH2F, {ptRecoAxis, etaRecoAxis}});

   }

   void process(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const& particleMC, aod::McCollisions const& mcCol)
   {
      for(auto& track : tracks){
         auto colId = track.collisionId();
         auto mftx = track.x();
         auto mfty = track.y();
         auto mftz = track.z();
         auto mftPhi = track.phi();
         auto mftEta = track.eta();
         auto mftpt = track.pt();
         registry.fill(HIST("hpTEta"), mftpt, mftEta);
         //if(mftpt<<0.5) continue;
         
      }
   }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexdistribution>(cfgc),
                      adaptAnalysisTask<MftTrackInfo>(cfgc),
  };
}