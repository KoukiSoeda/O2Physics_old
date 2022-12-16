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

struct MFTTrackInfo {
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
      registry.add("hTrackType", "hTrackType", {HistType::kTH1F, {trackTypeAxis}});
      registry.add("TrackEta", "; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis,{18,-4.6, -1.}}});

   }

   void process(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const& particleMC, aod::McCollisions const& mcCol, soa::Join<aod::FwdTracks, aod::FwdTracksDCA> const& fwdtracks)
   {
      //if(collision.has_collision()){
         auto pt = 0.;
         auto dcax = 0.;
         auto dcay = 0.;
         auto dca = 0.;
         auto chargeSign = 0.;
         auto zvtx = 0.;
         auto nFwdTracks = 0.;
         zvtx = collision.posZ();

         if(!useEvSel || (useEvSel && collision.sel8())){
            for(auto const& fwdtrack : fwdtracks){
               registry.fill(HIST("hTrackType"), fwdtrack.trackType());
               if(fwdtrack.has_collision()){
                  if(fwdtrack.trackType() == 0){
                     nFwdTracks++;
                     auto Eta = fwdtrack.eta();
                     pt = fwdtrack.pt();
                     dcax = fwdtrack.fwdDcaX();
                     dcay = fwdtrack.fwdDcaY();
                     dca = sqrt(dcax * dcax + dcay * dcay);
                     chargeSign = fwdtrack.sign();
                     registry.fill(HIST("TrackEta"), Eta);
                     registry.fill(HIST("hBasicDist"), pt, dcax, dcay, dca, zvtx);
                  }
               }
            }
         }
      //}
   }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexdistribution>(cfgc),
                      adaptAnalysisTask<MFTTrackInfo>(cfgc),
  };
}