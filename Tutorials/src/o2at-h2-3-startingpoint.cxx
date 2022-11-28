#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


//STEP 0
//This is an example of a conveient declaration of "using"
//WARNING: THIS IS NEW
using MyCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

//Starting point
struct partandfiltexample {
  Filter etaFilter = nabs(aod::track::eta) < 0.5f;
  Filter trackDCA = nabs(aod::track::dcaXY) <0.2f;
  using MyFilteredTracks = soa::Filtered<MyCompleteTracks>;

  Partition<MyFilteredTracks> leftTracks = aod::track::eta < 0.0f;
  Partition<MyFilteredTracks> rightTracks = aod::track::eta >= 0.0f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHistogramleft", "etaHistogramleft", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramleft", "ptHistogramleft", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"etaHistogramright", "etaHistogramright", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramright", "ptHistogramright", {HistType::kTH1F, {{nBins, 0., 10.0}}}}
    }
  };

  void process(aod::Collision const& collision, MyFilteredTracks const& tracks)
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    
    //partitions are not grouped by default
    auto leftTracksGrouped = leftTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto rightTracksGrouped = rightTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());

    for (auto& track : leftTracksGrouped) { // only for a subset
      if(track.tpcNClsCrossedRows() < 70) continue;
      registry.get<TH1>(HIST("etaHistogramleft"))->Fill(track.eta());
      registry.get<TH1>(HIST("ptHistogramleft"))->Fill(track.pt());
    }
    for (auto&track : rightTracksGrouped) {
      if(track.tpcNClsCrossedRows() < 70) continue;
      registry.get<TH1>(HIST("etaHistogramright"))->Fill(track.eta());
      registry.get<TH1>(HIST("ptHistogramright"))->Fill(track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<partandfiltexample>(cfgc)
  };
}