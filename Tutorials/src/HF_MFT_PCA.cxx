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
   void process(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const& particleMC, aod::McCollisions const& mcCol)
   {
      if(collision.has_Collision()){
         if(!useEvSel || (useEvSel && collision.sel8())){
            for(auto& track : tracks){
               
            }
         }
      }
   }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexdistribution>(cfgc),
                      adaptAnalysisTask<MFTTrackInfo>(cfgc),
  };
}