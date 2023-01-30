/// \author Koki Soeda
/// \since 30.01.2023

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
using namespace o2::track;
using namespace o2::soa;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct DCAandPCA {
   Service<TDatabasePDG> pdg;
   Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
   
   HistogramRegistry registry{
      "registry",
      {
         {"MFTTrack_z", " ; z(cm); ", {HistType::kTH1F, {{20000, -100, 100}}}},
      }
   };

   void processMCDCA(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                     MFTTracksLabeled const& tracks,
                     aod::McParticles const& particleMC,
                     aod::McCollisions const&)
   {
      if(collision.has_mcCollision()){
         if((collision.mcCollision().posZ() < zMax) && (collision.mcCollision().posZ() > -zMax)){
            if(!useEvSel || (useEvSel && collision.sel8())){
               float dcaXY;
               o2::track::TrackParCovFwd trackPar;
               for(auto& track : tracks){
                  if(!track.has_mcParticle()) continue;
                  auto particle = track.mcParticle();
                  auto t_collision = track.collision();
                  
                  if(fabs(particle.pdgCode())!=13) continue;
                  if(particle.mcCollisionId() == collision.mcCollision().globalIndex()){
                     if(particle.has_mothers()){
                        auto mcMom = particleMC.rawIteratorAt(particle.mothersIds()[0]);
                        //if(fabs(mcMom.pdgCode()==421){

                        //}
                     }
                  }
               }
            }
         }
      }      
   }
   PROCESS_SWITCH(DCAandPCA, processMCDCA, "Get the DCAxy of MC information for MFT", true);
};

WorkflowSpec
   defineDataProcessing(ConfigContext const& cfgc)
{
   return WorkflowSpec{adaptAnalysisTask<DCAandPCA>(cfgc)};
}
