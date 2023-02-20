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
using namespace o2::soa;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct DCAandPCA {
   Service<TDatabasePDG> pdg;
   Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
   Configurable<float> zMax{"zMax", 5., "value for Zvtx cut"};
   
   HistogramRegistry registry{
      "registry",
      {
         {"muTrack_dcax", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"muTrack_dcay", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"muTrack_mftx", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"muTrack_mfty", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"D0decayx", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"D0decayy", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"D0decayz", " ; z(cm); ", {HistType::kTH1F, {{5000, -50, 50}}}},
         {"pcax", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"pcay", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"pcaz", " ; z(cm); ", {HistType::kTH1F, {{5000, -50, 50}}}},
         {"PIDpurity", "", {HistType::kTH1F, {{102, -0.5, 100.5}}}},
         {"EstimatePDGcode", "", {HistType::kTH1F, {{2002, -1000.5, 1000.5}}}},
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
               for(auto& track0 : tracks){
                  if(!track0.has_mcParticle()) continue;
                  auto particle0 = track0.mcParticle();
                  auto t_collision = track0.collision();
                  
                  if(fabs(particle0.pdgCode())!=13) continue;
                  int estpdg;
                  int64_t truthK_ID,muonID,estID;
                  float vx_mu,vy_mu,vz_mu,px_mu,py_mu,pz_mu,t_mu,mft_mu_x,mft_mu_y,s_mu,dca_mu_x,dca_mu_y;
                  //float mft_det_x, mft_det_y, mft_det_z;
                  float vx_can,vy_can,vz_can,px_can,py_can,pz_can,t_can,mft_can_x,mft_can_y,s_can,dca_can_x,dca_can_y;
                  float s,t,a,b,c,d,e,f,g,h;
                  float mu_x,mu_y,mu_z,can_x,can_y,can_z,pca_x,pca_y,pca_z;
                  if(particle0.mcCollisionId() == collision.mcCollision().globalIndex()){
                     if(particle0.has_mothers()){
                        int signMUON = 0;
                        float closest_track = 10000;
                        auto mcMom = particleMC.rawIteratorAt(particle0.mothersIds()[0]);
                        auto Daughters = mcMom.daughters_as<aod::McParticles>();
                        registry.fill(HIST("D0decayx"), particle0.vx());
                        registry.fill(HIST("D0decayy"), particle0.vy());
                        registry.fill(HIST("D0decayz"), particle0.vz());

                        //Actual particles information
                        if(fabs(mcMom.pdgCode())==421){
                           signMUON++;
                           for(auto Daughter : Daughters){
                              if(fabs(Daughter.pdgCode())==321){
                                 truthK_ID = Daughter.globalIndex();
                              }
                              if(fabs(Daughter.pdgCode())==13){
                                 vx_mu = particle0.vx();
                                 vy_mu = particle0.vy();
                                 vz_mu = particle0.vz();
                                 px_mu = particle0.px();
                                 py_mu = particle0.py();
                                 pz_mu = particle0.pz();
                                 t_mu = (-46-vz_mu)/pz_mu;
                                 mft_mu_x = px_mu*t_mu;
                                 mft_mu_y = py_mu*t_mu;
                                 s_mu = vz_mu/(vz_mu+46);
                                 //mft_det_z = track0.z();
                                 //mft_det_x = track0.x();
                                 //mft_det_y = track0.y();
                                 dca_mu_x = vx_mu+s_mu*(mft_mu_x-vx_mu);
                                 dca_mu_y = vy_mu+s_mu*(mft_mu_y-vy_mu);
                                 muonID = particle0.globalIndex();
                                 
                                 //cout << "mftX: " << mft_det_x << "  mftY: " << mft_det_y << " mftZ: " << mft_det_z << endl;
                                 a = dca_mu_x;
                                 b = mft_mu_x;
                                 c = dca_mu_y;
                                 d = mft_mu_y;
                                 //cout << "mft_MC_x: " << b << "  mft_MC_y: " << d << endl;
                                 registry.fill(HIST("muTrack_dcax"), a);
                                 registry.fill(HIST("muTrack_dcay"), c);
                                 registry.fill(HIST("muTrack_mftx"), b);
                                 registry.fill(HIST("muTrack_mfty"), d);
                              }
                           }
                        }
                        if(signMUON==0) continue;
                        for(auto& track1 : tracks){
                           if(!track1.has_mcParticle()) continue;
                           auto particle1 = track1.mcParticle();
                           if(particle1.globalIndex()==muonID) continue;
                           vx_can = particle1.vx();
                           vy_can = particle1.vy();
                           vz_can = particle1.vz();
                           px_can = particle1.px();
                           py_can = particle1.py();
                           pz_can = particle1.pz();
                           t_can = (-46-vz_can)/pz_can;
                           mft_can_x = px_can*t_can;
                           mft_can_y = py_can*t_can;
                           s_can = vz_can/(vz_can+46);
                           dca_can_x = vx_can+s_can*(mft_can_x-vx_can);
                           dca_can_y = vy_can+s_can*(mft_can_y-vy_can);
                           e = dca_can_x;
                           f = mft_can_x;
                           g = dca_can_y;
                           h = mft_can_y;

                           s = ((-e*f+a*f-g*h+c*h)*(b*f+d*h-pow(46,2)) + (e*b-a*b+d*g-c*d)*(pow(f,2)+pow(h,2)+pow(46,2))) / ((pow(b,2)+pow(d,2)+pow(46,2))*(pow(f,2)+pow(h,2)+pow(46,2)) - pow(b*f+d*h-pow(46,2),2));
                           t = ((b*f+d*h-pow(46,2))*s - e*f + a*f - g*h + c*h) / (pow(f,2)+pow(h,2)+pow(46,2));

                           mu_x = dca_mu_x+s*mft_mu_x;
                           mu_y = dca_mu_y+s*mft_mu_y;
                           mu_z = -46*s;
                           can_x = dca_can_x+t*mft_can_x;
                           can_y = dca_can_y+t*mft_can_y;
                           can_z = -46*t;

                           float distance_track = pow(mu_x-can_x,2)+pow(mu_y-can_y,2)+pow(mu_z-can_z,2);
                           if(distance_track<=closest_track){
                              closest_track = distance_track;
                              pca_x = (mu_x+can_x)/2;
                              pca_y = (mu_y+can_y)/2;
                              pca_z = (mu_z+can_z)/2;
                              estpdg = particle1.pdgCode();
                              estID = particle1.globalIndex();
                           }
                        }
                        registry.fill(HIST("pcax"), pca_x);
                        registry.fill(HIST("pcay"), pca_y);
                        registry.fill(HIST("pcaz"), pca_z);
                        registry.fill(HIST("PIDpurity"), fabs(estID-truthK_ID));
                        registry.fill(HIST("EstimatePDGcode"), estpdg);
                        cout << "Distance of each track: " << closest_track << endl;
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
