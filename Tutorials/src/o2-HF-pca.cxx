/// \author Koki Soeda
/// \since 25/11/2022

#include <iostream>
#include <cmath>
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
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

struct AccessMcData {
  Service<TDatabasePDG> pdg;

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  ConfigurableAxis PtAxis{"PtAxis", {1001, -0.05, 10.05}, "pt axis for histograms"};
  Configurable<float> zMax{"zMax", 5., "value for Zvtx cut"};

  HistogramRegistry registry{
    "registry",
    {   
      {"TracksPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"TracksToPartPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"TracksToPartPtEtaFakeMcColl", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"TracksToPartPtEtaPrim", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}}, 
      //{"TracksPDGcode", " ; PDG code", {HistType::kTH1F, {{12000, -6000, 6000}}}},
      {"TracksPDGcode", " ; PDG code", {HistType::kTH1F, {{9, -0.5, 8.5}}}},
      {"MomPdgCode2", " ; Mom PDG code", {HistType::kTH1F, {{12000, -6000, 6000}}}},
      {"MomPdgCode", "", {HistType::kTH1F, {{19, -0.5, 18.5}}}},
      {"TracksPt", "; #pt; ", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"TracksEta", "; #eta; ", {HistType::kTH1F, {{35, -4.5, -1.}}}},
      {"TracksVx", "; #vx; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"TracksVy", "; #vy; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"TracksVz", "; #vz; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},

      {"RecMass", "; #mass; ", {HistType::kTH1F, {{1000, 0., 10.0}}}},
      
      {"PiToPartPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"PiToPartPtEtaFakeMcColl", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"PiToPartPtEtaPrim", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}}, 
      {"PiPDGcode", " ; PDG code", {HistType::kTH1F, {{12000, -6000, 6000}}}},
      {"MomOfPiPdgCode", "", {HistType::kTH1F, {{19, -0.5, 18.5}}}},
      {"mPiPt", "; #pt; ", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"mPiVx", "; #vx; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPiVy", "; #vy; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPiVz", "; #vz; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimPiPt", "; #pt; ", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"mPrimPiVx", "; #vx; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimPiVy", "; #vy; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimPiVz", "; #vz; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},

      {"KToPartPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"KToPartPtEtaFakeMcColl", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"KToPartPtEtaPrim", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}}, 
      {"KPDGcode", " ; PDG code", {HistType::kTH1F, {{12000, -6000, 6000}}}},
      {"MomOfKPdgCode", " ; Mom of kaon PDG code", {HistType::kTH1F, {{12000, -6000, 6000}}}},
      {"mKPt", "; #pt; ", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"mKVx", "; #vx; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mKVy", "; #vy; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mKVz", "; #vz; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimKPt", "; #pt; ", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"mPrimKVx", "; #vx; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimKVy", "; #vy; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimKVz", "; #vz; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},

      {"MuToPartPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"MuToPartPtEtaFakeMcColl", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}},
      {"MuToPartPtEtaPrim", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}}, 
      {"MuPDGcode", " ; PDG code", {HistType::kTH1F, {{12000, -6000, 6000}}}},
      {"MomOfMuPdgCode", "", {HistType::kTH1F, {{12000, -6000, 6000}}}},
      {"mMuPt", "; #pt; ", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"mMuVx", "; #vx; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mMuVy", "; #vy; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mMuVz", "; #vz; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimMuPt", "; #pt; ", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"mPrimMuVx", "; #vx; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimMuVy", "; #vy; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},
      {"mPrimMuVz", "; #vz; ", {HistType::kTH1F, {{2000, -100., 100.0}}}},

    }
  };

  void init(InitContext&)
  {
    auto hpdgcode = registry.get<TH1>(HIST("TracksPDGcode"));
    auto* x1 = hpdgcode->GetXaxis();
    x1->SetBinLabel(1, "Electron"); // 11
    x1->SetBinLabel(2, "Muon"); // 13
    x1->SetBinLabel(3, "Pion±"); // 211
    x1->SetBinLabel(4, "K±"); // 321
    x1->SetBinLabel(5, "Proton"); // 2212
    x1->SetBinLabel(6, "sigma-"); // 3112
    x1->SetBinLabel(7, "sigma+"); // 3222
    x1->SetBinLabel(8, "Xi"); // 3312

    auto hmpdgcode = registry.get<TH1>(HIST("MomPdgCode"));
    auto* x2 = hmpdgcode->GetXaxis();
    x2->SetBinLabel(1, "Primary"); // 1 - 37
    x2->SetBinLabel(2, "Pion0"); // 111
    x2->SetBinLabel(3, "Pion±"); // 211
    x2->SetBinLabel(4, "pho"); // 113
    x2->SetBinLabel(5, "K_l"); // 130
    x2->SetBinLabel(6, "Eta"); // 221
    x2->SetBinLabel(7, "Omega"); // 223
    x2->SetBinLabel(8, "K_s"); // 310
    x2->SetBinLabel(9, "K*0(892)"); // 313
    x2->SetBinLabel(10, "K±"); // 321
    x2->SetBinLabel(11, "K*±(892)"); // 323
    x2->SetBinLabel(12, "Eta_prim"); // 331
    x2->SetBinLabel(13, "Phi"); // 333
    x2->SetBinLabel(14, "D±"); // 411
    x2->SetBinLabel(15, "D*±"); // 413
    x2->SetBinLabel(16, "D0"); // 421
    x2->SetBinLabel(17, "D_s±"); // 431
    x2->SetBinLabel(18, "Beauty"); // 500-599
    x2->SetBinLabel(19, "Baryon"); // 1000 -

    auto hmpdgcodeToPi = registry.get<TH1>(HIST("MomOfPiPdgCode"));
    auto* x3 = hmpdgcodeToPi->GetXaxis();
    x3->SetBinLabel(1, "Primary"); // 1 - 37
    x3->SetBinLabel(2, "pho"); // 113
    x3->SetBinLabel(3, "K_l"); // 130
    x3->SetBinLabel(4, "Eta"); // 221
    x3->SetBinLabel(5, "Omega"); // 223
    x3->SetBinLabel(6, "K_s"); // 310
    x3->SetBinLabel(7, "K*0(892)"); // 313
    x3->SetBinLabel(8, "K±"); // 321
    x3->SetBinLabel(9, "K*±(892)"); // 323
    x3->SetBinLabel(10, "Eta_prim"); // 331
    x3->SetBinLabel(11, "Phi"); // 333
    x3->SetBinLabel(12, "D±"); // 411
    x3->SetBinLabel(13, "D*±"); // 413
    x3->SetBinLabel(14, "D0"); // 421
    x3->SetBinLabel(15, "D_s±"); // 431
    x3->SetBinLabel(16, "B0"); // 511
    x3->SetBinLabel(17, "B+"); // 521
    x3->SetBinLabel(18, "B0_s"); // 531
    x3->SetBinLabel(19, "Baryon"); // 1000 -    
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::MFTTracks const& tracks)
  {
    if (!useEvSel || (useEvSel && collision.sel8())) {

      if ((collision.posZ() < zMax) && (collision.posZ() > -zMax)) {
        for (auto& track : tracks) {
          registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
        }
      }
    }
  }

  // group according to McCollisions
  void processGen(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const& particleMC, aod::McCollisions const&)
  {
    int i=0;
    // In the MFT the measurement of pT is not precise, so we access it by using the particle's pT instead
    if (collision.has_mcCollision()) {
      if ((collision.mcCollision().posZ() < zMax) && (collision.mcCollision().posZ() > -zMax)){
        if (!useEvSel || (useEvSel && collision.sel8())) {
          for (auto& track : tracks) {
            if (!track.has_mcParticle()) continue;
            auto particle = track.mcParticle(); // this mcParticle doesn't necessarily come from the right mcCollision
            registry.fill(HIST("TracksToPartPtEta"), particle.pt(), particle.eta());

            if (particle.mcCollisionId() == collision.mcCollision().globalIndex()){ // mcParticles come from the right mcCollision
              registry.fill(HIST("TracksPDGcode"), particle.pdgCode());
              registry.fill(HIST("TracksPt"), particle.pt());
              registry.fill(HIST("TracksEta"), particle.eta());
              registry.fill(HIST("TracksVx"), particle.vx());
              registry.fill(HIST("TracksVy"), particle.vy());
              registry.fill(HIST("TracksVz"), particle.vz());

              if(fabs(particle.pdgCode()) ==  11 )  registry.fill(HIST("TracksPDGcode"), 0.0, 1);
              if(fabs(particle.pdgCode()) ==  13 )  registry.fill(HIST("TracksPDGcode"), 1.0, 1);
              if(fabs(particle.pdgCode()) == 211 )  registry.fill(HIST("TracksPDGcode"), 2.0, 1);
              if(fabs(particle.pdgCode()) == 321 )  registry.fill(HIST("TracksPDGcode"), 3.0, 1);
              if(fabs(particle.pdgCode()) == 2212 ) registry.fill(HIST("TracksPDGcode"), 4.0, 1);
              if(fabs(particle.pdgCode()) == 3112 ) registry.fill(HIST("TracksPDGcode"), 5.0, 1);
              if(fabs(particle.pdgCode()) == 3222 ) registry.fill(HIST("TracksPDGcode"), 6.0, 1);
              if(fabs(particle.pdgCode()) == 3312 ) registry.fill(HIST("TracksPDGcode"), 7.0, 1);


              if(particle.has_mothers()){
              auto mcMom = particleMC.rawIteratorAt(particle.mothersIds()[0]);
                if(mcMom.daughtersIds().back() == particle.globalIndex()){
                  if(fabs(mcMom.pdgCode()) < 38 )   registry.fill(HIST("MomPdgCode"), 0.0, 1);
                  if(fabs(mcMom.pdgCode()) == 111 ) registry.fill(HIST("MomPdgCode"), 1.0, 1);
                  if(fabs(mcMom.pdgCode()) == 211 ) registry.fill(HIST("MomPdgCode"), 2.0, 1);
                  if(fabs(mcMom.pdgCode()) == 113 ) registry.fill(HIST("MomPdgCode"), 3.0, 1);
                  if(fabs(mcMom.pdgCode()) == 130 ) registry.fill(HIST("MomPdgCode"), 4.0, 1);
                  if(fabs(mcMom.pdgCode()) == 221 ) registry.fill(HIST("MomPdgCode"), 5.0, 1);
                  if(fabs(mcMom.pdgCode()) == 223 ) registry.fill(HIST("MomPdgCode"), 6.0, 1);
                  if(fabs(mcMom.pdgCode()) == 310 ) registry.fill(HIST("MomPdgCode"), 7.0, 1);
                  if(fabs(mcMom.pdgCode()) == 313 ) registry.fill(HIST("MomPdgCode"), 8.0, 1);
                  if(fabs(mcMom.pdgCode()) == 321 ) registry.fill(HIST("MomPdgCode"), 9.0, 1);
                  if(fabs(mcMom.pdgCode()) == 323 ) registry.fill(HIST("MomPdgCode"), 10.0, 1);
                  if(fabs(mcMom.pdgCode()) == 331 ) registry.fill(HIST("MomPdgCode"), 11.0, 1);
                  if(fabs(mcMom.pdgCode()) == 333 ) registry.fill(HIST("MomPdgCode"), 12.0, 1);
                  if(fabs(mcMom.pdgCode()) == 411 ) registry.fill(HIST("MomPdgCode"), 13.0, 1);
                  if(fabs(mcMom.pdgCode()) == 413 ) registry.fill(HIST("MomPdgCode"), 14.0, 1);
                  if(fabs(mcMom.pdgCode()) == 421 ) registry.fill(HIST("MomPdgCode"), 15.0, 1);
                  if(fabs(mcMom.pdgCode()) == 431 ) registry.fill(HIST("MomPdgCode"), 16.0, 1);
                  if(fabs(mcMom.pdgCode()) > 499 && fabs(mcMom.pdgCode()) < 600) registry.fill(HIST("MomPdgCode"), 17.0, 1);
                  if(fabs(mcMom.pdgCode()) > 999 ) registry.fill(HIST("MomPdgCode"), 18.0, 1);
                  registry.fill(HIST("MomPdgCode2"), mcMom.pdgCode());
                }

              }
            }

            if ((particle.mcCollisionId() != collision.mcCollision().globalIndex()) || (!particle.isPhysicalPrimary())) {
              if (particle.mcCollisionId() != collision.mcCollision().globalIndex()) {
                registry.fill(HIST("TracksToPartPtEtaFakeMcColl"), particle.pt(), particle.eta());
              }
              continue;
            }
            // mcParticles come from the rght collisions and primary
            registry.fill(HIST("TracksToPartPtEtaPrim"), particle.pt(), particle.eta());
          } 

          for(auto& [t0, t1] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))){
            if (!t0.has_mcParticle() || !t1.has_mcParticle()) continue;
            auto p0 = t0.mcParticle();
            auto p1 = t1.mcParticle();
            if(!p0.has_mothers() || !p1.has_mothers()) continue;
            auto Mom = particleMC.rawIteratorAt(p0.mothersIds()[0]);
            if (p0.mcCollisionId() == collision.mcCollision().globalIndex()) continue;
            if (p1.mcCollisionId() == collision.mcCollision().globalIndex()) continue;
            // if (p0.mothersIds().size() != 1 || p1.mothersIds().size() != 1) continue;
            // if (p0.mothersIds()[0] != p1.mothersIds()[0]) continue;
            //if (Mom.pdgCode() != 113) continue;
            auto mass2 = (p0.e()+p1.e())*(p0.e()+p1.e()) - (p0.px()+p1.px())*(p0.px()+p1.px()) - (p0.py()+p1.py())*(p0.py()+p1.py()) - (p0.pz()+p1.pz())*(p0.pz()+p1.pz());
            auto mass = sqrt(mass2);

            registry.fill(HIST("RecMass"), mass);
            registry.fill(HIST("MomPdgCode"), Mom.pdgCode());
          }
        }
      }
    } 
  }
  PROCESS_SWITCH(AccessMcData, processGen, "Process particle-level info", true);

  void processPi(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const& particleMC, aod::McCollisions const&)
  {
    // In the MFT the measurement of pT is not precise, so we access it by using the particle's pT instead
    if (collision.has_mcCollision()){
      if ((collision.mcCollision().posZ() < zMax) && (collision.mcCollision().posZ() > -zMax)) {
        if (!useEvSel || (useEvSel && collision.sel8())) {
          for (auto& track : tracks) {
            if (!track.has_mcParticle()) continue;
              
            auto particle = track.mcParticle(); // this mcParticle doesn't necessarily come from the right mcCollision
                
            if(fabs(particle.pdgCode()) != 211) continue;

            registry.fill(HIST("PiToPartPtEta"), particle.pt(), particle.eta());

            if (particle.mcCollisionId() == collision.mcCollision().globalIndex()) { // mcParticles come from the right mcCollision
              registry.fill(HIST("PiPDGcode"), particle.pdgCode());
              registry.fill(HIST("mPiPt"), particle.pt());
              registry.fill(HIST("mPiVx"), particle.vx());
              registry.fill(HIST("mPiVy"), particle.vy());
              registry.fill(HIST("mPiVz"), particle.vz());

              if(particle.has_mothers()){
              auto mcMom = particleMC.rawIteratorAt(particle.mothersIds()[0]);
                if(mcMom.daughtersIds().back() == particle.globalIndex()){
                  if(fabs(mcMom.pdgCode()) < 38 )   registry.fill(HIST("MomOfPiPdgCode"), 0.0, 1);
                  if(fabs(mcMom.pdgCode()) == 113 ) registry.fill(HIST("MomOfPiPdgCode"), 1.0, 1);
                  if(fabs(mcMom.pdgCode()) == 130 ) registry.fill(HIST("MomOfPiPdgCode"), 2.0, 1);
                  if(fabs(mcMom.pdgCode()) == 221 ) registry.fill(HIST("MomOfPiPdgCode"), 3.0, 1);
                  if(fabs(mcMom.pdgCode()) == 223 ) registry.fill(HIST("MomOfPiPdgCode"), 4.0, 1);
                  if(fabs(mcMom.pdgCode()) == 310 ) registry.fill(HIST("MomOfPiPdgCode"), 5.0, 1);
                  if(fabs(mcMom.pdgCode()) == 313 ) registry.fill(HIST("MomOfPiPdgCode"), 6.0, 1);
                  if(fabs(mcMom.pdgCode()) == 321 ) registry.fill(HIST("MomOfPiPdgCode"), 7.0, 1);
                  if(fabs(mcMom.pdgCode()) == 323 ) registry.fill(HIST("MomOfPiPdgCode"), 8.0, 1);
                  if(fabs(mcMom.pdgCode()) == 331 ) registry.fill(HIST("MomOfPiPdgCode"), 9.0, 1);
                  if(fabs(mcMom.pdgCode()) == 333 ) registry.fill(HIST("MomOfPiPdgCode"), 10.0, 1);
                  if(fabs(mcMom.pdgCode()) == 411 ) registry.fill(HIST("MomOfPiPdgCode"), 11.0, 1);
                  if(fabs(mcMom.pdgCode()) == 413 ) registry.fill(HIST("MomOfPiPdgCode"), 12.0, 1);
                  if(fabs(mcMom.pdgCode()) == 421 ) registry.fill(HIST("MomOfPiPdgCode"), 13.0, 1);
                  if(fabs(mcMom.pdgCode()) == 431 ) registry.fill(HIST("MomOfPiPdgCode"), 14.0, 1);
                  if(fabs(mcMom.pdgCode()) == 511 ) registry.fill(HIST("MomOfPiPdgCode"), 15.0, 1);
                  if(fabs(mcMom.pdgCode()) == 521 ) registry.fill(HIST("MomOfPiPdgCode"), 16.0, 1);
                  if(fabs(mcMom.pdgCode()) == 531 ) registry.fill(HIST("MomOfPiPdgCode"), 17.0, 1);
                  if(fabs(mcMom.pdgCode()) > 999 ) registry.fill(HIST("MomOfPiPdgCode"), 18.0, 1);
                  //registry.fill(HIST("MomOfPiPdgCode"), mcMom.pdgCode());
                }
              }
            }
            if ((particle.mcCollisionId() != collision.mcCollision().globalIndex()) || (!particle.isPhysicalPrimary())) {
              if (particle.mcCollisionId() != collision.mcCollision().globalIndex()) {
                registry.fill(HIST("PiToPartPtEtaFakeMcColl"), particle.pt(), particle.eta());
              }
              continue;
            }
            // mcParticles come from the rght collisions and primary
            registry.fill(HIST("PiToPartPtEtaPrim"), particle.pt(), particle.eta());
            registry.fill(HIST("mPrimPiPt"), particle.pt());
            registry.fill(HIST("mPrimPiVx"), particle.vx());
            registry.fill(HIST("mPrimPiVy"), particle.vy());
            registry.fill(HIST("mPrimPiVz"), particle.vz());
          }
        }
      }
    }  
  }
  PROCESS_SWITCH(AccessMcData, processPi, "Process particle-level info for Pi", true);

  void processK(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const& particleMC, aod::McCollisions const&)
  {
    // In the MFT the measurement of pT is not precise, so we access it by using the particle's pT instead
    if (collision.has_mcCollision()){
      if ((collision.mcCollision().posZ() < zMax) && (collision.mcCollision().posZ() > -zMax)) {
        if (!useEvSel || (useEvSel && collision.sel8())) {
          for (auto& track : tracks) {
            if (!track.has_mcParticle()) continue;
              
            auto particle = track.mcParticle(); // this mcParticle doesn't necessarily come from the right mcCollision
                
            if(fabs(particle.pdgCode()) != 321) continue;

            registry.fill(HIST("KToPartPtEta"), particle.pt(), particle.eta());

            if (particle.mcCollisionId() == collision.mcCollision().globalIndex()) { // mcParticles come from the right mcCollision
              registry.fill(HIST("KPDGcode"), particle.pdgCode());
              registry.fill(HIST("mKPt"), particle.pt());
              registry.fill(HIST("mKVx"), particle.vx());
              registry.fill(HIST("mKVy"), particle.vy());
              registry.fill(HIST("mKVz"), particle.vz());

              if(particle.has_mothers()){
              auto mcMom = particleMC.rawIteratorAt(particle.mothersIds()[0]);
                if(mcMom.daughtersIds().back() == particle.globalIndex()){
                  registry.fill(HIST("MomOfKPdgCode"), mcMom.pdgCode());
                }
              }
            }
            if ((particle.mcCollisionId() != collision.mcCollision().globalIndex()) || (!particle.isPhysicalPrimary())) {
              if (particle.mcCollisionId() != collision.mcCollision().globalIndex()) {
                registry.fill(HIST("KToPartPtEtaFakeMcColl"), particle.pt(), particle.eta());
              }
              continue;
            }
            // mcParticles come from the rght collisions and primary
            registry.fill(HIST("KToPartPtEtaPrim"), particle.pt(), particle.eta());
            registry.fill(HIST("mPrimKPt"), particle.pt());
            registry.fill(HIST("mPrimKVx"), particle.vx());
            registry.fill(HIST("mPrimKVy"), particle.vy());
            registry.fill(HIST("mPrimKVz"), particle.vz());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(AccessMcData, processK, "Process particle-level info for K", true);

  void processMu(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const& particleMC, aod::McCollisions const&)
  {
    // In the MFT the measurement of pT is not precise, so we access it by using the particle's pT instead
    if (collision.has_mcCollision()){
      if ((collision.mcCollision().posZ() < zMax) && (collision.mcCollision().posZ() > -zMax)) {
        if (!useEvSel || (useEvSel && collision.sel8())) {
          for (auto& track : tracks) {
            if (!track.has_mcParticle()) continue;
              
            auto particle = track.mcParticle(); // this mcParticle doesn't necessarily come from the right mcCollision
                
            if(fabs(particle.pdgCode()) != 13) continue;

            registry.fill(HIST("MuToPartPtEta"), particle.pt(), particle.eta());

            if (particle.mcCollisionId() == collision.mcCollision().globalIndex()) { // mcParticles come from the right mcCollision
              registry.fill(HIST("MuPDGcode"), particle.pdgCode());
              registry.fill(HIST("mMuPt"), particle.pt());
              registry.fill(HIST("mMuVx"), particle.vx());
              registry.fill(HIST("mMuVy"), particle.vy());
              registry.fill(HIST("mMuVz"), particle.vz());

              if(particle.has_mothers()){
              auto mcMom = particleMC.rawIteratorAt(particle.mothersIds()[0]);
                if(mcMom.daughtersIds().back() == particle.globalIndex()){
                  // if(fabs(mcMom.pdgCode()) < 38 )   registry.fill(HIST("MomOfMuPdgCode"), 0.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 113 ) registry.fill(HIST("MomOfMuPdgCode"), 1.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 130 ) registry.fill(HIST("MomOfMuPdgCode"), 2.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 221 ) registry.fill(HIST("MomOfMuPdgCode"), 3.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 223 ) registry.fill(HIST("MomOfMuPdgCode"), 4.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 310 ) registry.fill(HIST("MomOfMuPdgCode"), 5.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 313 ) registry.fill(HIST("MomOfMuPdgCode"), 6.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 321 ) registry.fill(HIST("MomOfMuPdgCode"), 7.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 323 ) registry.fill(HIST("MomOfMuPdgCode"), 8.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 331 ) registry.fill(HIST("MomOfMuPdgCode"), 9.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 333 ) registry.fill(HIST("MomOfMuPdgCode"), 10.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 411 ) registry.fill(HIST("MomOfMuPdgCode"), 11.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 413 ) registry.fill(HIST("MomOfMuPdgCode"), 12.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 421 ) registry.fill(HIST("MomOfMuPdgCode"), 13.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 431 ) registry.fill(HIST("MomOfMuPdgCode"), 14.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 511 ) registry.fill(HIST("MomOfMuPdgCode"), 15.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 521 ) registry.fill(HIST("MomOfMuPdgCode"), 16.0, 1);
                  // if(fabs(mcMom.pdgCode()) == 531 ) registry.fill(HIST("MomOfMuPdgCode"), 17.0, 1);
                  // if(fabs(mcMom.pdgCode()) > 999 ) registry.fill(HIST("MomOfMuPdgCode"), 18.0, 1);
                  registry.fill(HIST("MomOfMuPdgCode"), mcMom.pdgCode());
                }
              }
            }
            if ((particle.mcCollisionId() != collision.mcCollision().globalIndex()) || (!particle.isPhysicalPrimary())) {
              if (particle.mcCollisionId() != collision.mcCollision().globalIndex()) {
                registry.fill(HIST("MuToPartPtEtaFakeMcColl"), particle.pt(), particle.eta());
              }
              continue;
            }
            // mcParticles come from the rght collisions and primary
            registry.fill(HIST("MuToPartPtEtaPrim"), particle.pt(), particle.eta());
            registry.fill(HIST("mPrimMuPt"), particle.pt());
            registry.fill(HIST("mPrimMuVx"), particle.vx());
            registry.fill(HIST("mPrimMuVy"), particle.vy());
            registry.fill(HIST("mPrimMuVz"), particle.vz());
          }
        }
      }
    }  
  }
  PROCESS_SWITCH(AccessMcData, processMu, "Process particle-level info for Mu", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexdistribution>(cfgc),
                      adaptAnalysisTask<AccessMcData>(cfgc),
  };
}