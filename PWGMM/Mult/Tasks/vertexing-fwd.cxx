/// \author Koki Soeda
/// \since 16/12/2022

#include <cmath>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MathUtils/Utils.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using CollisionsLabeled = soa::Join<o2::aod::Collisions, aod::McCollisionLabels>;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

AxisSpec ZAxis = {301, -30.1, 30.1};

struct vertexingfwd {

  Configurable<float> maxDCAXY{"maxDCAXY", 6.0, "max allowed transverse DCA"}; // To be used when associating ambitrack to collision using best DCA
  std::vector<uint64_t> ambTrackIds;

  HistogramRegistry registry{
    "registry",
    {{"ParticleZR", "; #it{z}_{vtx} (cm); #it{R} (cm) ; count", {HistType::kTH2F, {ZAxis, {1001, -0.5, 100.5}}}},         //
     {"Primary/ParticleZR", "; #it{z}_{vtx} (cm); #it{R} (cm) ; count", {HistType::kTH2F, {ZAxis, {1001, -0.5, 100.5}}}}, //

     {"Truth/DCAxTruth", "; DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Truth/DCAyTruth", "; DCA_{y}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

     {"Primary/DCAxPrimNAmb", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAyPrimNAmb", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAxPtPrimNAmb", "; DCA_{x} (cm); #it{p}_{T}^{true} (GeV/#it{c}) ; counts", {HistType::kTH2F, {{2500, -5, 5}, {21, -0.05, 10.05}}}},

     {"Primary/DCAxPrimAmb1", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAyPrimAmb1", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

     {"Primary/DCAxPrimAmb2", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAyPrimAmb2", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

     // ambiguous tracks and DCA with all compatible collisions
     {"Ambiguous/TracksDCAXY", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"Ambiguous/TracksDCAX", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"Ambiguous/TracksDCAY", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     // all tracks (before amb reassociation)
     {"DCAXYAll", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"DCAXAll", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"DCAYAll", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     // ambiguous tracks before reassociation
     {"DCAXYAmb", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"DCAXAmb", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"DCAYAmb", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     // non ambiguous tracks
     {"DCAXYNAmb", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"DCAXNAmb", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"DCAYNAmb", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     {"AmbiguousTrackStatus", "; ; N_{trk}", {HistType::kTH1F, {{8, 0.5, 8.5}}}},

     // DCAxy, x and y distributions for reassociated ambiguous tracks
     // when it is a false reassociation and when it is true
     {"Ambiguous/TracksDCAXYBest", "; DCA_{xy}^{best} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"Ambiguous/TracksDCAXBest", "; DCA_{x}^{best} (cm); counts", {HistType::kTH1F, {{500, -10, 10}}}},
     {"Ambiguous/TracksDCAYBest", "; DCA_{y}^{best} (cm); counts", {HistType::kTH1F, {{500, -10, 10}}}},

     {"Ambiguous/TracksDCAXYBestTrue", "; DCA_{xy}^{best, true} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"Ambiguous/TracksDCAXYBestFalse", "; DCA_{xy}^{best, false} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},

     {"Ambiguous/TracksDCAXBestTrue", "; DCA_{x}^{best, true} (cm); counts", {HistType::kTH1F, {{500, -10, 10}}}},
     {"Ambiguous/TracksDCAYBestTrue", "; DCA_{y}^{best, true} (cm); counts", {HistType::kTH1F, {{500, -10, 10}}}},
     {"Ambiguous/TracksDCAXBestFalse", "; DCA_{x}^{best, false} (cm); counts", {HistType::kTH1F, {{500, -10, 10}}}},
     {"Ambiguous/TracksDCAYBestFalse", "; DCA_{y}^{best, false} (cm); counts", {HistType::kTH1F, {{500, -10, 10}}}}

    }};

  void init(InitContext&)
  {
    auto hstatus = registry.get<TH1>(HIST("AmbiguousTrackStatus"));
    auto* x2 = hstatus->GetXaxis();
    x2->SetBinLabel(1, "MFT tracks");
    x2->SetBinLabel(2, "MFT ambiguous tracks");
    x2->SetBinLabel(3, "Reassigned tracks");
    x2->SetBinLabel(4, "Extra tracks");
    x2->SetBinLabel(5, "orig=true (re)");
    x2->SetBinLabel(6, "best=true (re)");
    x2->SetBinLabel(7, "not reassigned");
    x2->SetBinLabel(8, "not reassigned and true");
  }

  void processDCAamb(MFTTracksLabeled const&,
                     CollisionsLabeled const&, ExtBCs const& bcs,
                     aod::AmbiguousMFTTracks const& atracks,
                     aod::McParticles const&,
                     aod::McCollisions const&)
  {

    if (bcs.size() == 0) {
      return;
    }
    if (atracks.size() == 0) {
      return;
    }
    // initCCDB(bcs.begin()); if Bz is needed
    ambTrackIds.clear();

    float dcaXY;
    float bestDCA, bestDCAX, bestDCAY;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto& atrack : atracks) {

      dcaXY = 999; // DCAxy
      bestDCA = 999;

      auto track = atrack.mfttrack_as<MFTTracksLabeled>();
      if (track.has_collision()) {
        ambTrackIds.push_back(atrack.mfttrackId());
      }
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();

      auto bestCol = track.has_collision() ? track.collisionId() : -1;
      int bestMCCol = -1;
      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      auto compatibleBCs = atrack.bc_as<ExtBCs>();

      for (auto& bc : compatibleBCs) {
        if (!bc.has_collisions()) {
          continue;
        }
        auto collisions = bc.collisions_as<CollisionsLabeled>(); // compatible collisions
        for (auto const& collision : collisions) {

          // trackPar.propagateToZhelix(collision.posZ(), Bz); // track parameters propagation to the position of the z vertex
          trackPar.propagateToZlinear(collision.posZ());

          const auto dcaX(trackPar.getX() - collision.posX());
          const auto dcaY(trackPar.getY() - collision.posY());
          dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

          registry.fill(HIST("Ambiguous/TracksDCAXY"), dcaXY);
          registry.fill(HIST("Ambiguous/TracksDCAX"), dcaX);
          registry.fill(HIST("Ambiguous/TracksDCAY"), dcaY);
          // registry.fill(HIST("NumberOfContributors"), collision.numContrib());

          if ((dcaXY < bestDCA)) {
            bestCol = collision.globalIndex();
            bestMCCol = collision.mcCollisionId();
            bestDCA = dcaXY;
            bestDCAX = dcaX;
            bestDCAY = dcaY;
            bestTrackPar = trackPar;
          }
        }
      }

      // other option for the truth : collision.mcCollision().posZ();

      registry.fill(HIST("Ambiguous/TracksDCAXYBest"), bestDCA);
      registry.fill(HIST("Ambiguous/TracksDCAXBest"), bestDCAX);
      registry.fill(HIST("Ambiguous/TracksDCAYBest"), bestDCAY);

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("Primary/DCAxPrimAmb2"), bestDCAX);
        registry.fill(HIST("Primary/DCAyPrimAmb2"), bestDCAY);
      }

      auto mcCollID = particle.mcCollisionId();
      registry.fill(HIST("AmbiguousTrackStatus"), 2);
      if (!track.has_collision()) {
        // this is a track extra
        registry.fill(HIST("AmbiguousTrackStatus"), 4);
        if (bestMCCol == mcCollID) // correctly assigned to bestCol
        {
          registry.fill(HIST("AmbiguousTrackStatus"), 6);
          registry.fill(HIST("Ambiguous/TracksDCAXYBestTrue"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestTrue"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestTrue"), bestDCAY);
        } else {
          registry.fill(HIST("Ambiguous/TracksDCAXYBestFalse"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestFalse"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestFalse"), bestDCAY);
        }

      } else if (track.collisionId() != bestCol) {
        // this track has been reassigned
        auto collOrig = track.collision_as<CollisionsLabeled>();
        registry.fill(HIST("AmbiguousTrackStatus"), 3);
        if (bestMCCol == mcCollID) // correctly reassigned
        {
          registry.fill(HIST("AmbiguousTrackStatus"), 6);
          registry.fill(HIST("Ambiguous/TracksDCAXYBestTrue"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestTrue"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestTrue"), bestDCAY);
        } else { // uncorrectly reassigned
          registry.fill(HIST("Ambiguous/TracksDCAXYBestFalse"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestFalse"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestFalse"), bestDCAY);
        }

        if (collOrig.mcCollisionId() == mcCollID) { // initially correctly assigned
          registry.fill(HIST("AmbiguousTrackStatus"), 5);
        }
      } else // the track has a collision and track.collisionId() == bestCol
      {
        if (track.collisionId() != bestCol) {
          printf("------------------- PROBLEM HERE track.collisionId() %d, bestCollid %d\n", track.collisionId(), bestCol);
        }

        registry.fill(HIST("AmbiguousTrackStatus"), 7);
        if (bestMCCol == mcCollID) // correctly assigned
        {
          registry.fill(HIST("AmbiguousTrackStatus"), 8);
        }
      }
    }
  }
  PROCESS_SWITCH(vertexingfwd, processDCAamb, "get the DCAxy of MFT ambiguous tracks", true);

  void processDCA(MFTTracksLabeled const& tracks,
                  aod::Collisions const&,
                  aod::McParticles const&,
                  aod::McCollisions const&)
  {
    if (tracks.size() == 0) {
      return;
    }
    registry.fill(HIST("AmbiguousTrackStatus"), 1, tracks.size());
    float dcaXY;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto& track : tracks) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      if (!track.has_collision()) {
        continue;
      }
      auto collision = track.collision();

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      trackPar.propagateToZlinear(collision.posZ());

      const auto dcaX(trackPar.getX() - collision.posX());
      const auto dcaY(trackPar.getY() - collision.posY());
      dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

      float Rv = std::sqrt(pow(particle.vx(), 2) + pow(particle.vy(), 2));
      registry.fill(HIST("ParticleZR"), particle.vz(), Rv);

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("Primary/ParticleZR"), particle.vz(), Rv);
      }

      const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
      const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
      // auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth); // this is DCA_xy truth

      registry.fill(HIST("Truth/DCAxTruth"), dcaXtruth);
      registry.fill(HIST("Truth/DCAyTruth"), dcaYtruth);

      registry.fill(HIST("DCAXYAll"), dcaXY);
      registry.fill(HIST("DCAXAll"), dcaX);
      registry.fill(HIST("DCAYAll"), dcaY);

      if (find(ambTrackIds.begin(), ambTrackIds.end(), track.globalIndex()) != ambTrackIds.end()) {
        registry.fill(HIST("DCAXYAmb"), dcaXY);
        registry.fill(HIST("DCAXAmb"), dcaX);
        registry.fill(HIST("DCAYAmb"), dcaY);
        if (particle.isPhysicalPrimary()) {
          registry.fill(HIST("Primary/DCAxPrimAmb1"), dcaX);
          registry.fill(HIST("Primary/DCAyPrimAmb1"), dcaY);
        }

        continue; // this track has already been reassigned to bestcollision, don't double count
      }
      registry.fill(HIST("DCAXYNAmb"), dcaXY);
      registry.fill(HIST("DCAXNAmb"), dcaX);
      registry.fill(HIST("DCAYNAmb"), dcaY);

      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      // only tracks coming from primary particles here
      // non ambiguous tracks
      registry.fill(HIST("Primary/DCAxPrimNAmb"), dcaX);
      registry.fill(HIST("Primary/DCAyPrimNAmb"), dcaY);

      registry.fill(HIST("Primary/DCAxPtPrimNAmb"), dcaX, particle.pt());
    }
  }
  PROCESS_SWITCH(vertexingfwd, processDCA, "get the DCAxy of MFT tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexingfwd>(cfgc)};
}



/*
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// \file   vertexing-fwd.cxx
// \author Robin Caron <robin.caron@cern.ch>
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code loops over every ambiguous MFT tracks and associates
// them to a collision that has the smallest DCAxy

#include <cmath>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MathUtils/Utils.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using FullBCs = soa::Join<aod::BCs, aod::MatchedBCCollisionsSparse>;
using FullCollision = soa::Join<aod::Collisions, aod::McCollisionLabels>;

namespace o2::aod
{
DECLARE_SOA_TABLE(AmbiguousTracksMFT, "AOD", "AMBIGUOUSTRMFT", //! Table for MFT tracks which are not uniquely associated with a collision
                  o2::soa::Index<>, o2::aod::ambiguous::MFTTrackId, o2::aod::ambiguous::BCIdSlice, o2::soa::Marker<2>);
} // namespace o2::aod

struct vertexingfwd {

  /// Could be TEMPORARY: store the vertex, collision, dca, and ambiguous tracks information
  /// into different std::vector to easily handle them later outside the loops
  std::vector<int> vecCollForAmb;        // vector for collisions associated to an ambiguous track
  std::vector<double> vecDCACollForAmb;  // vector for dca collision associated to an ambiguous track
  std::vector<double> vecAmbTrack;       // vector for ambiguous track quantities like chi2 and z
  std::vector<double> vecZposCollForAmb; // vector for z vertex of collisions associated to an ambiguous track

  Configurable<float> maxDCAXY{"maxDCAXY", 6.0, "max allowed transverse DCA"}; // To be used when associating ambitrack to collision using best DCA

  HistogramRegistry registry{
    "registry",
    {{"TracksDCAXY", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"TracksDCAX", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{100, -10, 10}}}},
     {"TracksDCAY", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{100, -10, 10}}}},
     {"AmbiguousTracksStatus", "; Status; counts", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
     {"NbCollComp", "; NbCollComp", {HistType::kTH1F, {{20, -0.5, 19.5}}}},
     {"NumberOfContributors", "; N_{tr} for vertexing; counts", {HistType::kTH1F, {{100, 0, 100}}}},
     {"CollisionsMatchIndicesMC", "; Rec. minDCA ambitrack coll.ID; Gen. ambitrack coll.ID", {HistType::kTH2F, {{401, -0.5, 1000.5}, {401, -0.5, 1000.5}}}},
     {"TracksDCAXYBest", "; DCA_{xy}^{best} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"TracksDCAXYBestFalse", "; DCA_{xy}^{best, false} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"EfficiencyZvtx", "; z vertex; Efficiency", {HistType::kTProfile, {{100, -30, 30}}}},
     {"DeltaZvtx", "; #delta z (cm); counts", {HistType::kTH1F, {{400, -20, 20}}}},
     {"DeltaZvtxBest", "; #delta z = z_{best} - z_{true} (cm); counts", {HistType::kTH1F, {{400, -20, 20}}}},
     {"CorrectMatch", "; Matching value; counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}}}};

  void init(InitContext&)
  {
    auto hstat = registry.get<TH1>(HIST("CorrectMatch"));
    auto* x1 = hstat->GetXaxis();
    x1->SetBinLabel(1, "Incorrect match");
    x1->SetBinLabel(2, "Correct match");
    x1->SetBinLabel(3, "N_{ambitrack} associable");
    x1->SetBinLabel(4, "N_{ambitrack} w N_{coll} > 0");
    x1->SetBinLabel(5, "N_{ambitrack} total");

    auto hstatus = registry.get<TH1>(HIST("AmbiguousTracksStatus"));
    auto* x2 = hstatus->GetXaxis();
    x2->SetBinLabel(1, "MFT tracks ");
    x2->SetBinLabel(2, "MFT ambiguous tracks ");
    x2->SetBinLabel(3, "All ambiguous primary");
    x2->SetBinLabel(4, "Orphan + primary");
    x2->SetBinLabel(5, "All ambiguous secondary");
    x2->SetBinLabel(6, "Orphan + secondary");
  }

  int getIndexBestCollision(std::vector<double> vecOfDCA, int method = 0)
  {
    int indice = 0;
    if (vecOfDCA.size() == 0) {
      return -1;
    }
    if (method == 0) {
      indice = std::distance(vecOfDCA.begin(), std::min_element(vecOfDCA.begin(), vecOfDCA.end()));
    }
    return indice;
  }

  template <typename T>
  void doProcess(T const& ambitracks, aod::BCs const& bcs, MFTTracksLabeled const& tracks, FullCollision const& collisions, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    int ntracks = tracks.size();
    int nambitracks = ambitracks.size();

    registry.fill(HIST("AmbiguousTracksStatus"), 0.0, ntracks);
    registry.fill(HIST("AmbiguousTracksStatus"), 1.0, nambitracks);

    for (auto& ambitrack : ambitracks) {
      vecCollForAmb.clear();
      vecDCACollForAmb.clear();
      vecAmbTrack.clear();
      vecZposCollForAmb.clear();

      double value = 0.0;    // matching value for collision association to an ambiguous track
      double zVtxMCAmbi = 0; // z vertex associated to the mc collision
      int mcCollAmbiID = -1; // mc value for the collision containing the ambiguous track

      //auto track = ambitrack.mfttrack_as<MFTTracksLabeled>(); // Obtain the MFT ambiguous track with the MC labels
      auto track = ambitrack.template mfttrack_as<MFTTracksLabeled>(); //This is the synthax when calling a templated funcction on a template

      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      mcCollAmbiID = particle.mcCollisionId();
      zVtxMCAmbi = particle.mcCollision().posZ();

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("AmbiguousTracksStatus"), 2.0);
      } else {
        registry.fill(HIST("AmbiguousTracksStatus"), 4.0);
      }

      // Fill the std::vector for ambiguous tracks with the quantities needed
      vecAmbTrack.push_back(track.x());
      vecAmbTrack.push_back(track.y());
      vecAmbTrack.push_back(track.phi());
      vecAmbTrack.push_back(track.tgl());
      vecAmbTrack.push_back(track.signed1Pt());
      vecAmbTrack.push_back(track.z());
      vecAmbTrack.push_back(track.chi2());

      auto bcambis = ambitrack.bc();

      int collCounter = 0;
      for (auto& collision : collisions) {
        uint64_t mostProbableBC = collision.bc().globalBC();
        //uint64_t meanBC = mostProbableBC - std::lround(collision.collisionTime() / (o2::constants::lhc::LHCBunchSpacingNS / 1000));
        //int deltaBC = std::ceil(collision.collisionTimeRes() / (o2::constants::lhc::LHCBunchSpacingNS / 1000) * 4);

        for (auto& bcambi : bcambis) {

          if (bcambi.globalBC() != mostProbableBC) {
            continue;
          }
          //here the bc of the ambitrack is the bc of the collision we are looking at
          collCounter++;

          // We compute the DCAxy of this track wrt the primary vertex of the current collision
          SMatrix5 tpars(vecAmbTrack[0], vecAmbTrack[1], vecAmbTrack[2], vecAmbTrack[3], vecAmbTrack[4]);
          //            std::vector<double> v1{extAmbiTrack.cXX(), extAmbiTrack.cXY(), extAmbiTrack.cYY(), extAmbiTrack.cPhiX(), extAmbiTrack.cPhiY(),
          //                                   extAmbiTrack.cPhiPhi(), extAmbiTrack.cTglX(), extAmbiTrack.cTglY(), extAmbiTrack.cTglPhi(), extAmbiTrack.cTglTgl(),
          //                                   extAmbiTrack.c1PtX(), extAmbiTrack.c1PtY(), extAmbiTrack.c1PtPhi(), extAmbiTrack.c1PtTgl(), extAmbiTrack.c1Pt21Pt2()};

          std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
          SMatrix55 tcovs(v1.begin(), v1.end());
          o2::track::TrackParCovFwd pars1{vecAmbTrack[5], tpars, tcovs, vecAmbTrack[6]};

          // o2::track::TrackParCovFwd pars1{extAmbiTrack.z(), tpars, tcovs, chi2};
          pars1.propagateToZlinear(collision.posZ()); // track parameters propagation to the position of the z vertex

          const auto dcaX(pars1.getX() - collision.posX());
          const auto dcaY(pars1.getY() - collision.posY());
          auto dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

          registry.fill(HIST("TracksDCAXY"), dcaXY);
          registry.fill(HIST("TracksDCAX"), dcaX);
          registry.fill(HIST("TracksDCAY"), dcaY);
          registry.fill(HIST("NumberOfContributors"), collision.numContrib());

          if (dcaXY > maxDCAXY) {
            continue;
          }

          vecDCACollForAmb.push_back(dcaXY);

          if (!collision.has_mcCollision()) {
            continue;
          }

          int mcCollindex = collision.mcCollision().globalIndex();
          vecCollForAmb.push_back(mcCollindex);

          vecZposCollForAmb.push_back(collision.mcCollision().posZ());

          registry.fill(HIST("DeltaZvtx"), collision.mcCollision().posZ() - zVtxMCAmbi);
          break; //once we found a bc that corresponds to the bc of the collision we are working on, go out of the bcambi loop
          //We then look at the next collision
        }
      }

      registry.fill(HIST("NbCollComp"), collCounter);
      registry.fill(HIST("CorrectMatch"), 4.0); // counting for ambiguous track with N collisions >=0

      if (collCounter == 0) {                //these are orphan tracks
        if (!particle.isPhysicalPrimary()) { //orphan and secondary
          registry.fill(HIST("AmbiguousTracksStatus"), 5.0);
        }
        if (particle.isPhysicalPrimary()) { //orphan and primary
          registry.fill(HIST("AmbiguousTracksStatus"), 3.0);
        }
        continue;
      }

      registry.fill(HIST("CorrectMatch"), 3.0); // counting for ambiguous track with N collisions >0

      int indexMinDCA = getIndexBestCollision(vecDCACollForAmb, 0); // obtain min value in the stored vector of DCAs
      int indexMCcoll = -1;
      if (indexMinDCA == -1) {
        continue; //if no DCAxy was < maxDCAXY
      }
      indexMCcoll = vecCollForAmb[indexMinDCA];

      registry.fill(HIST("CollisionsMatchIndicesMC"), mcCollAmbiID, indexMCcoll);

      if (collCounter == 1) { //This shouldn't happen after Ruben's fix
        printf("strange ambiguous track of mfttrackId %d\n", ambitrack.mfttrackId());
        if (mcCollAmbiID == indexMCcoll) {
          printf("and this is a correct match for the ambiguous track of mfttrackid %d\n", ambitrack.mfttrackId());
        }
      }

      registry.fill(HIST("DeltaZvtxBest"), vecZposCollForAmb[indexMinDCA] - zVtxMCAmbi);

      if (mcCollAmbiID == indexMCcoll) { //correct match
        value = 1.0;
        // LOGF(info, " --> Ambitrack correctly associated to collision, dca= %f", vecDCACollForAmb[indexMinDCA]);
      }
      registry.fill(HIST("TracksDCAXYBest"), vecDCACollForAmb[indexMinDCA]);
      registry.fill(HIST("CorrectMatch"), value);
      registry.fill(HIST("EfficiencyZvtx"), zVtxMCAmbi, value);
      registry.fill(HIST("CorrectMatch"), 2.0); // Counting for amibuous track with N collisions > 0

      if (value == 0.0) {                                                           //incorrect match
        registry.fill(HIST("TracksDCAXYBestFalse"), vecDCACollForAmb[indexMinDCA]); // Incorrect association with min DCA
      }

    } // ambitracks loop
  }   //end of doProcess

  void processNew(aod::AmbiguousMFTTracks const& ambitracks, aod::BCs const& bcs, MFTTracksLabeled const& tracks, FullCollision const& collisions, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    doProcess(ambitracks, bcs, tracks, collisions, mcParticles, mcCollisions);
  }
  PROCESS_SWITCH(vertexingfwd, processNew, "Process ambiguous track DCA", true);

  void processOld(aod::AmbiguousTracksMFT const& ambitracks, aod::BCs const& bcs, MFTTracksLabeled const& tracks, FullCollision const& collisions, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    doProcess(ambitracks, bcs, tracks, collisions, mcParticles, mcCollisions);
  }
  PROCESS_SWITCH(vertexingfwd, processOld, "Process ambiguous track DCA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexingfwd>(cfgc)};
}
*/