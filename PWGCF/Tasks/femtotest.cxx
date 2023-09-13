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
//
/// \brief Use the HistogramRegistry to create histograms with distributions for track selection.
/// \author Sofia Tomassini
/// \since 31 May 2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <TParameter.h>
#include <TH1F.h>

#include "Framework/ASoA.h"
#include "MathUtils/Utils.h"
#include "Framework/DataTypes.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/StaticFor.h"
#include "PWGCF/DataModel/singletrackselector.h"

#include <vector>
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "PWGCF/DataModel/femtotest.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

double particle_mass(int PDGcode)
{
  // if(PDGcode == 2212) return TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  if (PDGcode == 1000010020)
    return 1.87561294257;
  else
    return TDatabasePDG::Instance()->GetParticle(PDGcode)->Mass();
}

struct FemtoCorrelations {
  // using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> min_P{"min_P", 0.0, ""};
  Configurable<float> max_P{"max_P", 100.0, ""};
  Configurable<float> eta{"eta", 100.0, ""};
  Configurable<int> tpcNClsFound{"tpcNClsFound", 0, ""};
  Configurable<float> tpcChi2NCl{"tpcChi2NCl", 100.0, ""};
  Configurable<float> dcaXY{"dcaXY", 10.0, ""};
  Configurable<float> dcaZ{"dcaZ", 10.0, ""};
  Configurable<int> itsNCls{"itsNCls", -1, ""};
  Configurable<int> VertexZ{"VertexZ", 10.0, ""};

  Configurable<int> sign_1{"sign_1", 1, ""};
  Configurable<int> particlePDG_1{"particlePDG_1", 2212, ""};
  Configurable<float> tpcNSigma_1{"tpcNSigma_1", 10.0, ""};
  Configurable<float> PIDtrshld_1{"PIDtrshld_1", 0.0, ""};
  Configurable<float> tpctofNSigma_1{"tpctofNSigma_1", 10.0, ""};

  Configurable<int> sign_2{"sign_2", 1, ""};
  Configurable<int> particlePDG_2{"particlePDG_2", 2212, ""};
  Configurable<float> tpcNSigma_2{"tpcNSigma_2", 10.0, ""};
  Configurable<float> PIDtrshld_2{"PIDtrshld_2", 0.0, ""};
  Configurable<float> tpctofNSigma_2{"tpctofNSigma_2", 10.0, ""};

  Configurable<float> deta{"deta", 0.01, ""};
  Configurable<float> dphi{"dphi", 0.01, ""};
  Configurable<float> radiusTPC{"radiusTPC", 1.2, ""};

  Configurable<int> multbinwidth{"multbinwidth", 50, ""};
  Configurable<int> vertexbinwidth{"vertexbinwidth", 2, ""};

  bool IsIdentical;

  std::map<std::string, float> trackcuts;
  std::map<std::string, double> PIDcuts_1;
  std::map<std::string, double> PIDcuts_2;

  void init(o2::framework::InitContext&)
  {
    registry.add("SE", "SE", kTH1F, {{500, 0.0, 5.0, "k*"}});
    registry.add("ME", "ME", kTH1F, {{500, 0.0, 5.0, "k*"}});

    IsIdentical = (particlePDG_1 == particlePDG_2) && (sign_1 == sign_2);

    trackcuts = {
      {"min_P", min_P}, {"max_P", max_P}, {"eta", eta}, {"tpcNClsFound", tpcNClsFound}, {"tpcChi2NCl", tpcChi2NCl}, {"dcaXY", dcaXY}, {"dcaZ", dcaXY}, {"itsNCls", itsNCls}};

    PIDcuts_1 = {
      {"sign", sign_1}, {"particlePDG", particlePDG_1}, {"tpcNSigma", tpcNSigma_1}, {"PIDtrshld", PIDtrshld_1}, {"tpctofNSigma", tpctofNSigma_1}};

    if (!IsIdentical) {
      PIDcuts_2 = {
        {"sign", sign_2}, {"particlePDG", particlePDG_2}, {"tpcNSigma", tpcNSigma_2}, {"PIDtrshld", PIDtrshld_2}, {"tpctofNSigma", tpctofNSigma_2}};
    }
  }

  void process(aod::SingleCollSel const& collisions, aod::SingleTrackSel const& tracks)
  {

    std::map<int, std::vector<FemtoParticle*>> selectedtracks_1;
    std::map<int, std::vector<FemtoParticle*>> selectedtracks_2;

    std::map<std::pair<int, int>, std::vector<int>> multbins;

    for (auto& track : tracks) {
      if (track.trackCuts(&trackcuts)) {

        if (track.pidCuts(&PIDcuts_1)) {
          // FemtoParticle* Particle1 = new FemtoParticle( track.energy(TDatabasePDG::Instance()->GetParticle(PIDcuts_1["particlePDG"])->Mass()), track.px(), track.py(), track.pz(), track.eta(), track.phi() );
          FemtoParticle* Particle1 = new FemtoParticle(track.energy(particle_mass(PIDcuts_1["particlePDG"])), track.px(), track.py(), track.pz(), track.eta(), track.phi());
          Particle1->SetSign(PIDcuts_1["sign"]);
          Particle1->SetMagField(0.5);
          selectedtracks_1[track.collisionId()].push_back(Particle1);
        }

        if (!IsIdentical && track.pidCuts(&PIDcuts_2)) {
          // FemtoParticle* Particle2 = new FemtoParticle( track.energy(TDatabasePDG::Instance()->GetParticle(PIDcuts_2["particlePDG"])->Mass()), track.px(), track.py(), track.pz(), track.eta(), track.phi() );
          FemtoParticle* Particle2 = new FemtoParticle(track.energy(particle_mass(PIDcuts_2["particlePDG"])), track.px(), track.py(), track.pz(), track.eta(), track.phi());
          Particle2->SetSign(PIDcuts_2["sign"]);
          Particle2->SetMagField(0.5);
          selectedtracks_2[track.collisionId()].push_back(Particle2);
        }
      }
    }

    for (auto& collision : collisions) {
      if (abs(collision.posZ()) > VertexZ)
        continue;
      if (selectedtracks_1.find(collision.globalIndex()) == selectedtracks_1.end()) {
        if (IsIdentical)
          continue;
        else if (selectedtracks_2.find(collision.globalIndex()) == selectedtracks_2.end())
          continue;
      }

      multbins[std::pair<int, int>{round(collision.posZ() / vertexbinwidth), floor(collision.mult() / multbinwidth)}].push_back(collision.globalIndex());
    }

    //============================================================================

    if (IsIdentical) {
      // same event identical

      for (auto i = selectedtracks_1.begin(); i != selectedtracks_1.end(); i++) {
        if ((i->second).size() < 2)
          continue;

        std::vector<std::vector<int>> comb_idx;
        comb(comb_idx, (i->second).size(), 2);

        for (auto indx : comb_idx) {

          FemtoPair* Pair = new FemtoPair((i->second)[indx[0]], (i->second)[indx[1]], IsIdentical);
          if (!Pair->IsClosePair(deta, dphi, radiusTPC))
            registry.fill(HIST("SE"), Pair->GetKstar()); // close pair separation doesn't work properly, need eta and phi at some position in the barrel
        }
      } // end of same event identical

      // mixed event identical

      for (auto i = multbins.begin(); i != multbins.end(); i++) {
        if ((i->second).size() < 2)
          continue;

        std::vector<std::vector<int>> comb_idx;
        comb(comb_idx, (i->second).size(), 2);

        for (auto indx : comb_idx) {
          int colId_1 = (i->second)[indx[0]];
          int colId_2 = (i->second)[indx[1]];

          for (auto ii : selectedtracks_1[colId_1]) {
            for (auto iii : selectedtracks_1[colId_2]) {
              FemtoPair* Pair = new FemtoPair(ii, iii, IsIdentical);
              registry.fill(HIST("ME"), Pair->GetKstar());
            }
          }
        }
      } // end of mixed event identical

    }

    else {

      // same event non-identical
      for (auto i = selectedtracks_1.begin(); i != selectedtracks_1.end(); i++) {

        auto j = selectedtracks_2.find(i->first);
        if (j == selectedtracks_2.end())
          continue;

        for (auto ii : i->second) {
          for (auto iii : j->second) {
            FemtoPair* Pair = new FemtoPair(ii, iii, IsIdentical);
            if (!Pair->IsClosePair(deta, dphi, radiusTPC))
              registry.fill(HIST("SE"), Pair->GetKstar()); // close pair separation doesn't work properly, need eta and phi at some position in the barrel
          }
        }
      } // end of same event non-identical

      for (auto i = multbins.begin(); i != multbins.end(); i++) {
        if ((i->second).size() < 2)
          continue;

        std::vector<std::vector<int>> comb_idx;
        comb(comb_idx, (i->second).size(), 2);

        for (auto indx : comb_idx) {
          int colId_1 = (i->second)[indx[0]];
          int colId_2 = (i->second)[indx[1]];

          for (auto ii : selectedtracks_1[colId_1]) {
            for (auto iii : selectedtracks_2[colId_2]) {
              FemtoPair* Pair = new FemtoPair(ii, iii, IsIdentical);
              registry.fill(HIST("ME"), Pair->GetKstar());
            }
          }
        }
      }
    }
  }
  // PROCESS_SWITCH(EtaPhiHistograms, processSelected, "process filtered track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FemtoCorrelations>(cfgc)};
}
