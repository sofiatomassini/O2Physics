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
/// \brief create a table applying some basic cuts on the ITS and DCA.
/// \author Sofia Tomassini
/// \since 31 May 2023

#include <Framework/AnalysisDataModel.h>
#include <fairlogger/Logger.h>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"

#include "PWGCF/DataModel/singletrackselector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;
//::singletrackselector; // the namespace defined in .h

struct singleTrackSelector {

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  // Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidEvTimeFlags, aod::TracksDCA,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                         aod::pidTPCFullDe, aod::pidTOFFullDe, aod::pidTOFbeta,
                         aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::FT0sCorrected>;

  Produces<o2::aod::SingleTrackSel> tableRow;

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (aod::evsel::sel8 == true));
  Filter vertexFilter = ((o2::aod::collision::posZ < 15.f) && (o2::aod::collision::posZ > -15.f));
  Filter trackFilter = ((o2::aod::track::itsChi2NCl <= 36.f) && (o2::aod::track::itsChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl <= 4.f));

  void process(soa::Filtered<Coll>::iterator const& collision, soa::Filtered<Trks> const& tracks)
  {

    tableRow.reserve(tracks.size());

    for (auto& track : tracks) {

      if (abs(track.tpcNSigmaPr()) > 4 || abs(track.tpcNSigmaDe()) > 4) {
        continue;
      } else {

        tableRow(track.collisionId(),
                 track.px(),
                 track.py(),
                 track.pz(),
                 track.p(),
                 track.pt(),
                 track.dcaXY(),
                 track.dcaZ(),
                 track.tpcNClsFound(),
                 track.tpcChi2NCl(),
                 // track.tpcSignal(),
                 // track.beta(),
                 track.itsNCls(),
                 track.itsChi2NCl(),
                 track.eta(),
                 track.phi(),
                 track.sign(),
                 //  track.tpcNClsCrossedRows(),
                 singletrackselector::packInTableInt<singletrackselector::storedcrossedrows::binning>(track.tpcNClsCrossedRows()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tofNSigmaPr()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tpcNSigmaPr()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tofNSigmaDe()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tpcNSigmaDe()));

        if ((singletrackselector::unPackInt<singletrackselector::storedcrossedrows::binning>(singletrackselector::packInTableInt<singletrackselector::storedcrossedrows::binning>(track.tpcNClsCrossedRows())) != track.tpcNClsCrossedRows()) && (singletrackselector::unPackInt<singletrackselector::storedcrossedrows::binning>(singletrackselector::packInTableInt<singletrackselector::storedcrossedrows::binning>(track.tpcNClsCrossedRows())) != (track.tpcNClsCrossedRows() - 1))) {
          LOG(info) << "crossedRows: " << track.tpcNClsCrossedRows()
                    << ", Unpacked crossedRows: " << singletrackselector::unPackInt<singletrackselector::storedcrossedrows::binning>(singletrackselector::packInTableInt<singletrackselector::storedcrossedrows::binning>(track.tpcNClsCrossedRows()));
          // << "deve essere uguale a" << singletrackselector::unPack<singletrackselector::nsigma::binning>;
        }
      }
    }

    // LOG(info) << "UnPacked crossedRows: " << singletrackselector::unPack(o2::aod::singletrackselector::storedcrossedrows::binning::binned_t(o2::aod::singletrackselector::storedcrossedrows::binning(tpcNClsCrossedRows)) )
    // LOG(info) << "Unpacked crossedRows: " << unpackedValue;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<singleTrackSelector>(cfgc)};
}
