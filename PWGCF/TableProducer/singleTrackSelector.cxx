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
using namespace o2::aod::singletrackselector; // the namespace defined in .h

struct singleTrackSelector {

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::TOFEvTime, aod::EvTimeTOFOnly, aod::TOFSignal, aod::pidEvTimeFlags,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                         aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::FT0sCorrected>;

  Produces<o2::aod::selectedtracks> tableRow;

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (o2::aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (o2::aod::evsel::sel8 == true));
  Filter vertexFilter = (o2::aod::‌collision::PosZ < 15.f && o2::aod::‌collision::PosZ > -15.f);
  Filter trackFilter = (o2::aod::track::ITSChi2NCl >= 0.f && o2::aod::track::ITSChi2NCl <= 36.f && o2::aod::track::TPCChi2NCl >= 0.f && o2::aod::track::TPCChi2NCl <= 4.f);

  void process(soa::Filtered<Coll>::iterator const& collision, soa::Filtered<Trks> const& tracks)
  {
    tableRow.reserve(tracks.size());

    for (auto& track : tracks) {

      tableRow(track.collisionId(),
               track.px(),
               track.py(),
               track.pz());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<singleTrackSelector>(cfgc)};
}
