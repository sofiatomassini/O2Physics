#ifndef COMMON_DATAMODEL_SINGLETRACKSELECTOR_H_
#define COMMON_DATAMODEL_SINGLETRACKSELECTOR_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod::singletrackselector {
    DECLARE_SOA_INDEX_COLUMN(Collision, collision); // Index to the collision
    DECLARE_SOA_COLUMN(Px, px, float); // Momentum of the track
    DECLARE_SOA_COLUMN(Py, py, float); // Momentum of the track
    DECLARE_SOA_COLUMN(Pz, pz, float); // Momentum of the track
}

DECLARE_SOA_TABLE(SingleTrackSel, "AOD", "STSEL", // Table of the delta TOF data format. One entry per track couple.
    o2::soa::Index<>,
    singletrackselector::CollisionId,
    singletrackselector::Px,
    singletrackselector::Py,
    singletrackselector::Pz);

#endif  


