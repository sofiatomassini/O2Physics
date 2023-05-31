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
/// \brief  create a table for single track selection.
/// \author Sofia Tomassini
/// \since 30 May 2023

#ifndef PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
#define PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace singletrackselector
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision); // Index to the collision
DECLARE_SOA_COLUMN(Px, px, float);              // Momentum of the track
DECLARE_SOA_COLUMN(Py, py, float);              // Momentum of the track
DECLARE_SOA_COLUMN(Pz, pz, float);              // Momentum of the track
} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleTrackSel, "AOD", "STSEL", // Table of the delta TOF data format. One entry per track couple.
                  o2::soa::Index<>,
                  singletrackselector::CollisionId,
                  singletrackselector::Px,
                  singletrackselector::Py,
                  singletrackselector::Pz);
} // namespace o2::aod
#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
