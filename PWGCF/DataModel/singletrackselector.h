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
#include "Common/DataModel/PIDResponse.h"

namespace o2::aod
{

namespace singletrackselector
{
// Function to pack a float into a binned value in table
template <typename binningType>
typename binningType::binned_t packInTable(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return binningType::underflowBin;
  } else if (valueToBin >= binningType::binned_max) {
    return binningType::overflowBin;
  } else if (valueToBin >= 0) {
    return static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) + 0.5f);
  } else {
    return static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) - 0.5f);
  }
}

namespace storedcrossedrows
{
struct binning {
 public:
  typedef int8_t binned_t;
  static constexpr int nbins = 160;
  //(1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 159.5;
  static constexpr float binned_min = -0.5;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};
} // namespace storedcrossedrows

namespace nsigma
{
struct binning {
 public:
  typedef int8_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 6.35;
  static constexpr float binned_min = -6.35;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};
} // namespace nsigma

DECLARE_SOA_INDEX_COLUMN(Collision, collision); // Index to the collision
DECLARE_SOA_COLUMN(Px, px, float);              // Momentum of the track
DECLARE_SOA_COLUMN(Py, py, float);              // Momentum of the track
DECLARE_SOA_COLUMN(Pz, pz, float);              // Momentum of the track
DECLARE_SOA_COLUMN(P, p, float);                // Momentum of the track
DECLARE_SOA_COLUMN(Pt, pt, float);              // Momentum of the track
// DECLARE_SOA_COLUMN(PosZ, posZ, float);          // vertex position along z
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float); // impact parameter of the track
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);   // impact parameter of the track
// DECLARE_SOA_COLUMN(StoredCrossedRows, tpcNClsCrossedRows, float);   // impact parameter of the track
DECLARE_SOA_COLUMN(StoredCrossedRows, tpcNClsCrossedRows, storedcrossedrows::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr, storedtofNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaPr, storedtpcNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return nsigma::binning::bin_width * static_cast<float>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return nsigma::binning::bin_width * static_cast<float>(nsigma_binned); });

} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleTrackSel, "AOD", "STSEL", // Table of the variables for single track selection.
                  o2::soa::Index<>,
                  singletrackselector::CollisionId,
                  singletrackselector::Px,
                  singletrackselector::Py,
                  singletrackselector::Pz,
                  singletrackselector::P,
                  singletrackselector::Pt,
                  singletrackselector::DcaXY,
                  singletrackselector::DcaZ,
                  singletrackselector::StoredCrossedRows,
                  singletrackselector::StoredTOFNSigmaPr,
                  singletrackselector::StoredTPCNSigmaPr,
                  singletrackselector::TOFNSigmaPr<singletrackselector::StoredTOFNSigmaPr>,
                  singletrackselector::TPCNSigmaPr<singletrackselector::StoredTPCNSigmaPr>);
} // namespace o2::aod
#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
