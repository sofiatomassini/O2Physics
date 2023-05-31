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

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct EtaPhiHistograms {
  using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if (doprocessStandard == true) {
      registry.add("eta", "eta", kTH1F, {{102, -2.01, 2.01, "eta"}});                                      //
      registry.add("phi", "phi", kTH1F, {{100, 0., 2. * M_PI, "phi"}});                                    //
      registry.add("pt", "pt", kTH1F, {{100, 0., 5.0, "pt"}});                                             //
      registry.add("NsigmaTPC", "NsigmaTPC", kTH2F, {{100, 0., 5.0, "pt"}, {100, -5., 5.0, "nsigmaTPC"}}); //
      registry.add("dEdxTPC", "dEdxTPC", kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}});
      registry.add("dEdxTPCinner", "dEdxTPCinner", kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}});
      registry.add("TOFSignal", "TOFSignal", kTH2F, {{200, 0., 5.0}, {100, 0., 1.5}});
      registry.add("NsigmaTOF_p", "NsigmaTOF_p", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("NsigmaTPCvsTOF", "n#sigma_{TPC}; n#sigma_{TOF}", kTH2F, {{100, -5., 5.}, {100, -5., 5.}}); // protons
      // registry.add("NsigmaTPCvsTOF^2", "n#sigma_{TPC}; n#sigma_{TOF}", kTH2F, {{100, 0., 2}, {100, 0., 2.}});  // protons^2
      registry.add("NsigmaTOF_K", "NsigmaTOF_K", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); // Kaons
      registry.add("NsigmaTPC_K", "NsigmaTPC_K", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); // Kaons
      registry.add("DCAxy", "#DCAxy", kTH1F, {{100, -0.2, 0.2}});
      registry.add("DCAz", "#DCAz", kTH1F, {{100, -0.2, 0.2}});
      registry.add("DCAxyvspt", "#DCAxyvspt", kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}});
      registry.add("DCAzvspt", "#DCAzvspt", kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}});
      registry.add("TPC_yess", "TPC_yess", kTH1F, {{100, 0., 20.}});
      registry.add("TPC_yess_eta", "TPC_yess_eta", kTH2F, {{100, 0., 20.0}, {200, -1., 1.}});
      //{"TOF_yess", "TOF_yess", {HistType::kTH1F, {{100, 0., 10.0}});
      registry.add("TPCandTOF", "TPCandTOF", kTH1F, {{100, 0., 20.0}});
      registry.add("TPCandTOF_eta", "TPCandTOF_eta", kTH2F, {{100, 0., 20.0}, {200, -1., 1.}});
      registry.add("hcrossedrows", ";track crossed rows;entries", kTH1F, {{160, -0.5, 159.5}});
      registry.add("charge_pt", "charge_pt", kTH2F, {{100, 0., 10.}, {10, -1.5, 1.5}});
      registry.add("charge", "charge", kTH1F, {{100, -1.5, 1.5}});
      registry.add("NsigmaTPC_d", "NsigmaTPC_d", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("NsigmaTOF_d", "NsigmaTOF_d", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("TOFSignalpr", ";TOFSignalpr", kTH2F, {{100, 0., 5.0}, {1000, 0., 1000.}});
      registry.add("TOFSignalK", "TOFSignalK", kTH2F, {{1000, 0., 5.0}, {1000, 0, 1000.}});
      registry.add("TOFSignald", "TOFSignald", kTH2F, {{1000, 0., 5.0}, {1000, 0., 1000.}});
      //{"separation_pk", "separation_pk", {HistType::kTH2F,{{100, 0., 5.}, {}}}}
      //
      // definisco di nuovo gli stessi istogrammi ma _cut
      registry.add("eta_cut", "#eta", kTH1F, {{102, -2.01, 2.01}});                             //
      registry.add("phi_cut", "#varphi", kTH1F, {{100, 0., 2. * M_PI}});                        //
      registry.add("pt_cut", "pt", kTH1F, {{100, 0., 5.0}});                                    //
      registry.add("NsigmaTPC_cut", "NsigmaTPC_cut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); //
      registry.add("dEdxTPC_cut", "dEdxTPC_cut", kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}});
      registry.add("dEdxTPCinner_cut", "dEdxTPCinner_cut", kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}});
      registry.add("TOFSignal_cut", "TOFSignal_cut", kTH2F, {{200, 0., 5.0}, {100, 0., 1.5}});
      registry.add("NsigmaTOF_p_cut", "NsigmaTOF_p_cut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}),
        registry.add("NsigmaTPCvsTOF_cut", "n#sigma_{TPC}; n#sigma_{TOF}", kTH2F, {{100, -5., 5.}, {100, -5., 5.}}); // protons
      // registry.add("NsigmaTPCvsTOF^2_cut", "n#sigma_{TPC}; n#sigma_{TOF}", kTH2F, {{100, 0., 2.}, {100, 0., 2.}});   // protons^2
      registry.add("NsigmaTOF_K_cut", "NsigmaTOF_K_cut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); // Kaons
      registry.add("NsigmaTPC_K_cut", "NsigmaTPC_K_cut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); // Kaons
      registry.add("DCAxy_cut", "#DCAxy_cut", kTH1F, {{100, -0.2, 0.2}});
      registry.add("DCAz_cut", "#DCAz_cut", kTH1F, {{100, -0.2, 0.2}});
      registry.add("DCAxyvspt_cut", "#DCAxyvspt_cut", kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}});
      registry.add("DCAzvspt_cut", "#DCAzvspt_cut", kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}});
      registry.add("TPC_yess_cut", "TPC_yess_cut", kTH1F, {{100, 0., 20.}});
      registry.add("TPC_yess_eta_cut", "TPC_yess_eta_cut", kTH2F, {{100, 0., 20.0}, {200, -1., 1.}});
      //{"TOF_yess", "TOF_yess", kTH1F, {{100, 0., 10.0}});
      registry.add("TPCandTOF_cut", "TPCandTOF_cut", kTH1F, {{100, 0., 20.0}});
      registry.add("TPCandTOF_eta_cut", "TPCandTOF_eta_cut", kTH2F, {{100, 0., 20.0}, {200, -1., 1.}});
      registry.add("hcrossedrows_cut", ";track crossed rows;entries", kTH1F, {{160, -0.5, 159.5}});
      registry.add("charge_pt_cut", "charge_pt_cut", kTH2F, {{100, 0., 10.}, {10, -1.5, 1.5}});
      registry.add("charge_cut", "charge_cut", kTH1F, {{100, -1.5, 1.5}});
      registry.add("NsigmaTPC_d_cut", "NsigmaTPC_d_cut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("NsigmaTOF_d_cut", "NsigmaTOF_d_cut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("TOFSignalpr_cut", ";TOFSignalpr_cut", kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}});
      registry.add("TOFSignalK_cut", "TOFSignalK_cut", kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}});
      registry.add("TOFSignald_cut", "TOFSignald_cut", kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}});
      ////standard cuts histos

      registry.add("eta_stdcut", "#eta_stdcut", kTH1F, {{102, -2.01, 2.01}});                         //
      registry.add("phi_stdcut", "#varphi_stdcut", kTH1F, {{100, 0., 2. * M_PI}});                    //
      registry.add("pt_stdcut", "pt_stdcut", kTH1F, {{100, 0., 5.0}});                                //
      registry.add("NsigmaTPC_stdcut", "NsigmaTPC_stdcut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); //
      registry.add("dEdxTPC_stdcut", "dEdxTPC_stdcut", kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}});
      registry.add("dEdxTPCinner_stdcut", "dEdxTPCinner_stdcut", kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}});
      registry.add("TOFSignal_stdcut", "TOFSignal_stdcut", kTH2F, {{200, 0., 5.0}, {100, 0., 1.5}});
      registry.add("NsigmaTOF_p_stdcut", "NsigmaTOF_p_stdcut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("NsigmaTPCvsTOF_stdcut", "n#sigma_{TPC}; n#sigma_{TOF}", kTH2F, {{100, -5., 5.}, {100, -5., 5.}}); // protons
      // registry.add("NsigmaTPCvsTOF^2_stdcut", "n#sigma_{TPC}; n#sigma_{TOF}", kTH2F, {{100, 0., 5.}, {100, 0., 5.}}); // protons^2
      registry.add("NsigmaTOF_K_stdcut", "NsigmaTOF_K_stdcut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); // Kaons
      registry.add("NsigmaTPC_K_stdcut", "NsigmaTPC_K_stdcut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}); // Kaons
      registry.add("DCAxy_stdcut", "#DCAxy_stdcut", kTH1F, {{100, -0.2, 0.2}});
      registry.add("DCAz_stdcut", "#DCAz_stdcut", kTH1F, {{100, -0.2, 0.2}});
      registry.add("DCAxyvspt_stdcut", "#DCAxyvspt_stdcut", kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}});
      registry.add("DCAzvspt_stdcut", "#DCAzvspt_stdcut", kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}});
      registry.add("TPC_yess_stdcut", "TPC_yess_stdcut", kTH1F, {{100, 0., 20.}});
      registry.add("TPC_yess_eta_stdcut", "TPC_yess_eta_stdcut", kTH2F, {{100, 0., 20.0}, {200, -1., 1.}});
      //{"TOF_yess", "TOF_yess", kTH1F, {{100, 0., 10.0}});
      registry.add("TPCandTOF_stdcut", "TPCandTOF_stdcut", kTH1F, {{100, 0., 20.0}});
      registry.add("TPCandTOF_eta_stdcut", "TPCandTOF_eta_stdcut", kTH2F, {{100, 0., 20.0}, {200, -1., 1.}});
      registry.add("hcrossedrows_stdcut", ";track crossed rows;entries", kTH1F, {{160, -0.5, 159.5}});
      registry.add("charge_pt_stdcut", "charge_pt_stdcut", kTH2F, {{100, 0., 10.}, {10, -1.5, 1.5}});
      registry.add("charge_stdcut", "charge_stdcut", kTH1F, {{100, -1.5, 1.5}});
      registry.add("NsigmaTPC_d_stdcut", "NsigmaTPC_d_stdcut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("NsigmaTOF_d_stdcut", "NsigmaTOF_d_stdcut", kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}});
      registry.add("TOFSignalpr_stdcut", ";TOFSignalpr_stdcut", kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}});
      registry.add("TOFSignalK_stdcut", "TOFSignalK_stdcut", kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}});
      registry.add("TOFSignald_stdcut", "TOFSignald_stdcut", kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}});
    }
    if (doprocessSelected == true) {
      registry.add("px", "px", kTH1F, {{100, 0., 5., "px"}});
      registry.add("py", "py", kTH1F, {{100, 0., 5., "py"}});
      registry.add("pz", "pz", kTH1F, {{100, 0., 5., "pz"}});
    }
  } //

  void processStandard(allinfo const& tracks)
  {
    for (auto& track : tracks) {
      registry.fill(HIST("eta"), track.eta());
      registry.fill(HIST("phi"), track.phi());
      registry.fill(HIST("pt"), track.pt());
      registry.fill(HIST("dEdxTPC"), track.pt(), track.tpcSignal());
      registry.fill(HIST("dEdxTPCinner"), track.tpcInnerParam(), track.tpcSignal());
      registry.fill(HIST("NsigmaTPC"), track.pt(), track.tpcNSigmaPr()); // tpcNSigmaStorePr
      registry.fill(HIST("TOFSignal"), track.pt(), track.beta());
      registry.fill(HIST("NsigmaTOF_p"), track.pt(), track.tofNSigmaPr());
      registry.fill(HIST("NsigmaTPCvsTOF"), track.tpcNSigmaPr(), track.tofNSigmaPr());
      // registry.fill(HIST("NsigmaTPCvsTOF^2"), track.tpcNSigmaPr() * track.tpcNSigmaPr(), track.tofNSigmaPr() * track.tofNSigmaPr());
      registry.fill(HIST("NsigmaTOF_K"), track.pt(), track.tofNSigmaKa());
      registry.fill(HIST("NsigmaTPC_K"), track.pt(), track.tpcNSigmaKa());
      registry.fill(HIST("DCAxy"), track.dcaXY());
      registry.fill(HIST("DCAz"), track.dcaZ());
      registry.fill(HIST("DCAxyvspt"), track.pt(), track.dcaXY());
      registry.fill(HIST("DCAzvspt"), track.pt(), track.dcaZ());
      registry.fill(HIST("hcrossedrows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("charge_pt"), track.pt(), track.sign());
      registry.fill(HIST("charge"), track.sign());
      registry.fill(HIST("NsigmaTOF_d"), track.pt(), track.tofNSigmaDe());
      registry.fill(HIST("NsigmaTPC_d"), track.pt(), track.tpcNSigmaDe());
      // registry.fill(HIST("TOF_yess"), track.pt(), track.hasTOF());
      // registry.fill(HIST("TOFSignalK"), track.pt(), track.tofExpSignalKa());
      registry.fill(HIST("TOFSignalpr"), track.pt(), track.tofExpSignalPr(track.tofSignal()));
      registry.fill(HIST("TOFSignalK"), track.pt(), track.tofExpSignalKa(track.tofSignal()));
      registry.fill(HIST("TOFSignald"), track.pt(), track.tofExpSignalDe(track.tofSignal()));

      if (track.hasTPC()) {
        registry.fill(HIST("TPC_yess"), track.pt());
        registry.fill(HIST("TPC_yess_eta"), track.pt(), track.eta());
        if (track.hasTOF()) {
          registry.fill(HIST("TPCandTOF"), track.pt());
          registry.fill(HIST("TPCandTOF_eta"), track.pt(), track.eta());
        }
      }

      // sempre sulle tracce faccio il fill degli istogrammi dopo i cut a pt>150 MEV/c e |eta|<0.8
      if (track.pt() > 0.15 && abs(track.eta()) < 0.8) {

        registry.fill(HIST("eta_cut"), track.eta());
        registry.fill(HIST("phi_cut"), track.phi());
        registry.fill(HIST("pt_cut"), track.pt());
        registry.fill(HIST("dEdxTPC_cut"), track.pt(), track.tpcSignal());
        registry.fill(HIST("dEdxTPCinner_cut"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("NsigmaTPC_cut"), track.pt(), track.tpcNSigmaPr()); // tpcNSigmaStorePr
        registry.fill(HIST("TOFSignal_cut"), track.pt(), track.beta());
        registry.fill(HIST("NsigmaTOF_p_cut"), track.pt(), track.tofNSigmaPr());
        registry.fill(HIST("NsigmaTPCvsTOF_cut"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        // registry.fill(HIST("NsigmaTPCvsTOF^2_cut"), track.tpcNSigmaPr() * track.tpcNSigmaPr(), track.tofNSigmaPr() * track.tofNSigmaPr());
        registry.fill(HIST("NsigmaTOF_K_cut"), track.pt(), track.tofNSigmaKa());
        registry.fill(HIST("NsigmaTPC_K_cut"), track.pt(), track.tpcNSigmaKa());
        registry.fill(HIST("DCAxy_cut"), track.dcaXY());
        registry.fill(HIST("DCAz_cut"), track.dcaZ());
        registry.fill(HIST("DCAxyvspt_cut"), track.pt(), track.dcaXY());
        registry.fill(HIST("DCAzvspt_cut"), track.pt(), track.dcaZ());
        registry.fill(HIST("hcrossedrows_cut"), track.tpcNClsCrossedRows());
        registry.fill(HIST("charge_pt_cut"), track.pt(), track.sign());

        registry.fill(HIST("charge_cut"), track.sign());
        registry.fill(HIST("NsigmaTOF_d_cut"), track.pt(), track.tofNSigmaDe());
        registry.fill(HIST("NsigmaTPC_d_cut"), track.pt(), track.tpcNSigmaDe());
        registry.fill(HIST("TOFSignalpr_cut"), track.pt(), track.tofExpSignalPr(track.tofSignal()));
        registry.fill(HIST("TOFSignalK_cut"), track.pt(), track.tofExpSignalKa(track.tofSignal()));
        registry.fill(HIST("TOFSignald_cut"), track.pt(), track.tofExpSignalDe(track.tofSignal()));

        if (track.hasTPC()) {
          registry.fill(HIST("TPC_yess_cut"), track.pt());
          registry.fill(HIST("TPC_yess_eta_cut"), track.pt(), track.eta());
          if (track.hasTOF()) {
            registry.fill(HIST("TPCandTOF_cut"), track.pt());
            registry.fill(HIST("TPCandTOF_eta_cut"), track.pt(), track.eta());
          }
        }
      }

      /// filling histos after standard cuts

      if (track.isPrimaryTrack() && (track.pt() > 0.15 && abs(track.eta()) < 0.8)) {
        registry.fill(HIST("eta_stdcut"), track.eta());
        registry.fill(HIST("phi_stdcut"), track.phi());
        registry.fill(HIST("pt_stdcut"), track.pt());
        registry.fill(HIST("dEdxTPC_stdcut"), track.pt(), track.tpcSignal());
        registry.fill(HIST("dEdxTPCinner_stdcut"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("NsigmaTPC_stdcut"), track.pt(), track.tpcNSigmaPr()); // tpcNSigmaStorePr
        registry.fill(HIST("TOFSignal_stdcut"), track.pt(), track.beta());
        registry.fill(HIST("NsigmaTOF_p_stdcut"), track.pt(), track.tofNSigmaPr());
        registry.fill(HIST("NsigmaTPCvsTOF_stdcut"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        // registry.fill(HIST("NsigmaTPCvsTOF^2_stdcut"),track.tpcNSigmaPr() * track.tpcNSigmaPr(), track.tofNSigmaPr() * track.tofNSigmaPr());
        registry.fill(HIST("NsigmaTOF_K_stdcut"), track.pt(), track.tofNSigmaKa());
        registry.fill(HIST("NsigmaTPC_K_stdcut"), track.pt(), track.tpcNSigmaKa());
        registry.fill(HIST("DCAxy_stdcut"), track.dcaXY());
        registry.fill(HIST("DCAz_stdcut"), track.dcaZ());
        registry.fill(HIST("DCAxyvspt_stdcut"), track.pt(), track.dcaXY());
        registry.fill(HIST("DCAzvspt_stdcut"), track.pt(), track.dcaZ());
        registry.fill(HIST("hcrossedrows_stdcut"), track.tpcNClsCrossedRows());
        registry.fill(HIST("charge_pt_stdcut"), track.pt(), track.sign());

        registry.fill(HIST("charge_stdcut"), track.sign());
        registry.fill(HIST("NsigmaTOF_d_stdcut"), track.pt(), track.tofNSigmaDe());
        registry.fill(HIST("NsigmaTPC_d_stdcut"), track.pt(), track.tpcNSigmaDe());
        registry.fill(HIST("TOFSignalpr_stdcut"), track.pt(), track.tofExpSignalPr(track.tofSignal()));
        registry.fill(HIST("TOFSignalK_stdcut"), track.pt(), track.tofExpSignalKa(track.tofSignal()));
        registry.fill(HIST("TOFSignald_stdcut"), track.pt(), track.tofExpSignalDe(track.tofSignal()));

        if (track.hasTPC()) {
          registry.fill(HIST("TPC_yess_stdcut"), track.pt());
          registry.fill(HIST("TPC_yess_eta_stdcut"), track.pt(), track.eta());
          if (track.hasTOF()) {
            registry.fill(HIST("TPCandTOF_stdcut"), track.pt());
            registry.fill(HIST("TPCandTOF_eta_stdcut"), track.pt(), track.eta());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(EtaPhiHistograms, processStandard, "process non filtered track", false);

  void processSelected(aod::SingleTrackSel const& tracks)
  {
    for (auto& track : tracks) {

      registry.fill(HIST("px"), track.px());
      registry.fill(HIST("py"), track.py());
      registry.fill(HIST("pz"), track.pz());
    }
  }
  PROCESS_SWITCH(EtaPhiHistograms, processSelected, "process filtered track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EtaPhiHistograms>(cfgc)};
}
