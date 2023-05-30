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
/// \brief Use the HistogramRegistry to manipulate histograms.
/// \author
/// \since

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

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;

struct EtaPhiHistograms {
  using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe >; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{
    "registry",
    {
      {"eta", "#eta", {HistType::kTH1F, {{102, -2.01, 2.01}}}},     //
      {"phi", "#varphi", {HistType::kTH1F, {{100, 0., 2. * M_PI}}}}, //
      {"pt", "pt", {HistType::kTH1F, {{100, 0., 5.0}}}}, //
      {"NsigmaTPC", "NsigmaTPC", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}}, //
      {"dEdxTPC", "dEdxTPC", {HistType::kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}}}},
      {"dEdxTPCinner", "dEdxTPCinner", {HistType::kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}}}},
      {"TOFSignal", "TOFSignal", {HistType::kTH2F, {{200, 0., 5.0}, {100, 0., 1.5}}}},
      {"NsigmaTOF_p", "NsigmaTOF_p", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"NsigmaTPCvsTOF","n#sigma_{TPC}; n#sigma_{TOF}", {HistType::kTH2F, {{100,-5., 5.}, {100, -5., 5.}}}}, //protons
      {"NsigmaTPCvsTOF^2","n#sigma_{TPC}; n#sigma_{TOF}", {HistType::kTH2F, {{100, 0., 2}, {100, 0., 2.}}}}, //protons^2
      {"NsigmaTOF_K", "NsigmaTOF_K", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},//Kaons
      {"NsigmaTPC_K", "NsigmaTPC_K", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}}, //Kaons
      {"DCAxy", "#DCAxy", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
      {"DCAz", "#DCAz", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
      {"DCAxyvspt", "#DCAxyvspt", {HistType::kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}}}},
      {"DCAzvspt", "#DCAzvspt", {HistType::kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}}}},
      {"TPC_yess", "TPC_yess", {HistType::kTH1F, {{100, 0., 20.}}}},
      {"TPC_yess_eta", "TPC_yess_eta", {HistType::kTH2F, {{100, 0., 20.0}, {200, -1., 1.}}}},
      //{"TOF_yess", "TOF_yess", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"TPCandTOF", "TPCandTOF", {HistType::kTH1F, {{100, 0., 20.0}}}},
      {"TPCandTOF_eta", "TPCandTOF_eta", {HistType::kTH2F, {{100, 0., 20.0}, {200, -1., 1.}}}},
      {"hcrossedrows", ";track crossed rows;entries", {HistType::kTH1F, {{160, -0.5, 159.5}}}},
      {"charge_pt", "charge_pt", {HistType::kTH2F, {{100, 0., 10.}, {10, -1.5, 1.5}}}},
      {"charge", "charge", {HistType::kTH1F, {{100, -1.5, 1.5 }}}},
      {"NsigmaTPC_d", "NsigmaTPC_d",{HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"NsigmaTOF_d", "NsigmaTOF_d",{HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"TOFSignalpr", ";TOFSignalpr", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 1000.}}}},
      {"TOFSignalK", "TOFSignalK", {HistType::kTH2F, {{1000, 0., 5.0}, {1000, 0, 1000.}}}},
      {"TOFSignald","TOFSignald",{HistType::kTH2F, {{1000, 0., 5.0}, {1000, 0., 1000.}}}},
      //{"separation_pk", "separation_pk", {HistType::kTH2F,{{100, 0., 5.}, {}}}}
      //
      // definisco di nuovo gli stessi istogrammi ma _cut
      {"eta_cut", "#eta", {HistType::kTH1F, {{102, -2.01, 2.01}}}},     //
      {"phi_cut", "#varphi", {HistType::kTH1F, {{100, 0., 2. * M_PI}}}}, //
      {"pt_cut", "pt", {HistType::kTH1F, {{100, 0., 5.0}}}}, //
      {"NsigmaTPC_cut", "NsigmaTPC_cut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}}, //
      {"dEdxTPC_cut", "dEdxTPC_cut", {HistType::kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}}}},
      {"dEdxTPCinner_cut", "dEdxTPCinner_cut", {HistType::kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}}}},
      {"TOFSignal_cut", "TOFSignal_cut", {HistType::kTH2F, {{200, 0., 5.0}, {100, 0., 1.5}}}},
      {"NsigmaTOF_p_cut", "NsigmaTOF_p_cut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"NsigmaTPCvsTOF_cut","n#sigma_{TPC}; n#sigma_{TOF}", {HistType::kTH2F, {{100,-5., 5.}, {100, -5., 5.}}}}, //protons
      {"NsigmaTPCvsTOF^2_cut","n#sigma_{TPC}; n#sigma_{TOF}", {HistType::kTH2F, {{100, 0., 2.}, {100, 0., 2.}}}}, //protons^2
      {"NsigmaTOF_K_cut", "NsigmaTOF_K_cut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},//Kaons
      {"NsigmaTPC_K_cut", "NsigmaTPC_K_cut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}}, //Kaons
      {"DCAxy_cut", "#DCAxy_cut", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
      {"DCAz_cut", "#DCAz_cut", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
      {"DCAxyvspt_cut", "#DCAxyvspt_cut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}}}},
      {"DCAzvspt_cut", "#DCAzvspt_cut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}}}},
      {"TPC_yess_cut", "TPC_yess_cut", {HistType::kTH1F, {{100, 0., 20.}}}},
      {"TPC_yess_eta_cut", "TPC_yess_eta_cut", {HistType::kTH2F, {{100, 0., 20.0}, {200, -1., 1.}}}},
      //{"TOF_yess", "TOF_yess", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"TPCandTOF_cut", "TPCandTOF_cut", {HistType::kTH1F, {{100, 0., 20.0}}}},
      {"TPCandTOF_eta_cut", "TPCandTOF_eta_cut", {HistType::kTH2F, {{100, 0., 20.0}, {200, -1., 1.}}}},
      {"hcrossedrows_cut", ";track crossed rows;entries", {HistType::kTH1F, {{160, -0.5, 159.5}}}},
      {"charge_pt_cut", "charge_pt_cut", {HistType::kTH2F, {{100, 0., 10.}, {10, -1.5, 1.5}}}},
      {"charge_cut", "charge_cut", {HistType::kTH1F, {{100, -1.5, 1.5 }}}},
      {"NsigmaTPC_d_cut", "NsigmaTPC_d_cut",{HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"NsigmaTOF_d_cut", "NsigmaTOF_d_cut",{HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"TOFSignalpr_cut", ";TOFSignalpr_cut", {HistType::kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}}}},
      {"TOFSignalK_cut", "TOFSignalK_cut", {HistType::kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}}}},
      {"TOFSignald_cut","TOFSignald_cut",{HistType::kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}}}},
      ////standard cuts histos

      {"eta_stdcut", "#eta_stdcut", {HistType::kTH1F, {{102, -2.01, 2.01}}}},     //
      {"phi_stdcut", "#varphi_stdcut", {HistType::kTH1F, {{100, 0., 2. * M_PI}}}}, //
      {"pt_stdcut", "pt_stdcut", {HistType::kTH1F, {{100, 0., 5.0}}}}, //
      {"NsigmaTPC_stdcut", "NsigmaTPC_stdcut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}}, //
      {"dEdxTPC_stdcut", "dEdxTPC_stdcut", {HistType::kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}}}},
      {"dEdxTPCinner_stdcut", "dEdxTPCinner_stdcut", {HistType::kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}}}},
      {"TOFSignal_stdcut", "TOFSignal_stdcut", {HistType::kTH2F, {{200, 0., 5.0}, {100, 0., 1.5}}}},
      {"NsigmaTOF_p_stdcut", "NsigmaTOF_p_stdcut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"NsigmaTPCvsTOF_stdcut","n#sigma_{TPC}; n#sigma_{TOF}", {HistType::kTH2F, {{100,-5., 5.}, {100, -5., 5.}}}}, //protons
      {"NsigmaTPCvsTOF^2_stdcut","n#sigma_{TPC}; n#sigma_{TOF}", {HistType::kTH2F, {{100, 0., 5.}, {100, 0., 5.}}}}, //protons^2
      {"NsigmaTOF_K_stdcut", "NsigmaTOF_K_stdcut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},//Kaons
      {"NsigmaTPC_K_stdcut", "NsigmaTPC_K_stdcut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}}, //Kaons
      {"DCAxy_stdcut", "#DCAxy_stdcut", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
      {"DCAz_stdcut", "#DCAz_stdcut", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
      {"DCAxyvspt_stdcut", "#DCAxyvspt_stdcut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}}}},
      {"DCAzvspt_stdcut", "#DCAzvspt_stdcut", {HistType::kTH2F, {{100, 0., 5.0}, {100, -0.2, 0.2}}}},
      {"TPC_yess_stdcut", "TPC_yess_stdcut", {HistType::kTH1F, {{100, 0., 20.}}}},
      {"TPC_yess_eta_stdcut", "TPC_yess_eta_stdcut", {HistType::kTH2F, {{100, 0., 20.0}, {200, -1., 1.}}}},
      //{"TOF_yess", "TOF_yess", {HistType::kTH1F, {{100, 0., 10.0}}}},
      {"TPCandTOF_stdcut", "TPCandTOF_stdcut", {HistType::kTH1F, {{100, 0., 20.0}}}},
      {"TPCandTOF_eta_stdcut", "TPCandTOF_eta_stdcut", {HistType::kTH2F, {{100, 0., 20.0}, {200, -1., 1.}}}},
      {"hcrossedrows_stdcut", ";track crossed rows;entries", {HistType::kTH1F, {{160, -0.5, 159.5}}}},
      {"charge_pt_stdcut", "charge_pt_stdcut", {HistType::kTH2F, {{100, 0., 10.}, {10, -1.5, 1.5}}}},
      {"charge_stdcut", "charge_stdcut", {HistType::kTH1F, {{100, -1.5, 1.5 }}}},
      {"NsigmaTPC_d_stdcut", "NsigmaTPC_d_stdcut",{HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"NsigmaTOF_d_stdcut", "NsigmaTOF_d_stdcut",{HistType::kTH2F, {{100, 0., 5.0}, {100, -5., 5.0}}}},
      {"TOFSignalpr_stdcut", ";TOFSignalpr_stdcut", {HistType::kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}}}},
      {"TOFSignalK_stdcut", "TOFSignalK_stdcut", {HistType::kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}}}},
      {"TOFSignald_stdcut","TOFSignald_stdcut",{HistType::kTH2F, {{1000, 0., 5.0}, {200, 0., 200.}}}}

    }                                                               //
  };

  void process(aod::Collision const&, allinfo const& tracks)
  {
    for (auto& track : tracks) {
      registry.get<TH1>(HIST("eta"))->Fill(track.eta());
      registry.get<TH1>(HIST("phi"))->Fill(track.phi());
      registry.get<TH1>(HIST("pt"))->Fill(track.pt());
      registry.get<TH2>(HIST("dEdxTPC"))->Fill(track.pt(), track.tpcSignal());
      registry.get<TH2>(HIST("dEdxTPCinner"))->Fill(track.tpcInnerParam(), track.tpcSignal());
      registry.get<TH2>(HIST("NsigmaTPC"))->Fill(track.pt(), track.tpcNSigmaPr()); // tpcNSigmaStorePr
      registry.get<TH2>(HIST("TOFSignal"))->Fill(track.pt(), track.beta());
      registry.get<TH2>(HIST("NsigmaTOF_p"))->Fill(track.pt(), track.tofNSigmaPr());
      registry.get<TH2>(HIST("NsigmaTPCvsTOF"))->Fill(track.tpcNSigmaPr(), track.tofNSigmaPr());
      registry.get<TH2>(HIST("NsigmaTPCvsTOF^2"))->Fill(track.tpcNSigmaPr()*track.tpcNSigmaPr(), track.tofNSigmaPr()*track.tofNSigmaPr());
      registry.get<TH2>(HIST("NsigmaTOF_K"))->Fill(track.pt(), track.tofNSigmaKa());
      registry.get<TH2>(HIST("NsigmaTPC_K"))->Fill(track.pt(), track.tpcNSigmaKa());
      registry.get<TH1>(HIST("DCAxy"))->Fill(track.dcaXY());
      registry.get<TH1>(HIST("DCAz"))->Fill(track.dcaZ());
      registry.get<TH2>(HIST("DCAxyvspt"))->Fill(track.pt(), track.dcaXY());
      registry.get<TH2>(HIST("DCAzvspt"))->Fill(track.pt(), track.dcaZ());
      registry.fill(HIST("hcrossedrows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("charge_pt"), track.pt(), track.sign());

      registry.fill(HIST("charge"),track.sign());
      registry.get<TH2>(HIST("NsigmaTOF_d"))->Fill(track.pt(), track.tofNSigmaDe());
      registry.get<TH2>(HIST("NsigmaTPC_d"))->Fill(track.pt(), track.tpcNSigmaDe());
      //registry.fill(HIST("TOF_yess"), track.pt(), track.hasTOF());
      //registry.fill(HIST("TOFSignalK"), track.pt(), track.tofExpSignalKa());
      registry.fill(HIST("TOFSignalpr"), track.pt(), track.tofExpSignalPr(track.tofSignal()));
      registry.fill(HIST("TOFSignalK"), track.pt(), track.tofExpSignalKa(track.tofSignal()));
      registry.fill(HIST("TOFSignald"), track.pt(), track.tofExpSignalDe(track.tofSignal()));


        if (track.hasTPC()) {
          registry.fill(HIST("TPC_yess"), track.pt());
          registry.fill(HIST("TPC_yess_eta"), track.pt(),track.eta());
          if (track.hasTOF()){registry.fill(HIST("TPCandTOF"), track.pt());
                              registry.fill(HIST("TPCandTOF_eta"), track.pt(),track.eta());
                              }
          }

        

// sempre sulle tracce faccio il fill degli istogrammi dopo i cut a pt>150 MEV/c e |eta|<0.8
        if (track.pt()>0.15 && abs(track.eta())<0.8)  {

          registry.get<TH1>(HIST("eta_cut"))->Fill(track.eta());
          registry.get<TH1>(HIST("phi_cut"))->Fill(track.phi());
          registry.get<TH1>(HIST("pt_cut"))->Fill(track.pt());
          registry.get<TH2>(HIST("dEdxTPC_cut"))->Fill(track.pt(), track.tpcSignal());
          registry.get<TH2>(HIST("dEdxTPCinner_cut"))->Fill(track.tpcInnerParam(), track.tpcSignal());
          registry.get<TH2>(HIST("NsigmaTPC_cut"))->Fill(track.pt(), track.tpcNSigmaPr()); // tpcNSigmaStorePr
          registry.get<TH2>(HIST("TOFSignal_cut"))->Fill(track.pt(), track.beta());
          registry.get<TH2>(HIST("NsigmaTOF_p_cut"))->Fill(track.pt(), track.tofNSigmaPr());
          registry.get<TH2>(HIST("NsigmaTPCvsTOF_cut"))->Fill(track.tpcNSigmaPr(), track.tofNSigmaPr());
          registry.get<TH2>(HIST("NsigmaTPCvsTOF^2_cut"))->Fill(track.tpcNSigmaPr()*track.tpcNSigmaPr(), track.tofNSigmaPr()*track.tofNSigmaPr());
          registry.get<TH2>(HIST("NsigmaTOF_K_cut"))->Fill(track.pt(), track.tofNSigmaKa());
          registry.get<TH2>(HIST("NsigmaTPC_K_cut"))->Fill(track.pt(), track.tpcNSigmaKa());
          registry.get<TH1>(HIST("DCAxy_cut"))->Fill(track.dcaXY());
          registry.get<TH1>(HIST("DCAz_cut"))->Fill(track.dcaZ());
          registry.get<TH2>(HIST("DCAxyvspt_cut"))->Fill(track.pt(), track.dcaXY());
          registry.get<TH2>(HIST("DCAzvspt_cut"))->Fill(track.pt(), track.dcaZ());
          registry.fill(HIST("hcrossedrows_cut"), track.tpcNClsCrossedRows());
          registry.fill(HIST("charge_pt_cut"), track.pt(), track.sign());

          registry.fill(HIST("charge_cut"),track.sign());
          registry.get<TH2>(HIST("NsigmaTOF_d_cut"))->Fill(track.pt(), track.tofNSigmaDe());
          registry.get<TH2>(HIST("NsigmaTPC_d_cut"))->Fill(track.pt(), track.tpcNSigmaDe());
          registry.fill(HIST("TOFSignalpr_cut"), track.pt(), track.tofExpSignalPr(track.tofSignal()));
          registry.fill(HIST("TOFSignalK_cut"), track.pt(), track.tofExpSignalKa(track.tofSignal()));
          registry.fill(HIST("TOFSignald_cut"), track.pt(), track.tofExpSignalDe(track.tofSignal()));

          

          if (track.hasTPC()) {
          registry.fill(HIST("TPC_yess_cut"), track.pt());
          registry.fill(HIST("TPC_yess_eta_cut"), track.pt(),track.eta());
          if (track.hasTOF()){registry.fill(HIST("TPCandTOF_cut"), track.pt());
                              registry.fill(HIST("TPCandTOF_eta_cut"), track.pt(),track.eta());
                              }
          }
          
          }


          ///filling histos after standard cuts

          if(track.isPrimaryTrack() && (track.pt()>0.15 && abs(track.eta())<0.8)){
              registry.get<TH1>(HIST("eta_stdcut"))->Fill(track.eta());
              registry.get<TH1>(HIST("phi_stdcut"))->Fill(track.phi());
              registry.get<TH1>(HIST("pt_stdcut"))->Fill(track.pt());
              registry.get<TH2>(HIST("dEdxTPC_stdcut"))->Fill(track.pt(), track.tpcSignal());
              registry.get<TH2>(HIST("dEdxTPCinner_stdcut"))->Fill(track.tpcInnerParam(), track.tpcSignal());
              registry.get<TH2>(HIST("NsigmaTPC_stdcut"))->Fill(track.pt(), track.tpcNSigmaPr()); // tpcNSigmaStorePr
              registry.get<TH2>(HIST("TOFSignal_stdcut"))->Fill(track.pt(), track.beta());
              registry.get<TH2>(HIST("NsigmaTOF_p_stdcut"))->Fill(track.pt(), track.tofNSigmaPr());
              registry.get<TH2>(HIST("NsigmaTPCvsTOF_stdcut"))->Fill(track.tpcNSigmaPr(), track.tofNSigmaPr());
              registry.get<TH2>(HIST("NsigmaTPCvsTOF^2_stdcut"))->Fill(track.tpcNSigmaPr()*track.tpcNSigmaPr(), track.tofNSigmaPr()*track.tofNSigmaPr());
              registry.get<TH2>(HIST("NsigmaTOF_K_stdcut"))->Fill(track.pt(), track.tofNSigmaKa());
              registry.get<TH2>(HIST("NsigmaTPC_K_stdcut"))->Fill(track.pt(), track.tpcNSigmaKa());
              registry.get<TH1>(HIST("DCAxy_stdcut"))->Fill(track.dcaXY());
              registry.get<TH1>(HIST("DCAz_stdcut"))->Fill(track.dcaZ());
              registry.get<TH2>(HIST("DCAxyvspt_stdcut"))->Fill(track.pt(), track.dcaXY());
              registry.get<TH2>(HIST("DCAzvspt_stdcut"))->Fill(track.pt(), track.dcaZ());
              registry.fill(HIST("hcrossedrows_stdcut"), track.tpcNClsCrossedRows());
              registry.fill(HIST("charge_pt_stdcut"), track.pt(), track.sign());

              registry.fill(HIST("charge_stdcut"),track.sign());
              registry.get<TH2>(HIST("NsigmaTOF_d_stdcut"))->Fill(track.pt(), track.tofNSigmaDe());
              registry.get<TH2>(HIST("NsigmaTPC_d_stdcut"))->Fill(track.pt(), track.tpcNSigmaDe());
              registry.fill(HIST("TOFSignalpr_stdcut"), track.pt(), track.tofExpSignalPr(track.tofSignal()));
              registry.fill(HIST("TOFSignalK_stdcut"), track.pt(), track.tofExpSignalKa(track.tofSignal()));
              registry.fill(HIST("TOFSignald_stdcut"), track.pt(), track.tofExpSignalDe(track.tofSignal()));
              
     

              if (track.hasTPC()) {
              registry.fill(HIST("TPC_yess_stdcut"), track.pt());
              registry.fill(HIST("TPC_yess_eta_stdcut"), track.pt(),track.eta());
              if (track.hasTOF()){registry.fill(HIST("TPCandTOF_stdcut"), track.pt());
                                  registry.fill(HIST("TPCandTOF_eta_stdcut"), track.pt(),track.eta()); 
                                  }
              }
            }




    }

    
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EtaPhiHistograms>(cfgc)
  };
}
