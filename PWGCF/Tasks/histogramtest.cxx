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
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <TParameter.h>
#include <TH1F.h>
#include <iostream>

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
#include "TDatabasePDG.h"

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

struct QAHistograms {
  // using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> sign{"sign", 1, ""};
  Configurable<float> min_P{"min_P", 0.0, ""};
  Configurable<float> max_P{"max_P", 100.0, ""};
  Configurable<float> eta{"eta", 10.0, ""};
  Configurable<int> tpcNClsFound{"tpcNClsFound", 0, ""};
  Configurable<float> tpcChi2NCl{"tpcChi2NCl", 100.0, ""};
  Configurable<float> crossedRowsfindableCls{"crossedRows/findableCls", 1.0, ""};
  Configurable<float> dcaXY{"dcaXY", 10.0, ""};
  Configurable<float> dcaZ{"dcaZ", 10.0, ""};
  Configurable<int> itsNCls{"itsNCls", -1, ""};
  Configurable<float> PIDtrshld{"PIDtrshld", 0.75, ""};
  Configurable<int> particlePDG{"particlePDG", 2212, ""};
  Configurable<float> tpcNSigma{"tpcNSigma", 999.0, ""};
  Configurable<float> tpctofNSigma{"tpctofNSigma", 999.0, ""};

  std::map<std::string, float> trackcuts;
  std::map<std::string, double> PIDcuts;

  void init(o2::framework::InitContext&)
  {
    registry.add("TPCSignal_nocuts", "TPC signal without cuts", kTH2F, {{{200, 0., 5.0, "#it{p}_{inner} (GeV/#it{c})"}, {1000, 0., 1000.0, "dE/dx in TPC (arbitrary units)"}}});
    registry.add("TOFSignal_nocuts", "TOF signal without cuts", kTH2F, {{{200, 0., 5.0, "#it{p} (GeV/#it{c})"}, {100, 0., 1.5, "#beta"}}});

    trackcuts = {
      {"min_P", min_P}, {"max_P", max_P}, {"eta", eta}, {"tpcNClsFound", tpcNClsFound}, {"crossedRows/findableCls", crossedRowsfindableCls}, {"tpcChi2NCl", tpcChi2NCl}, {"dcaXY", dcaXY}, {"dcaZ", dcaXY}, {"itsNCls", itsNCls}};

    PIDcuts = {
      {"sign", sign}, {"PIDtrshld", PIDtrshld}, {"particlePDG", particlePDG}, {"tpcNSigma", tpcNSigma}, {"tpctofNSigma", tpctofNSigma}};

    registry.add("eta", "Eta; #eta; counts", kTH1F, {{200, -2.5, 2.5}});
    registry.add("phi", "Phi; #phi; counts", kTH1F, {{200, 0., 2. * M_PI}});
    registry.add("etaphi", "eta vs phi", kTH2F, {{100, -1.5, 1.5, "#eta"}, {200, 0., 2. * M_PI, "#phi"}});
    registry.add("px", "px; p_{x} (GeV/#it{c}); counts ", kTH1F, {{100, 0., 5.}});
    registry.add("py", "py; p_{y} (GeV/#it{c}); counts ", kTH1F, {{100, 0., 5.}});
    registry.add("pz", "pz; p_{z} (GeV/#it{c}); counts", kTH1F, {{100, 0., 5.}});
    registry.add("p", "Momentum; p (GeV/#it{c}); counts", kTH1F, {{100, 0., 5.}});
    registry.add("pt", "Transverse momentum; #it{p}_{T} (GeV/#it{c}); counts", kTH1F, {{100, 0., 5.}});
    registry.add("sign", "Sign; sign: counts", kTH1F, {{3, -1.5, 1.5}});
    registry.add("TPCSignal", "TPC Signal", kTH2F, {{{200, 0., 5.0, "#it{p}_{inner} (GeV/#it{c})"}, {1000, 0., 1000.0, "dE/dx in TPC (arbitrary units)"}}});
    registry.add("TOFSignal", "TOF Signal", kTH2F, {{200, 0., 5.0, "#it{p} (GeV/#it{c})"}, {100, 0., 1.5, "#beta"}});
    registry.add("dcaxy_to_p", "DCAxy vs p", kTH2F, {{100, 0., 5.0, "p (GeV/#it{c})"}, {500, -1.5, 1.5, "DCAxy (cm)"}});
    registry.add("dcaxy_to_pt", "DCAxy vs #it{p}_{T} ", kTH2F, {{100, 0., 5., "#it{p}_{T} (GeV/#it{c})"}, {500, -1.5, 1.5, "DCAxy (cm)"}});
    registry.add("dcaz_to_p", "DCAz vs p", kTH2F, {{100, 0., 5., "p (GeV/#it{c})"}, {500, -1.5, 1.5, "DCAz (cm)"}});
    registry.add("dcaz_to_pt", "DCAz vs #it{p}_{T}", kTH2F, {{100, 0., 5., "#it{p}_{T} (GeV/#it{c})"}, {500, -1.5, 1.5, "DCAz (cm)"}});
    registry.add("dcaxy", "DCAxy; DCAxy (cm); counts", kTH1F, {{100, -1., 1.}});
    registry.add("dcaz", "DCAz; DCAz (cm); counts", kTH1F, {{100, -1., 1.}});
    registry.add("crossed_rows", "Number of crossed rows; Number of crossed rows; counts", kTH1F, {{160, -0.5, 159.5}});
    registry.add("TPCClusters", "Number of TPC Clusters; Number of TPC Clusters; counts", kTH1F, {{163, -0.5, 162.5}});
    registry.add("ITSClusters", "Number of ITS Clusters; Number of ITS Clusters; counts", kTH1F, {{10, -0.5, 9.5}});
    registry.add("ITSchi2", "ITS chi2; ITS #chi^{2}; counts", kTH1F, {{100, 0.0, 40.}});
    registry.add("TPCchi2", "TPC chi2; TPC #chi^{2} counts;", kTH1F, {{100, 0.0, 6.}});
    //registry.add("collisionId", "Collisions; collision index; counts", kTH1F, {{1000, 0., 12000}});
    if (PIDcuts["particlePDG"] == 2212) {
      registry.add("nsigmaTOFPr", "Proton n_{#sigma} TOF", kTH2F, {{100, 0., 5., "#it{p} (GeV/#it{c})"}, {100, -10., 10., "n^{TOF}_{#sigma, p}"}});
      registry.add("nsigmaTPCPr", "Proton n_{#sigma} TPC", kTH2F, {{100, 0., 5., "#it{p} (GeV/#it{c})"}, {100, -10., 10., "n^{TPC}_{#sigma, p}"}});
      registry.add("nsigmaTPCTOF", "Proton n_{#sigma} TPC vs n_{#sigma} TOF", kTH3F, {{100, 0., 5., "#it{p} (GeV/#it{c})"}, {100, -10., 10., "n^{TPC}_{#sigma, p}"}, {100, -10., 10., "n^{TOF}_{#sigma, p}"}});
    }
    if (PIDcuts["particlePDG"] == 1000010020) {
      registry.add("nsigmaTOFDe", "Deuteron n_{#sigma} TOF ", kTH2F, {{100, 0., 5., "#it{p} (GeV/#it{c})"}, {100, -10., 10., "n^{TOF}_{#sigma, d}"}});
      registry.add("nsigmaTPCDe", "Deuteron n_{#sigma} TPC", kTH2F, {{100, 0., 5., "#it{p} (GeV/#it{c})"}, {100, -10., 10., "n^{TPC}_{#sigma, d}"}});
      registry.add("nsigmaTPCTOF", "Deuteron n_{#sigma} TPC vs n_{#sigma} TOF", kTH3F, {{100, 0., 5., "#it{p} (GeV/#it{c})"}, {100, -10., 10., "n^{TPC}_{#sigma, d}"}, {100, -10., 10., "n^{TOF}_{#sigma, d}"}});
    }
    // if(singletrackselector::HasTOF)
    registry.add("events", "vertex position along z; v_z (cm); counts", kTH1F, {{20, -20., 20.}});

  }

  void process(aod::SingleCollSels const& collisions, aod::SingleTrackSel const& tracks)
  {

    for (auto& track : tracks) {
      registry.fill(HIST("TPCSignal_nocuts"), track.tpcInnerParam(), track.tpcSignal());
      registry.fill(HIST("TOFSignal_nocuts"), track.p(), track.beta());
      if (!track.trackCuts(&trackcuts) || !track.pidCuts(&PIDcuts))
        continue;
      else {

        registry.fill(HIST("eta"), track.eta());
        registry.fill(HIST("phi"), track.phi());
        registry.fill(HIST("etaphi"), track.eta(), track.phi());
        registry.fill(HIST("px"), track.px());
        registry.fill(HIST("py"), track.py());
        registry.fill(HIST("pz"), track.pz());
        registry.fill(HIST("p"), track.p());
        registry.fill(HIST("pt"), track.pt());
        registry.fill(HIST("sign"), track.sign());
        registry.fill(HIST("TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("TOFSignal"), track.p(), track.beta());
        registry.fill(HIST("dcaxy_to_p"), track.p(), track.dcaXY());
        registry.fill(HIST("dcaxy_to_pt"), track.pt(), track.dcaXY());
        registry.fill(HIST("dcaz_to_p"), track.p(), track.dcaZ());
        registry.fill(HIST("dcaz_to_pt"), track.pt(), track.dcaZ());
        registry.fill(HIST("dcaxy"), track.dcaXY());
        registry.fill(HIST("dcaz"), track.dcaZ());
        registry.fill(HIST("TPCClusters"), track.tpcNClsFound());
        registry.fill(HIST("ITSClusters"), track.itsNCls());
        registry.fill(HIST("ITSchi2"), track.itsChi2NCl());
        registry.fill(HIST("TPCchi2"), track.tpcChi2NCl());
        registry.fill(HIST("crossed_rows"), track.tpcNClsCrossedRows());
        // registry.fill(HIST("collisionId"), track.collisionId());
        if (PIDcuts["particlePDG"] == 2212) {
          registry.fill(HIST("nsigmaTOFPr"), track.p(), track.tofNSigmaPr());
          registry.fill(HIST("nsigmaTPCPr"), track.p(), track.tpcNSigmaPr());
          registry.fill(HIST("nsigmaTPCTOF"), track.p(), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
        if (PIDcuts["particlePDG"] == 1000010020) {
          registry.fill(HIST("nsigmaTOFDe"), track.p(), track.tofNSigmaDe());
          registry.fill(HIST("nsigmaTPCDe"), track.p(), track.tpcNSigmaDe());
          registry.fill(HIST("nsigmaTPCTOF"), track.p(), track.tpcNSigmaDe(), track.tofNSigmaDe());
        }
      }
    }

    for (auto& collision : collisions) {
      registry.fill(HIST("events"), collision.posZ());
    }
  }
  // PROCESS_SWITCH(EtaPhiHistograms, processSelected, "process filtered track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QAHistograms>(cfgc)};
}
