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

// jet finder d0 mcp charged task
//
// Authors: Nima Zardoshti

#include "PWGJE/TableProducer/jetfinderhf.cxx"

using JetFinderD0MCParticleLevelCharged = JetFinderHFTask<CandidatesD0Data, CandidatesD0MCD, ParticlesD0, o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents, o2::aod::D0ChargedMCParticleLevelJetConstituentsSub>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderD0MCParticleLevelCharged>(cfgc,
                                                                          SetDefaultProcesses{},
                                                                          TaskName{"jet-finder-d0-mcp-charged"}));

  return WorkflowSpec{tasks};
}
