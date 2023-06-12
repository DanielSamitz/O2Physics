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

/// \file treeCreatorLcToK0sP.cxx
/// \brief Writer of the cascade candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///        Modified version of treeCreatorLcToPKPi.cxx
///
/// \author Daniel Samitz <daniel.samitz@cern.ch>

//#include "Common/Core/trackUtilities.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
//#include "ReconstructionDataFormats/DCA.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

/// Writes the full information in an output TTree
struct HfTreeCreatorLcToK0sP {
  Produces<o2::aod::HfCandCascFull> rowCandidateFull;
  Produces<o2::aod::HfCandCascFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandCascFullParticles> rowCandidateFullParticles;

  Configurable<double> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of candidates to store in the tree"};

  void init(InitContext const&)
  {
  }

  template <typename T, typename U>
  void fillCandidate(const T& candidate, const U& bach, int8_t flagMc, int8_t originMcRec)
  {
    rowCandidateFull(
      bach.collision().bcId(),
      bach.collision().numContrib(),
      candidate.posX(),
      candidate.posY(),
      candidate.posZ(),
      candidate.xSecondaryVertex(),
      candidate.ySecondaryVertex(),
      candidate.zSecondaryVertex(),
      candidate.errorDecayLength(),
      candidate.errorDecayLengthXY(),
      candidate.chi2PCA(),
      candidate.rSecondaryVertex(),
      candidate.decayLength(),
      candidate.decayLengthXY(),
      candidate.decayLengthNormalised(),
      candidate.decayLengthXYNormalised(),
      candidate.impactParameterNormalised0(),
      candidate.ptProng0(),
      RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
      candidate.impactParameterNormalised1(),
      candidate.ptProng1(),
      RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
      candidate.pxProng0(),
      candidate.pyProng0(),
      candidate.pzProng0(),
      candidate.pxProng1(),
      candidate.pyProng1(),
      candidate.pzProng1(),
      candidate.impactParameter0(),
      candidate.impactParameter1(),
      candidate.errorImpactParameter0(),
      candidate.errorImpactParameter1(),
      candidate.v0x(),
      candidate.v0y(),
      candidate.v0z(),
      candidate.v0radius(),
      candidate.v0cosPA(),
      candidate.mLambda(),
      candidate.mAntiLambda(),
      candidate.mK0Short(),
      candidate.mGamma(),
      o2::aod::hf_cand_casc::ctV0K0s(candidate),
      o2::aod::hf_cand_casc::ctV0Lambda(candidate),
      candidate.dcaV0daughters(),
      candidate.pxpos(),
      candidate.pypos(),
      candidate.pzpos(),
      candidate.ptV0Pos(),
      candidate.dcapostopv(),
      candidate.pxneg(),
      candidate.pyneg(),
      candidate.pzneg(),
      candidate.ptV0Neg(),
      candidate.dcanegtopv(),
      bach.tpcNSigmaPr(),
      bach.tofNSigmaPr(),
      o2::aod::hf_cand_casc::invMassLcToK0sP(candidate),
      candidate.pt(),
      candidate.p(),
      candidate.cpa(),
      candidate.cpaXY(),
      o2::aod::hf_cand_3prong::ctLc(candidate),
      candidate.eta(),
      candidate.phi(),
      o2::aod::hf_cand_3prong::yLc(candidate),
      o2::aod::hf_cand_3prong::eLc(candidate),
      flagMc,
      originMcRec);
  }

  template <typename T>
  void fillEvent(const T& collision)
  {
    rowCandidateFullEvents(
      collision.bcId(),
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ());
  }

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const& mccollisions,
                 soa::Join<aod::HfCandCascade, aod::HfCandCascadeMcRec, aod::HfSelLcToK0sP> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& particles,
                 aod::BigTracksPID const& tracks)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto const& collision : collisions) {
      fillEvent(collision);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto const& candidate : candidates) {
      auto bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor
      double pseudoRndm = bach.pt() * 1000. - (int16_t)(bach.pt() * 1000);
      if (candidate.isSelLcToK0sP() >= 1 && pseudoRndm < downSampleBkgFactor) {
        fillCandidate(candidate, bach, candidate.flagMcMatchRec(), candidate.originMcRec());
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto const& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(array{particle.px(), particle.py(), particle.pz()},
                       RecoDecay::getMassPDG(particle.pdgCode())),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorLcToK0sP, processMc, "Process MC tree writer", true);

  void processData(aod::Collisions const& collisions,
                   soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP> const& candidates,
                   aod::BigTracksPID const& tracks)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto const& collision : collisions) {
      fillEvent(collision);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto const& candidate : candidates) {
      auto bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor
      double pseudoRndm = bach.pt() * 1000. - (int16_t)(bach.pt() * 1000);
      if (candidate.isSelLcToK0sP() >= 1 && pseudoRndm < downSampleBkgFactor) {
        fillCandidate(candidate, bach, 0, 0);
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorLcToK0sP, processData, "Process data tree writer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorLcToK0sP>(cfgc),
  };
}
