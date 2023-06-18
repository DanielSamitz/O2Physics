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

/// \file candidateSelectorLcToK0sP.cxx
/// \brief Lc --> K0s+p selection task.
/// \note based on candidateSelectorD0.cxx
///
/// \author Chiara Zampolli <Chiara.Zampolli@cern.ch>, CERN
///         Daniel Samitz, <daniel.samitz@cern.ch>, Vienna

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Tools/ML/model.h"
//#include "Common/Core/TrackSelectorPID.h"
//#include "PWGHF/Utils/utilsDebugLcToK0sP.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc;
using namespace o2::analysis::hf_cuts_lc_to_k0s_p;
using namespace o2::ml;

// candidate rejection types
enum CandidateRejection {
  ptCand = 0,
  mK0s,
  mLambda,
  mGamma,
  ptBach,
  ptV0Dau,
  ptV0,
  d0Bach,
  d0V0,
  tpc,
  tof,
  v0Radius,
  v0DauDCAToPV,
  v0DauDCA,
  decLengthMin,
  decLengthMax,
  v0CosPA,
  fiducialY,
  BDT,
  NCandidateRejection
};


struct HfCandidateSelectorLcToK0sP {
  Produces<aod::HfSelLcToK0sP> hfSelLcToK0sPCandidate;
  Produces<aod::HfCandCascFull2> reducedTree;


  Configurable<bool> applyML{"applyML", false, "Flag to enable or disable ML application"};
  Configurable<std::string> onnxFile{"onnxFile", "", "ONNX file for ML model for Lc+ candidates"};
  Configurable<double> thresholdBDTScore{"thresholdBDTScore", 0., "Threshold value for BDT output scores of Lc+ candidates"};

  Configurable<bool> doMc{"processMc", false, "fill histograms using MC information"};

  Configurable<bool> writeReducedTree{"writeReducedTree", false, "write the selected candidates to a new tree"};

  Configurable<bool> applyYCut{"applyYCut", true, "apply fiducial cut on candidate rapidity"};


  Configurable<double> pPidThreshold{"pPidThreshold", 1.0, "Threshold to switch between low and high p TrackSelectors"};
  Configurable<double> nSigmaTpcMaxLowP{"nSigmaTpcMaxLowP", 2.0, "Max nSigam in TPC for bachelor at low p"};
  Configurable<double> nSigmaTpcMaxHighP{"nSigmaTpcMaxHighP", 9999., "Max nSigam in TPC for bachelor at high p"};
  Configurable<double> nSigmaTofMaxLowP{"nSigmaTofMaxLowP", 9999., "Max nSigam in TOF for bachelor at low p"};
  Configurable<double> nSigmaTofMaxHighP{"nSigmaTofMaxHighP", 3.0, "Max nSigam in TOF for bachelor at high p"};
  Configurable<bool> requireTofLowP{"requireTofLowP", false, "require TOF information for bachelor at low p"};
  Configurable<bool> requireTofHighP{"requireTofHighP", true, "require TOF information for bachelor at high p"};

  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_k0s_p::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Lc candidate selection per pT bin"};

  Configurable<double> v0RadiusMin{"v0RadiusMin", 0.9, "min. V0 radius"};
  Configurable<double> v0CosPAMin{"v0CosPAMin", 0.995, "min. cos pointing angle"};
  Configurable<double> v0DauDCAMax{"v0DauDCAMax", 1.0, "max. DCA between V0 daughters"};
  Configurable<double> v0DauDCAToPVMin{"v0DauDCAToPVMin", 0.1, "min. DCA of V0 dauhgters to PV"};


  HistogramRegistry registry{"registry", {}};

  int dataTypeML;
  OnnxModel model;

  void init(InitContext const&){
    const int nBinsCandidates = 2 + CandidateRejection::NCandidateRejection;
    std::string labels[nBinsCandidates];
    labels[0] = "processed";
    labels[1] = "selected";
    labels[2 + CandidateRejection::ptCand] = "rej. ptCand";
    labels[2 + CandidateRejection::mK0s] = "rej. mK0s";
    labels[2 + CandidateRejection::mLambda] = "rej. mLambda";
    labels[2 + CandidateRejection::mGamma] = "rej. mGamma";
    labels[2 + CandidateRejection::ptBach] = "rej. ptBach";
    labels[2 + CandidateRejection::ptV0Dau] = "rej. ptV0Dau";
    labels[2 + CandidateRejection::ptV0] = "rej. ptV0";
    labels[2 + CandidateRejection::d0Bach] = "rej. d0Bach";
    labels[2 + CandidateRejection::d0V0] = "rej. d0V0";
    labels[2 + CandidateRejection::tpc] = "rej. TPC";
    labels[2 + CandidateRejection::tof] = "rej. TOF";
    labels[2 + CandidateRejection::v0CosPA] = "rej. v0CosPA";
    labels[2 + CandidateRejection::v0DauDCA] = "rej. v0DauDCA";
    labels[2 + CandidateRejection::v0DauDCAToPV] = "rej. v0DauDCAToPV";
    labels[2 + CandidateRejection::v0Radius] = "rej. v0Radius";
    labels[2 + CandidateRejection::decLengthMin] = "rej. decLengthMin";
    labels[2 + CandidateRejection::decLengthMax] = "rej. decLengthMax";
    labels[2 + CandidateRejection::fiducialY] = "rej. fiducialY";
    labels[2 + CandidateRejection::fiducialY] = "rej. BDT";
    AxisSpec axisCandidates = {nBinsCandidates, 0, nBinsCandidates, ""};
    AxisSpec axisBinsPt = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    registry.add("hCandidates", "Candidates;;entries", HistType::kTH1I, {axisCandidates});
    registry.add("hCandidatesVsPtCand", "Candidates;;entries", HistType::kTH2I, {axisCandidates,axisBinsPt});
    for (int iBin = 0; iBin < nBinsCandidates; iBin++) {
      registry.get<TH1>(HIST("hCandidates"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      registry.get<TH2>(HIST("hCandidatesVsPtCand"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    if (doMc){
      registry.add("hCandidatesRecSig", "CandidatesSig;;entries", HistType::kTH1I, {axisCandidates});
      registry.add("hCandidatesVsPtCandRecSig", "Candidates;;entries", HistType::kTH2I, {axisCandidates,axisBinsPt});
      registry.add("hCandidatesRecBg", "CandidatesBg;;entries", HistType::kTH1I, {axisCandidates});
      registry.add("hCandidatesVsPtCandRecBg", "CandidatesBg;;entries", HistType::kTH2I, {axisCandidates,axisBinsPt});
      for (int iBin = 0; iBin < nBinsCandidates; iBin++) {
        registry.get<TH1>(HIST("hCandidatesRecSig"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
        registry.get<TH2>(HIST("hCandidatesVsPtCandRecSig"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
        registry.get<TH1>(HIST("hCandidatesRecBg"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
        registry.get<TH2>(HIST("hCandidatesVsPtCandRecBg"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }

    if (applyML) {
      registry.add<TH1>("hBDTScore", "BDT score distribution for Lc;BDT score;counts", HistType::kTH1F, {{100,0,1,"score"}});
      if (doMc){
        registry.add<TH1>("hBDTScoreRecSig", "BDT score distribution for Lc;BDT score;counts", HistType::kTH1F, {{100,0,1,"score"}});
        registry.add<TH1>("hBDTScoreRecBg", "BDT score distribution for Lc;BDT score;counts", HistType::kTH1F, {{100,0,1,"score"}});
      }
      model.initModel(onnxFile.value, false, 1);
      auto session = model.getSession();
      auto inputShapes = session->GetInputShapes();
      if (inputShapes[0][0] < 0) {
        LOGF(warning, "Model for Lc with negative input shape likely because converted with hummingbird, setting it to 1.");
        inputShapes[0][0] = 1;
      }
      std::vector<float> dummyInput(model.getNumInputNodes(), 1.);
      model.evalModel(dummyInput); // Init the model evaluations
      dataTypeML = session->GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetElementType();
    }

  }


  template <typename T>
  void selection(const T& hfCandCascade, uint32_t &status, float &bdtScore)
  {
    auto candPt = hfCandCascade.pt();
    int ptBin = findBin(binsPt, candPt);
    if (ptBin == -1) {
       SETBIT(status, CandidateRejection::ptCand);
       return;
    }


    if (std::abs(hfCandCascade.v0MK0Short() - RecoDecay::getMassPDG(kK0Short)) > cuts->get(ptBin, "mK0s")) {
      SETBIT(status, CandidateRejection::mK0s);
    }

    if ((std::abs(hfCandCascade.v0MLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda")) || (std::abs(hfCandCascade.v0MAntiLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda"))) {
      SETBIT(status, CandidateRejection::mLambda);
    }

    if (std::abs(hfCandCascade.v0MGamma() - RecoDecay::getMassPDG(kGamma)) < cuts->get(ptBin, "mGamma")) {
      SETBIT(status, CandidateRejection::mGamma);
    }

    if (hfCandCascade.ptProng0() < cuts->get(ptBin, "ptBach")) {
      SETBIT(status, CandidateRejection::ptBach);
    }

    if ((hfCandCascade.ptV0Pos() < cuts->get(ptBin, "ptV0Dau")) || (hfCandCascade.ptV0Neg() < cuts->get(ptBin, "ptV0Dau"))) {
      SETBIT(status, CandidateRejection::ptV0Dau);
    }

    if (hfCandCascade.ptProng1() < cuts->get(ptBin, "ptV0")) {
      SETBIT(status, CandidateRejection::ptV0);
    }

    if (std::abs(hfCandCascade.impactParameter0()) > cuts->get(ptBin, "d0Bach")) {
      SETBIT(status, CandidateRejection::d0Bach);
    }

    if (std::abs(hfCandCascade.impactParameter1()) > cuts->get(ptBin, "d0V0")) {
      SETBIT(status, CandidateRejection::d0V0);
    }

    if (hfCandCascade.v0Radius() < v0RadiusMin) {
      SETBIT(status, CandidateRejection::v0Radius);
    }

    if (hfCandCascade.v0CosPA() < v0CosPAMin) {
      SETBIT(status, CandidateRejection::v0CosPA);
    }

    if (hfCandCascade.dcaV0daughters() > v0DauDCAMax) {
      SETBIT(status, CandidateRejection::v0DauDCA);
    }

    if ((std::abs(hfCandCascade.dcapostopv()) < v0DauDCAToPVMin) || (std::abs(hfCandCascade.dcanegtopv()) < v0DauDCAToPVMin)) {
      SETBIT(status, CandidateRejection::v0DauDCAToPV);
    }

    double nSigmaTpcMax = nSigmaTpcMaxLowP;
    if (hfCandCascade.pProng0()>=pPidThreshold){
      nSigmaTpcMax = nSigmaTpcMaxHighP;
    }
    if (std::abs(hfCandCascade.nSigmaTPCPr0())>nSigmaTpcMax){
      SETBIT(status, CandidateRejection::tpc);
    }

    bool requireTof = requireTofLowP;;
    double nSigmaTofMax = nSigmaTofMaxLowP;
    bool hasTOF = false;
    if (hfCandCascade.nSigmaTOFPr0()>-998.){
      hasTOF = true;
    }
    if (hfCandCascade.pProng0()>=pPidThreshold){
      requireTof = requireTofHighP;
      nSigmaTofMax = nSigmaTofMaxHighP;
    }
    bool tof = true;
    if ((!hasTOF)){
      if (requireTof){
        tof = false;
      }
    } else {
      if (std::abs(hfCandCascade.nSigmaTOFPr0())>nSigmaTofMax){
        tof = false;
      }
    }
    if (!tof){
      SETBIT(status, CandidateRejection::tof);
    }

    if (hfCandCascade.decayLength() < cuts->get(ptBin, "decLengthMin")) {
      SETBIT(status, CandidateRejection::decLengthMin);
    }
    if (hfCandCascade.decayLength() > cuts->get(ptBin, "decLengthMax")) {
      SETBIT(status, CandidateRejection::decLengthMax);
    }

    if (applyYCut){
      double ycut = 0.8;
      if (candPt < 5.){
        ycut = -(0.2/15.)*pow(candPt,2)+(1.9/15.)*candPt+0.5;
      }
      if (std::abs(hfCandCascade.y()) > ycut) {
        SETBIT(status, CandidateRejection::fiducialY);
      }
    }

    if (applyML && status==0){
      std::vector<float> inputFeaturesF{hfCandCascade.ptProng0(), hfCandCascade.ptProng1(), hfCandCascade.impactParameter0(), hfCandCascade.impactParameter1(), hfCandCascade.v0Radius(), hfCandCascade.v0CosPA(), hfCandCascade.v0MLambda(), hfCandCascade.v0MAntiLambda(), hfCandCascade.v0MK0Short(), hfCandCascade.v0MGamma(), hfCandCascade.dcaV0daughters(), hfCandCascade.ptV0Pos(), hfCandCascade.dcapostopv(), hfCandCascade.ptV0Neg(), hfCandCascade.dcanegtopv(), hfCandCascade.nSigmaTPCPr0(), hfCandCascade.nSigmaTOFPr0()};
      std::vector<double> inputFeaturesD{hfCandCascade.ptProng0(), hfCandCascade.ptProng1(), hfCandCascade.impactParameter0(), hfCandCascade.impactParameter1(), hfCandCascade.v0Radius(), hfCandCascade.v0CosPA(), hfCandCascade.v0MLambda(), hfCandCascade.v0MAntiLambda(), hfCandCascade.v0MK0Short(), hfCandCascade.v0MGamma(), hfCandCascade.dcaV0daughters(), hfCandCascade.ptV0Pos(), hfCandCascade.dcapostopv(), hfCandCascade.ptV0Neg(), hfCandCascade.dcanegtopv(), hfCandCascade.nSigmaTPCPr0(), hfCandCascade.nSigmaTOFPr0()};

      if (dataTypeML == 1) {
        auto scoresRaw = model.evalModel(inputFeaturesF);
        bdtScore = scoresRaw[1];
      } else if (dataTypeML == 11) {
        auto scoresRaw = model.evalModel(inputFeaturesD);
        bdtScore = scoresRaw[1];
      } else {
        LOG(error) << "Error running model inference for Lc: Unexpected input data type.";
      }
      if (bdtScore < thresholdBDTScore){
        SETBIT(status, CandidateRejection::BDT);
      }
    }


  }

  
  void process(aod::HfCandCascFull const& candidates)
  {
    uint32_t statusLc = 0; // final selection flag

    for (const auto& candidate : candidates) {               // looping over cascade candidates


      statusLc = 0;
      float bdtScore = -100.;
      selection(candidate,statusLc, bdtScore);
      registry.fill(HIST("hBDTScore"), bdtScore);
      if (std::abs(candidate.flagMc()) == 1){
        registry.fill(HIST("hBDTScoreRecSig"), bdtScore);
      }
      else {
        registry.fill(HIST("hBDTScoreRecBg"), bdtScore);
      }
      //SETBIT(statusLc,0);
      hfSelLcToK0sPCandidate(statusLc);

      double pt = candidate.pt();

      registry.fill(HIST("hCandidates"), 0.5);
      registry.fill(HIST("hCandidatesVsPtCand"), 0.5, pt);
      if (statusLc == 0){
        registry.fill(HIST("hCandidates"), 1.5);
        registry.fill(HIST("hCandidatesVsPtCand"), 1.5, pt);
      }
      int bin=2;
      for (uint32_t i = 1; i<TMath::Power(2,CandidateRejection::NCandidateRejection)-1; i *= 2){
        if (i&statusLc){
          registry.fill(HIST("hCandidates"),bin+0.5);
          registry.fill(HIST("hCandidatesVsPtCand"),bin+0.5, pt);
        }
        bin++;
      }


      if (doMc){
        if (std::abs(candidate.flagMc()) == 1){
          registry.fill(HIST("hCandidatesRecSig"), 0.5);
          registry.fill(HIST("hCandidatesVsPtCandRecSig"), 0.5, pt);
          if (statusLc == 0){
            registry.fill(HIST("hCandidatesRecSig"), 1.5);
            registry.fill(HIST("hCandidatesVsPtCandRecSig"), 1.5, pt);
          }
          int bin=2;
          for (uint32_t i = 1; i<TMath::Power(2,CandidateRejection::NCandidateRejection)-1; i *= 2){
            if (i&statusLc){
              registry.fill(HIST("hCandidatesRecSig"),bin+0.5);
              registry.fill(HIST("hCandidatesVsPtCandRecSig"),bin+0.5, pt);
            }
            bin++;
          }
        }
        else {
          registry.fill(HIST("hCandidatesRecBg"), 0.5);
          registry.fill(HIST("hCandidatesVsPtCandRecBg"), 0.5, pt);
          if (statusLc == 0){
            registry.fill(HIST("hCandidatesRecBg"), 1.5);
            registry.fill(HIST("hCandidatesVsPtCandRecBg"), 1.5, pt);
          }
          int bin=2;
          for (uint32_t i = 1; i<TMath::Power(2,CandidateRejection::NCandidateRejection)-1; i *= 2){
            if (i&statusLc){
              registry.fill(HIST("hCandidatesRecBg"),bin+0.5);
              registry.fill(HIST("hCandidatesVsPtCandRecBg"),bin+0.5, pt);
            }
            bin++;
          }
        }
      }
      if (statusLc==0 && writeReducedTree){
        reducedTree(
          candidate.bcId(),
          candidate.numContrib(),
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
          candidate.pProng0(),
          candidate.impactParameterNormalised1(),
          candidate.ptProng1(),
          candidate.pProng1(),
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
          candidate.v0Radius(),
          candidate.v0CosPA(),
          candidate.v0MLambda(),
          candidate.v0MAntiLambda(),
          candidate.v0MK0Short(),
          candidate.v0MGamma(),
          candidate.v0CtK0Short(),
          candidate.v0CtLambda(),
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
          candidate.nSigmaTPCPr0(),
          candidate.nSigmaTOFPr0(),
          candidate.m(),
          candidate.pt(),
          candidate.p(),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.ct(),
          candidate.eta(),
          candidate.phi(),
          candidate.y(),
          candidate.e(),
          candidate.flagMc(),
          candidate.originMcRec(),
          false // dummy
          ); 
      }
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfcg)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorLcToK0sP>(cfcg));
  return workflow;
}
