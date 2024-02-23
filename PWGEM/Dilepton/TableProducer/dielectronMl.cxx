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
/// \file dielectronMl.cxx
/// \task for testing ML application for dielectron analyses
/// \author Daniel Samitz, <daniel.samitz@cern.ch>, SMI Vienna
///         Elisa Meninno, <elisa.meninno@cern.ch>, SMI Vienna

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/DataModel/SelectionTables.h"
#include "PWGEM/Dilepton/Utils/MlResponseDielectronSingleTrack.h"
#include "PWGEM/Dilepton/Utils/MlResponseDielectronPair.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::aod;

using MySkimmedTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov>;
using MySkimmedTracksWithPID = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelCov>;
using MyTracksWithPID = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                  aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                  aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;

// define some default values for single track analysis
namespace dielectron_ml_cuts_single_track
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};

static constexpr int nBinsPt = 1;
static constexpr int nCutScores = 2;
// default values for the pT bin edges, offset by 1 from the bin numbers in cuts array
constexpr double binsPt[nBinsPt + 1] = {
  0.,
  999.};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

// default values for the ML model paths, one model per pT bin
static const std::vector<std::string> modelPaths = {
  ""};

// default values for the cut directions
constexpr int cutDir[nCutScores] = {CutSmaller, CutGreater};
auto vecCutDir = std::vector<int>{cutDir, cutDir + nCutScores};

// default values for the cuts
constexpr double cuts[nBinsPt][nCutScores] = {
  {0.5, 0.5}};

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0"};

// column labels
static const std::vector<std::string> labelsCutScore = {"Signal", "Background"};
} // namespace dielectron_ml_cuts_single_track

// define some default values for pair analysis
namespace dielectron_ml_cuts_pair
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};

static constexpr int nBinsM = 1;
static constexpr int nCutScores = 2;
// default values for the mass bin edges, offset by 1 from the bin numbers in cuts array
constexpr double binsM[nBinsM + 1] = {
  0.,
  999.};
auto vecBinsM = std::vector<double>{binsM, binsM + nBinsM + 1};

// default values for the ML model paths, one model per mass bin
static const std::vector<std::string> modelPaths = {
  ""};

// default values for the cut directions
constexpr int cutDir[nCutScores] = {CutSmaller, CutGreater};
auto vecCutDir = std::vector<int>{cutDir, cutDir + nCutScores};

// default values for the cuts
constexpr double cuts[nBinsM][nCutScores] = {
  {0.5, 0.5}};

// row labels
static const std::vector<std::string> labelsM = {
  "mass bin 0"};

// column labels
static const std::vector<std::string> labelsCutScore = {"Signal", "Background"};
} // namespace dielectron_ml_cuts_pair

struct DielectronMlSingleTrack {
  Produces<aod::dqMlSelTrack> singleTrackSelection;

  // ML inference
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{dielectron_ml_cuts_single_track::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{dielectron_ml_cuts_single_track::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {dielectron_ml_cuts_single_track::cuts[0], dielectron_ml_cuts_single_track::nBinsPt, dielectron_ml_cuts_single_track::nCutScores, dielectron_ml_cuts_single_track::labelsPt, dielectron_ml_cuts_single_track::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(dielectron_ml_cuts_single_track::nCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{""}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{""}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  // quality cuts (from treeCreatorElectronMl.cxx)
  Configurable<bool> doQualityCuts{"doQualityCuts", false, "enable track quality cuts"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};

  o2::analysis::MlResponseDielectronSingleTrack<float> mlResponse;
  o2::ccdb::CcdbApi ccdbApi;
  std::vector<std::shared_ptr<TH1>> hModelScore;
  std::vector<std::shared_ptr<TH2>> hModelScoreVsPt;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (doprocessSkimmedSingleTrack || doprocessAO2DSingleTrack) {
      mlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        mlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        mlResponse.setModelPathsLocal(onnxFileNames);
      }
      mlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      mlResponse.init();

      // load histograms for ML score
      AxisSpec axisScore = {100, 0, 1, "score"};
      AxisSpec axisBinsPt = {binsPtMl, "#it{p}_{T} (GeV/#it{c})"};
      for (int classMl = 0; classMl < nClassesMl; classMl++) {
        hModelScore.push_back(registry.add<TH1>("hMlScore" + TString(cutsMl->getLabelsCols()[classMl]), "Model score distribution;Model score;counts", HistType::kTH1F, {axisScore}));
        hModelScoreVsPt.push_back(registry.add<TH2>("hMlScore" + TString(cutsMl->getLabelsCols()[classMl]) + "VsPt", "Model score distribution;Model score;counts", HistType::kTH2F, {axisScore, axisBinsPt}));
      }

      registry.add("hMlSelection", "ML selection status;counts", {HistType::kTH1F, {{2, -0.5, 1.5}}});
      registry.add("hMlSelectionVsPt", "ML selection status;counts", {HistType::kTH2F, {{2, -0.5, 1.5}, axisBinsPt}});
    }
  }

  template <typename T>
  bool applyQualityCuts(T const& track)
  {
    // consistent with treeCreatorElectronMl.cxx
    if (!track.hasITS()) {
      return false;
    }
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }
    if (track.itsChi2NCl() < -1) { // if tracks are not reconstructed properly, chi2/ITSncls is set to -999;
      return false;
    }
    if (abs(track.eta()) > maxeta) {
      return false;
    }
    if (track.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }
    if (abs(track.dcaXY()) > 1.) {
      return false;
    }
    if (abs(track.dcaZ() > 1.)) {
      return false;
    }
    return true;
  }

  template <typename T>
  void runSingleTracks(T const& tracks)
  {
    for (const auto& track : tracks) {
      auto pt = track.pt();
      if (doQualityCuts && !applyQualityCuts(track)) {
        singleTrackSelection(false);
        registry.fill(HIST("hMlSelection"), 0);
        registry.fill(HIST("hMlSelectionVsPt"), 0, pt);
        continue;
      }
      std::vector<float> inputFeatures = mlResponse.getInputFeatures(track);
      std::vector<float> outputMl = {};

      bool isSelected = mlResponse.isSelectedMl(inputFeatures, pt, outputMl);
      for (int classMl = 0; classMl < nClassesMl; classMl++) {
        hModelScore[classMl]->Fill(outputMl[classMl]);
        hModelScoreVsPt[classMl]->Fill(outputMl[classMl], pt);
      }
      singleTrackSelection(isSelected);
      registry.fill(HIST("hMlSelection"), isSelected);
      registry.fill(HIST("hMlSelectionVsPt"), isSelected, pt);
    }
  }

  void processSkimmedSingleTrack(MySkimmedTracksWithPID const& tracks)
  {
    runSingleTracks(tracks);
  }
  PROCESS_SWITCH(DielectronMlSingleTrack, processSkimmedSingleTrack, "Apply ML selection on skimmed output on single tracks", true);

  void processAO2DSingleTrack(MyTracksWithPID const& tracks)
  {
    runSingleTracks(tracks);
  }
  PROCESS_SWITCH(DielectronMlSingleTrack, processAO2DSingleTrack, "Apply ML selection on skimmed output on single tracks", false);

  void processDummyAO2D(aod::Collisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(DielectronMlSingleTrack, processDummyAO2D, "Dummy", false);

  void processDummySkimmed(aod::ReducedEvents const&)
  {
    // dummy
  }
  PROCESS_SWITCH(DielectronMlSingleTrack, processDummySkimmed, "Dummy", false);
};

struct DielectronMlPair {
  Produces<aod::dqMlSelDielectron> pairSelection;

  // ML inference
  Configurable<std::vector<double>> binsMMl{"binsMMl", std::vector<double>{dielectron_ml_cuts_pair::vecBinsM}, "Mass bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{dielectron_ml_cuts_pair::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {dielectron_ml_cuts_pair::cuts[0], dielectron_ml_cuts_pair::nBinsM, dielectron_ml_cuts_pair::nCutScores, dielectron_ml_cuts_pair::labelsM, dielectron_ml_cuts_pair::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", static_cast<int>(dielectron_ml_cuts_pair::nCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{""}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{""}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::MlResponseDielectronPair<float> mlResponse;
  o2::ccdb::CcdbApi ccdbApi;
  std::vector<std::shared_ptr<TH1>> hModelScore;
  std::vector<std::shared_ptr<TH2>> hModelScoreVsM;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (doprocessPair) {
      mlResponse.configure(binsMMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        mlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        mlResponse.setModelPathsLocal(onnxFileNames);
      }
      mlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      mlResponse.init();

      // load histograms for ML score
      AxisSpec axisScore = {100, 0, 1, "score"};
      AxisSpec axisBinsM = {binsMMl, "#it{M} (GeV/#it{c^{2}})"};
      for (int classMl = 0; classMl < nClassesMl; classMl++) {
        hModelScore.push_back(registry.add<TH1>("hMlScore" + TString(cutsMl->getLabelsCols()[classMl]), "Model score distribution;Model score;counts", HistType::kTH1F, {axisScore}));
        hModelScoreVsM.push_back(registry.add<TH2>("hMlScore" + TString(cutsMl->getLabelsCols()[classMl]) + "VsM", "Model score distribution;Model score;counts", HistType::kTH2F, {axisScore, axisBinsM}));
      }
      registry.add("hMlSelection", "ML selection status;counts", {HistType::kTH1F, {{2, -0.5, 1.5}}});
      registry.add("hMlSelectionVsM", "ML selection status;counts", {HistType::kTH2F, {{2, -0.5, 1.5}, axisBinsM}});
    }
  }

  template <typename E, typename T>
  void doPair(E const& event, T const& tracks1, T const& tracks2)
  {
    int pos = 0;
    std::vector<int> vecSel = {};
    for (const auto& track1 : tracks1) {
      if ((!event.isEventSelected() == 1) || (!track1.isBarrelSelected() > 0)) {
        pairSelection(vecSel, -1);
        continue;
      }
      for (const auto& track2 : tracks2) {
        if (track2.globalIndex() <= track1.globalIndex()) {
          continue;
        }
        if (!track2.isBarrelSelected() > 0) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(track1.pt(), track1.eta(), track1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(track2.pt(), track2.eta(), track2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        auto m = v12.M();
        std::vector<float> inputFeatures = mlResponse.getInputFeatures(track1, track2);
        std::vector<float> outputMl = {};

        bool isSelected = mlResponse.isSelectedMl(inputFeatures, m, outputMl);
        for (int classMl = 0; classMl < nClassesMl; classMl++) {
          hModelScore[classMl]->Fill(outputMl[classMl]);
          hModelScoreVsM[classMl]->Fill(outputMl[classMl], m);
        }
        registry.fill(HIST("hMlSelection"), isSelected);
        registry.fill(HIST("hMlSelectionVsM"), isSelected, m);
        vecSel.push_back(isSelected);
      }
      pairSelection(vecSel, pos);
      pos++;
      vecSel = {};
    }
  }

  void processPair(MyEventsSelected::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks)
  {
    // dummy value for magentic field. ToDo: take it from ccdb!
    float d_bz = 1.;
    mlResponse.setBz(d_bz);
    doPair(event, tracks, tracks);
  }
  PROCESS_SWITCH(DielectronMlPair, processPair, "Apply ML selection at pair level", false);

  void processDummyAO2D(aod::Collisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(DielectronMlPair, processDummyAO2D, "Dummy", false);

  void processDummySkimmed(aod::ReducedEvents const&)
  {
    // dummy
  }
  PROCESS_SWITCH(DielectronMlPair, processDummySkimmed, "Dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DielectronMlSingleTrack>(cfgc),
    adaptAnalysisTask<DielectronMlPair>(cfgc)};
}
