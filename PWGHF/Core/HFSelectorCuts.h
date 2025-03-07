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

/// \file HFSelectorCuts.h
/// \brief Default pT bins and cut arrays for heavy-flavour selectors and analysis tasks

#ifndef HF_SELECTOR_CUTS_H_
#define HF_SELECTOR_CUTS_H_

#include "Framework/Configurable.h"
#include <vector>
#include <string>

namespace o2::analysis
{
namespace pdg
{
enum Code {
  kD0 = 421,
  kD0bar = -421,
  kDPlus = 411,
  kDs = 431,
  kLambdaCPlus = 4122,
  kXiCPlus = 4232,
  kXiCCPlusPlus = 4422,
  kLambdaB0 = 5122,
  kJpsi = 443,
  kChic1 = 20443,
  kBPlus = 521,
  kX3872 = 9920443
};
} // namespace pdg

/// Finds pT bin in an array.
/// \param bins  array of pT bins
/// \param value  pT
/// \return index of the pT bin
/// \note Accounts for the offset so that pt bin array can be used to also configure a histogram axis.
template <typename T1, typename T2>
int findBin(T1 const& bins, T2 value)
{
  if (value < bins->front()) {
    return -1;
  }
  if (value >= bins->back()) {
    return -1;
  }
  return std::distance(bins->begin(), std::upper_bound(bins->begin(), bins->end(), value)) - 1;
}

// namespace per channel

namespace hf_cuts_single_track
{
static constexpr int npTBinsTrack = 6;
static constexpr int nCutVarsTrack = 2;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any candidate type (2-prong, 3-prong)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBinsTrack[npTBinsTrack + 1] = {
  0,
  0.5,
  1.0,
  1.5,
  2.0,
  3.0,
  1000.0};
auto pTBinsTrack_v = std::vector<double>{pTBinsTrack, pTBinsTrack + npTBinsTrack + 1};

// default values for the cuts
constexpr double cutsTrack[npTBinsTrack][nCutVarsTrack] = {{0.0025, 10.},  /* 0   < pt < 0.5 */
                                                           {0.0025, 10.},  /* 0.5 < pt < 1 */
                                                           {0.0025, 10.},  /* 1   < pt < 1.5 */
                                                           {0.0025, 10.},  /* 1.5 < pt < 2 */
                                                           {0.0000, 10.},  /* 2   < pt < 3 */
                                                           {0.0000, 10.}}; /* 3   < pt < 1000 */

// row labels
static const std::vector<std::string> pTBinLabelsTrack{};

// column labels
static const std::vector<std::string> cutVarLabelsTrack = {"min_dcaxytoprimary", "max_dcaxytoprimary"};
} // namespace hf_cuts_single_track

namespace hf_cuts_bdt_multiclass
{
static constexpr int npTBins = 1;
static constexpr int nCutBDTScores = 3;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any charm candidate
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0.,
  1000.0};
auto pTBinsVec = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cutsBDT[npTBins][nCutBDTScores] = {{0.1, 0.5, 0.5}};

// row labels
static const std::vector<std::string> pTBinLabels{};

// column labels
static const std::vector<std::string> cutBDTLabels = {"BDTbkg", "BDTprompt", "BDTnonprompt"};
} // namespace hf_cuts_bdt_multiclass

namespace hf_cuts_presel_2prong
{
static constexpr int npTBins = 2;
static constexpr int nCutVars = 4;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any 2-prong candidate
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  1.,
  5.,
  1000.0};
auto pTBinsVec = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cuts[npTBins][nCutVars] = {{1.65, 2.15, 0.5, 100.},  /* 1 < pt < 5 */
                                            {1.65, 2.15, 0.5, 100.}}; /* 5 < pt < 1000 */

// row labels
static const std::vector<std::string> pTBinLabels{};

// column labels
static const std::vector<std::string> cutVarLabels = {"massMin", "massMax", "cosp", "d0d0"};
} // namespace hf_cuts_presel_2prong

namespace hf_cuts_presel_3prong
{
static constexpr int npTBins = 2;
static constexpr int nCutVars = 4;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any 3-prong candidate
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  1.,
  5.,
  1000.0};
auto pTBinsVec = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cuts[npTBins][nCutVars] = {{1.75, 2.05, 0.7, 0.02},  /* 1 < pt < 5 */
                                            {1.75, 2.05, 0.5, 0.02}}; /* 5 < pt < 1000 */

// row labels
static const std::vector<std::string> pTBinLabels{};

// column labels
static const std::vector<std::string> cutVarLabels = {"massMin", "massMax", "cosp", "decL"};
} // namespace hf_cuts_presel_3prong

namespace hf_cuts_d0_topik
{
static constexpr int npTBins = 25;
static constexpr int nCutVars = 14;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0,
  0.5,
  1.0,
  1.5,
  2.0,
  2.5,
  3.0,
  3.5,
  4.0,
  4.5,
  5.0,
  5.5,
  6.0,
  6.5,
  7.0,
  7.5,
  8.0,
  9.0,
  10.0,
  12.0,
  16.0,
  20.0,
  24.0,
  36.0,
  50.0,
  100.0};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cuts[npTBins][nCutVars] = {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06},   /* 0   < pT < 0.5 */
                                            {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06},   /* 0.5 < pT < 1   */
                                            {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06},  /* 1   < pT < 1.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06},  /* 1.5 < pT < 2   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06},  /* 2   < pT < 2.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06},  /* 2.5 < pT < 3   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},  /* 3   < pT < 3.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},  /* 3.5 < pT < 4   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 4   < pT < 4.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 4.5 < pT < 5   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 5   < pT < 5.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 5.5 < pT < 6   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 6   < pT < 6.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 6.5 < pT < 7   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 7   < pT < 7.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 7.5 < pT < 8   */
                                            {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 8   < pT < 9   */
                                            {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 9   < pT < 10  */
                                            {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 10  < pT < 12  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},   /* 12  < pT < 16  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},  /* 16  < pT < 20  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},  /* 20  < pT < 24  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},  /* 24  < pT < 36  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06},  /* 36  < pT < 50  */
                                            {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0., 10., 10., 0.06}}; /* 50  < pT < 100 */

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11",
  "pT bin 12",
  "pT bin 13",
  "pT bin 14",
  "pT bin 15",
  "pT bin 16",
  "pT bin 17",
  "pT bin 18",
  "pT bin 19",
  "pT bin 20",
  "pT bin 21",
  "pT bin 22",
  "pT bin 23",
  "pT bin 24"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "DCA", "cos theta*", "pT K", "pT Pi", "d0K", "d0pi", "d0d0", "cos pointing angle", "cos pointing angle xy", "normalized decay length XY", "decay length", "decay length XY", "minimum decay length"};
} // namespace hf_cuts_d0_topik

namespace hf_cuts_lc_topkpi
{
static constexpr int npTBins = 10;
static constexpr int nCutVars = 7;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0.,
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  24.,
  36.};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cuts[npTBins][nCutVars] = {{0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 0  < pT < 1  */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 1  < pT < 2  */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 2  < pT < 3  */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 3  < pT < 4  */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 4  < pT < 5  */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 5  < pT < 6  */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 6  < pT < 8  */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 8  < pT < 12 */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.},  /* 12 < pT < 24 */
                                            {0.400, 0.4, 0.4, 0.4, 0., 0.005, 0.}}; /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "pT p", "pT K", "pT Pi", "Chi2PCA", "decay length", "cos pointing angle"};
} // namespace hf_cuts_lc_topkpi

namespace hf_cuts_lc_tok0sp
{
static constexpr int npTBins = 8;
static constexpr int nCutVars = 8;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  24.};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
// mK0s(GeV)     mLambdas(GeV)    mGammas(GeV)    ptp     ptK0sdau     ptK0s     d0p     d0K0
constexpr double cuts[npTBins][nCutVars] = {{0.008, 0.005, 0.1, 0.5, 0.3, 0.6, 0.05, 999999.},  // 1 < pt < 2
                                            {0.008, 0.005, 0.1, 0.5, 0.4, 1.3, 0.05, 999999.},  // 2 < pt < 3
                                            {0.009, 0.005, 0.1, 0.6, 0.4, 1.3, 0.05, 999999.},  // 3 < pt < 4
                                            {0.011, 0.005, 0.1, 0.6, 0.4, 1.4, 0.05, 999999.},  // 4 < pt < 5
                                            {0.013, 0.005, 0.1, 0.6, 0.4, 1.4, 0.06, 999999.},  // 5 < pt < 6
                                            {0.013, 0.005, 0.1, 0.9, 0.4, 1.6, 0.09, 999999.},  // 6 < pt < 8
                                            {0.016, 0.005, 0.1, 0.9, 0.4, 1.7, 0.10, 999999.},  // 8 < pt < 12
                                            {0.019, 0.005, 0.1, 1.0, 0.4, 1.9, 0.20, 999999.}}; // 12 < pt < 24

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7"};

// column labels
static const std::vector<std::string> cutVarLabels = {"mK0s", "mLambda", "mGamma", "ptBach", "ptV0Dau", "ptV0", "d0Bach", "d0V0"};
} // namespace hf_cuts_lc_tok0sp

namespace hf_cuts_dplus_topikpi
{
static const int npTBins = 12;
static const int nCutVars = 8;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  7.,
  8.,
  10.,
  12.,
  16.,
  24.,
  36.};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
// selections from pp at 5 TeV 2017 analysis https://alice-notes.web.cern.ch/node/808
constexpr double cuts[npTBins][nCutVars] = {{0.2, 0.3, 0.3, 0.07, 6., 0.96, 0.985, 2.5},  /* 1  < pT < 2  */
                                            {0.2, 0.3, 0.3, 0.07, 5., 0.96, 0.985, 2.5},  /* 2  < pT < 3  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.980, 2.5},  /* 3  < pT < 4  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 4  < pT < 5  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 5  < pT < 6  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 6  < pT < 7  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 7  < pT < 8  */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 8  < pT < 10 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 10 < pT < 12 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 12 < pT < 16 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 16 < pT < 24 */
                                            {0.2, 0.3, 0.3, 0.20, 5., 0.94, 0.000, 2.5}}; /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> cutVarLabels = {"deltaM", "pT Pi", "pT K", "decay length", "normalized decay length XY", "cos pointing angle", "cos pointing angle XY", "max normalized deltaIP"};
} // namespace hf_cuts_dplus_topikpi

namespace hf_cuts_ds_tokkpi
{
static const int npTBins = 12;
static const int nCutVars = 8;
// momentary cuts
constexpr double pTBins[npTBins + 1] = {
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  7.,
  8.,
  10.,
  12.,
  16.,
  24.,
  36.};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

constexpr double cuts[npTBins][nCutVars] = {{0.2, 0.3, 0.3, 0.07, 6., 0.96, 0.985, 2.5},  /* 1  < pT < 2  */
                                            {0.2, 0.3, 0.3, 0.07, 5., 0.96, 0.985, 2.5},  /* 2  < pT < 3  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.980, 2.5},  /* 3  < pT < 4  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 4  < pT < 5  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 5  < pT < 6  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 6  < pT < 7  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 7  < pT < 8  */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 8  < pT < 10 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 10 < pT < 12 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 12 < pT < 16 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 16 < pT < 24 */
                                            {0.2, 0.3, 0.3, 0.20, 5., 0.94, 0.000, 2.5}}; /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "pT Pi", "pT K", "decay length", "normalized decay length XY", "cos pointing angle", "cos pointing angle XY", "max normalized deltaIP"};
} // namespace hf_cuts_ds_tokkpi

namespace hf_cuts_xic_topkpi
{
static const int npTBins = 10;
static const int nCutVars = 7;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0.,
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  24.,
  36.};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cuts[npTBins][nCutVars] = {{0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 0  < pT < 1  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 1  < pT < 2  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 2  < pT < 3  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 3  < pT < 4  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 4  < pT < 5  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 5  < pT < 6  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 6  < pT < 8  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 8  < pT < 12 */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8},  /* 12 < pT < 24 */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8}}; /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "pT p", "pT K", "pT Pi", "chi2PCA", "decay length", "cos pointing angle"};
} // namespace hf_cuts_xic_topkpi

namespace hf_cuts_xicc_topkpipi
{
static const int npTBins = 10;
static const int nCutVars = 14;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0.,
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  24.,
  36.};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cuts[npTBins][nCutVars] = {{0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 0  < pT < 1  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 1  < pT < 2  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 2  < pT < 3  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 3  < pT < 4  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 4  < pT < 5  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 5  < pT < 6  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 6  < pT < 8  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 8  < pT < 12 */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 12 < pT < 24 */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8}}; /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "pT Xic", "pT Pi", "min d0 Xic", "max d0 Xic", "min d0 Pi", "max d0 Pi", "d0d0", "chi2PCA", "min decay length", "max decay length", "max decay length XY", "cos pointing angle", "cos pointing angle XY"};
} // namespace hf_cuts_xicc_topkpipi

namespace hf_cuts_jpsi_toee
{
static constexpr int npTBins = 9;
static constexpr int nCutVars = 5;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  15.0,
};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
constexpr double cuts[npTBins][nCutVars] = {{0.5, 0.2, 0.4, 1, 1.},  /* 0   < pT < 0.5 */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 0.5 < pT < 1   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 1   < pT < 2   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 2   < pT < 3   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 3   < pT < 4   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 4   < pT < 5   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 5   < pT < 7   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 7   < pT < 10  */
                                            {0.5, 0.2, 0.4, 1, 1.}}; /* 10  < pT < 15  */

// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "DCA_xy", "DCA_z", "pT El", "chi2PCA"};
} // namespace hf_cuts_jpsi_toee

namespace hf_cuts_bplus_tod0pi
{
static constexpr int npTBins = 12;
static constexpr int nCutVars = 11;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  13.0,
  16.0,
  20.0,
  24.0};

auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
// DeltaM CPA d0D0 d0Pi pTD0 pTPi BDecayLength BDecayLengthXY IPProd DeltaMD0 CthetaStr
constexpr double cuts[npTBins][nCutVars] = {{1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0 < pt < 0.5 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0.5 < pt < 1 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 1 < pt < 2 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 2 < pt < 3 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 3 < pt < 4 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 4 < pt < 5 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 5 < pt < 7 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 7 < pt < 10 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 10 < pt < 13 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 13 < pt < 16 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 16 < pt < 20 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8}}; /* 20 < pt < 24 */
// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "CPA", "d0 D0", "d0 Pi", "pT D0", "pT Pi", "B decLen", "B decLenXY", "Imp. Par. Product", "DeltaMD0", "Cos ThetaStar"};
} // namespace hf_cuts_bplus_tod0pi

namespace hf_cuts_lb_tolcpi
{
static constexpr int npTBins = 12;
static constexpr int nCutVars = 12;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  13.0,
  16.0,
  20.0,
  24.0};

auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
// DeltaM CPA chi2PCA d0Lc d0Pi pTLc pTPi LbDecayLength LbDecayLengthXY IPProd DeltaMLc CthetaStr
constexpr double cuts[npTBins][nCutVars] = {{1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0 < pt < 0.5 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0.5 < pt < 1 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 1 < pt < 2 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 2 < pt < 3 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 3 < pt < 4 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 4 < pt < 5 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 5 < pt < 7 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 7 < pt < 10 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 10 < pt < 13 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 13 < pt < 16 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 16 < pt < 20 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8}}; /* 20 < pt < 24 */
// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> cutVarLabels = {"m", "CPA", "Chi2PCA", "d0 Lc+", "d0 Pi", "pT Lc+", "pT Pi", "Lb decLen", "Lb decLenXY", "Imp. Par. Product", "DeltaMLc", "Cos ThetaStar"};
} // namespace hf_cuts_lb_tolcpi

namespace hf_cuts_x_tojpsipipi
{
static constexpr int npTBins = 9;
static constexpr int nCutVars = 7;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  15.0,
};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
//                                            m   CPA  d0Jpsi  d0Pi pTJpsi pTPi chi2PCA
constexpr double cuts[npTBins][nCutVars] = {{0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 0<pt<0.5 */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 0.5<pt<1 */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 1<pt<2   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 2<pt<3   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 3<pt<4   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 4<pt<5   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 5<pt<7   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 7<pt<10  */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.}}; /* 10<pt<15 */
// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8"};
// column labels
static const std::vector<std::string> cutVarLabels = {"m", "CPA", "d0 Jpsi", "d0 Pi", "pT Jpsi", "pT Pi", "chi2PCA"};
} // namespace hf_cuts_x_tojpsipipi

namespace hf_cuts_chic_tojpsigamma
{
// dummy selections for chic --> TO BE IMPLEMENTED
static constexpr int npTBins = 9;
static constexpr int nCutVars = 7;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double pTBins[npTBins + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  15.0,
};
auto pTBins_v = std::vector<double>{pTBins, pTBins + npTBins + 1};

// default values for the cuts
//                                            m   CPA  d0Jpsi  d0gamma pTJpsi pTgamma chi2PCA
constexpr double cuts[npTBins][nCutVars] = {{3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 0<pt<0.5 */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 0.5<pt<1 */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 1<pt<2   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 2<pt<3   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 3<pt<4   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 4<pt<5   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 5<pt<7   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 7<pt<10  */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.}}; /* 10<pt<15 */
// row labels
static const std::vector<std::string> pTBinLabels = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8"};
// column labels
static const std::vector<std::string> cutVarLabels = {"m", "CPA", "d0 Jpsi", "d0 Gamma", "pT Jpsi", "pT Gamma", "chi2PCA"};
} // namespace hf_cuts_chic_tojpsigamma

} // namespace o2::analysis

#endif // HF_SELECTOR_CUTS_H_
