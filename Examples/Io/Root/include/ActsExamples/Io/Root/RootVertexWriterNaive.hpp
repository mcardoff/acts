// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include "Acts/Vertexing/Vertex.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

using VertexCollection = std::vector<Acts::Vertex>;
using TruthVertexContainer = std::vector<Acts::Vector4>;

struct HardScatterTrackInfo {
  // need to be vectors, multiple tracks per event
  // assuming one TrackInfoContainer per event
  std::vector<float> trackPt;
  std::vector<float> trackEta;
  std::vector<float> trackZ0;
  std::vector<float> dZ0PVZ;
  // one hs vertex per event, only need one int
  int nHSTracks = 0;
};

namespace ActsExamples {
struct AlgorithmContext;

/// Write out vertices as a flat TTree.
///
/// Each entry in the TTree corresponds to one vertex for optimum writing
/// speed. The event number is part of the written data.
///
/// Safe to use from multiple writer threads. To avoid thread-saftey issues,
/// the writer must be the sole owner of the underlying file. Thus, the
/// output file pointer can not be given from the outside.
class RootVertexWriterNaive final : public WriterT<VertexCollection> {
 public:
  struct Config {
    /// Input vertex collection to write.
    std::string inputVertices;
    /// Input truth vertex collection
    std::string inputTruthVertices;
    /// Input info on Tracks in truth HS Vertex
    std::string inputTruthHSTrackInfo;
    /// Input reco vertex collection (athena)
    std::string inputRecoVertices;
    /// Input info on Tracks in reco HS Vertex
    std::string inputRecoHSTrackInfo;
    /// Path to the output file.
    std::string filePath;
    /// Output file access mode.
    std::string fileMode = "RECREATE";
    /// Name of the tree within the output file.
    std::string treeName = "vertices";
  };

  /// Construct the vertex writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  RootVertexWriterNaive(const Config& cfg, Acts::Logging::Level lvl);

  /// Ensure underlying file is closed.
  ~RootVertexWriterNaive() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] vertices are the vertices to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const VertexCollection& vertices) override;

 private:
  Config m_cfg;

  std::mutex m_writeMutex;

  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;
  
  ReadDataHandle<TruthVertexContainer> m_inputTruthVertices{
    this, "InputTruthVertices"};
  ReadDataHandle<TruthVertexContainer> m_inputRecoVertices{
    this, "InputRecoVertices"};
  ReadDataHandle<HardScatterTrackInfo> m_inputRecoHSTrackInfo{
    this, "InputRecoTrackHSInfo"};
  ReadDataHandle<HardScatterTrackInfo> m_inputTruthHSTrackInfo{
    this, "InputTruthTrackHSInfo"};

  /// Event identifier.
  std::uint32_t m_eventId = 0;

  // ACTS vertex info
  int m_acts_nRecoVtx = -1;
  int m_acts_nTracksHS = -1;
  std::vector<float> m_acts_vx;
  std::vector<float> m_acts_vy;
  std::vector<float> m_acts_vz;
  std::vector<float> m_acts_vt;
  std::vector<float> m_acts_resX;
  std::vector<float> m_acts_resY;
  std::vector<float> m_acts_resZ;
  std::vector<float> m_acts_resT;
  std::vector<float> m_acts_hs_resX;
  std::vector<float> m_acts_hs_resY;
  std::vector<float> m_acts_hs_resZ;
  std::vector<float> m_acts_hs_resT;
  std::vector<float> m_acts_dzNearest;
  std::vector<float> m_acts_dzNearestAbs;
  std::vector<float> m_acts_track_z0;
  std::vector<float> m_acts_z0_pvz;
  std::vector<float> m_acts_track_eta;
  std::vector<float> m_acts_track_pt;

  // Athena vertex info
  int m_athe_nRecoVtx = -1;
  int m_athe_nTracksHS = -1;
  std::vector<float> m_athe_vx;
  std::vector<float> m_athe_vy;
  std::vector<float> m_athe_vz;
  std::vector<float> m_athe_vt;
  std::vector<float> m_athe_resX;
  std::vector<float> m_athe_resY;
  std::vector<float> m_athe_resZ;
  std::vector<float> m_athe_resT;
  std::vector<float> m_athe_hs_resX;
  std::vector<float> m_athe_hs_resY;
  std::vector<float> m_athe_hs_resZ;
  std::vector<float> m_athe_hs_resT;
  std::vector<float> m_athe_dzNearest;
  std::vector<float> m_athe_dzNearestAbs;
  std::vector<float> m_athe_track_z0;
  std::vector<float> m_athe_z0_pvz;
  std::vector<float> m_athe_track_eta;
  std::vector<float> m_athe_track_pt;

  // Truth Vertex info
  int m_nTrueVtx = -1;
  int m_true_nTracksHS;
  std::vector<float> m_truthX;
  std::vector<float> m_truthY;
  std::vector<float> m_truthZ;
  std::vector<float> m_truthT;
  std::vector<float> m_truth_dzNearest;
  std::vector<float> m_truth_dzNearestAbs;
  std::vector<float> m_truth_track_z0;
  std::vector<float> m_truth_z0_pvz;
  std::vector<float> m_truth_track_eta;
  std::vector<float> m_truth_track_pt;

  // Stats
  std::vector<float> m_chiSquared;
  std::vector<float> m_nDoF;
};

}  // namespace ActsExamples
