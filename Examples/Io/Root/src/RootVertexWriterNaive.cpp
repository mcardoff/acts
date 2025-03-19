// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootVertexWriterNaive.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include "Acts/Vertexing/Vertex.hpp"

#include <algorithm>
#include <cstdint>
#include <ios>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <unordered_map>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

RootVertexWriterNaive::RootVertexWriterNaive(const RootVertexWriterNaive::Config& cfg,
                                   Acts::Logging::Level lvl)
    : WriterT(cfg.inputVertices, "RootVertexWriterNaive", lvl), m_cfg(cfg) {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (m_cfg.inputTruthVertices.empty()) {
    throw std::invalid_argument("Collection with truth vertices missing");
  }

  m_inputTruthVertices.initialize(m_cfg.inputTruthVertices);
  m_inputTruthHSTrackInfo.initialize(m_cfg.inputTruthHSTrackInfo);
  m_inputRecoVertices.initialize(m_cfg.inputRecoVertices);
  m_inputRecoHSTrackInfo.initialize(m_cfg.inputRecoHSTrackInfo);
  
  // open root file and create the tree
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // setup the branches
  m_outputTree->Branch("event_id", &m_eventId);

  // ACTS
  m_outputTree->Branch("ACTS_nVertex"     , &m_acts_nRecoVtx);
  m_outputTree->Branch("ACTS_nTracksHS"   , &m_acts_nTracksHS);
  m_outputTree->Branch("ACTS_Vertex_x"    , &m_acts_vx);
  m_outputTree->Branch("ACTS_Vertex_y"    , &m_acts_vy);
  m_outputTree->Branch("ACTS_Vertex_z"    , &m_acts_vz);
  m_outputTree->Branch("ACTS_Vertex_t"    , &m_acts_vt);
  m_outputTree->Branch("ACTS_res_x"       , &m_acts_resX);
  m_outputTree->Branch("ACTS_res_y"       , &m_acts_resY);
  m_outputTree->Branch("ACTS_res_z"       , &m_acts_resZ);
  m_outputTree->Branch("ACTS_res_t"       , &m_acts_resT);
  m_outputTree->Branch("ACTS_HS_res_x"    , &m_acts_hs_resX);
  m_outputTree->Branch("ACTS_HS_res_y"    , &m_acts_hs_resY);
  m_outputTree->Branch("ACTS_HS_res_z"    , &m_acts_hs_resZ);
  m_outputTree->Branch("ACTS_HS_res_t"    , &m_acts_hs_resT);
  m_outputTree->Branch("ACTS_dzNearest"   , &m_acts_dzNearest);
  m_outputTree->Branch("ACTS_dzNearestAbs", &m_acts_dzNearestAbs);
  m_outputTree->Branch("ACTS_HS_track_z0" , &m_acts_track_z0);
  m_outputTree->Branch("ACTS_HS_z0_pvz"   , &m_acts_z0_pvz);
  m_outputTree->Branch("ACTS_HS_track_eta", &m_acts_track_eta);
  m_outputTree->Branch("ACTS_HS_track_pt" , &m_acts_track_pt);

  // Athena
  m_outputTree->Branch("Athena_nVertex"     , &m_athe_nRecoVtx);
  m_outputTree->Branch("Athena_nTracksHS"   , &m_athe_nTracksHS);
  m_outputTree->Branch("Athena_Vertex_x"    , &m_athe_vx);
  m_outputTree->Branch("Athena_Vertex_y"    , &m_athe_vy);
  m_outputTree->Branch("Athena_Vertex_z"    , &m_athe_vz);
  m_outputTree->Branch("Athena_Vertex_t"    , &m_athe_vt);
  m_outputTree->Branch("Athena_res_x"       , &m_athe_resX);
  m_outputTree->Branch("Athena_res_y"       , &m_athe_resY);
  m_outputTree->Branch("Athena_res_z"       , &m_athe_resZ);
  m_outputTree->Branch("Athena_res_t"       , &m_athe_resT);
  m_outputTree->Branch("Athena_HS_res_x"    , &m_athe_hs_resX);
  m_outputTree->Branch("Athena_HS_res_y"    , &m_athe_hs_resY);
  m_outputTree->Branch("Athena_HS_res_z"    , &m_athe_hs_resZ);
  m_outputTree->Branch("Athena_HS_res_t"    , &m_athe_hs_resT);
  m_outputTree->Branch("Athena_dzNearest"   , &m_athe_dzNearest);
  m_outputTree->Branch("Athena_dzNearestAbs", &m_athe_dzNearestAbs);
  m_outputTree->Branch("Athena_HS_track_z0" , &m_athe_track_z0);
  m_outputTree->Branch("Athena_HS_z0_pvz"   , &m_athe_z0_pvz);
  m_outputTree->Branch("Athena_HS_track_eta", &m_athe_track_eta);
  m_outputTree->Branch("Athena_HS_track_pt" , &m_athe_track_pt);


  // Truth
  m_outputTree->Branch("Truth_nVertex"     , &m_nTrueVtx);
  m_outputTree->Branch("Truth_nTracksHS"   , &m_true_nTracksHS);
  m_outputTree->Branch("Truth_Vertex_x"    , &m_truthX);
  m_outputTree->Branch("Truth_Vertex_y"    , &m_truthY);
  m_outputTree->Branch("Truth_Vertex_z"    , &m_truthZ);
  m_outputTree->Branch("Truth_Vertex_t"    , &m_truthT);
  m_outputTree->Branch("Truth_dzNearest"   , &m_truth_dzNearest);
  m_outputTree->Branch("Truth_dzNearestAbs", &m_truth_dzNearestAbs);
  m_outputTree->Branch("Truth_HS_track_z0" , &m_truth_track_z0);
  m_outputTree->Branch("Truth_HS_z0_pvz"   , &m_truth_z0_pvz);
  m_outputTree->Branch("Truth_HS_track_eta", &m_truth_track_eta);
  m_outputTree->Branch("Truth_HS_track_pt" , &m_truth_track_pt);
  

  // Stats
  m_outputTree->Branch("chiSquared", &m_chiSquared);
  m_outputTree->Branch("numberDoF" , &m_nDoF);
}

RootVertexWriterNaive::~RootVertexWriterNaive() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootVertexWriterNaive::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote vertices to tree '" << m_cfg.treeName << "' in '"
                                       << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootVertexWriterNaive::writeT(const AlgorithmContext& ctx,
					  const VertexCollection& vertices) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Read truth vertex input collection
  const TruthVertexContainer& truthVertices = m_inputTruthVertices(ctx);
  const HardScatterTrackInfo& truthHSTrackInfo = m_inputTruthHSTrackInfo(ctx);
  
  const TruthVertexContainer& recoVertices = m_inputRecoVertices(ctx);
  const HardScatterTrackInfo& recoHSTrackInfo = m_inputRecoHSTrackInfo(ctx);

  m_acts_nRecoVtx = vertices.size();     // from ACTS
  m_athe_nRecoVtx = recoVertices.size(); // from Athena
  m_nTrueVtx = truthVertices.size();     // truth
  m_eventId = ctx.eventNumber;           // provided by ACTS

  // find ACTS HS vertex (highest sumpt2)
  auto calcSumPt2 = [](const Acts::Vertex& vtx) {
    double sumPt2 = 0;
    for (const auto& trk : vtx.tracks()) {
      if (trk.trackWeight > 0.1) {
        double pt = trk.originalParams.as<Acts::BoundTrackParameters>()
                        ->transverseMomentum();
        sumPt2 += pt * pt;
      }
    }
    return sumPt2;
  };
  
  unsigned long hs_idx;
  double highest_sumpt2 = -1;
  for (const auto& [idx, vtx] : Acts::enumerate(vertices)) {
    double sumpt2 = calcSumPt2(vtx);

    if (sumpt2 > highest_sumpt2) {
      hs_idx = idx;
      highest_sumpt2 = sumpt2;
    }
  }
  
  for (const auto& [vtxIndex, thisVtx] : Acts::enumerate(vertices)) {
    // position
    auto thisPosition = thisVtx.fullPosition();
    m_acts_vx.push_back(Acts::clampValue<float>(thisPosition.x()));
    m_acts_vy.push_back(Acts::clampValue<float>(thisPosition.y()));
    m_acts_vz.push_back(Acts::clampValue<float>(thisPosition.z()));
    m_acts_vt.push_back(Acts::clampValue<float>(thisPosition.w()));
    
    m_chiSquared.push_back(Acts::clampValue<float>(thisVtx.fitQuality().first));
    m_nDoF.push_back(Acts::clampValue<float>(thisVtx.fitQuality().second));

    //find closest (non-same) reco vertex in z
    double nearest_dz=1000000, abs_dz=1000000;
    for(const auto& [idx, thatVtx] : Acts::enumerate(vertices)) {
      if(vtxIndex == idx) {
	continue; // skip if same
      }
      auto diff     = thisPosition.z()-thatVtx.fullPosition().z();
      auto abs_diff = std::abs(diff);
      if(abs_diff < abs_dz) {
	nearest_dz=diff;
	abs_dz = abs_diff;
      }
    }
    m_acts_dzNearest.push_back(nearest_dz);
    m_acts_dzNearestAbs.push_back(abs_dz);

    auto bestTruth = truthVertices[0];
    double bestAbsDiff = std::abs(thisPosition.z()-bestTruth.z());
    for(const auto& [idx, truthVtx] : Acts::enumerate(truthVertices)) {
      auto absDiff = std::abs(thisPosition.z()-truthVtx.z());
      if(absDiff < bestAbsDiff) {
	bestAbsDiff = absDiff;
	bestTruth = truthVtx;
      }
    }

    // Now use this to fill resolutions
    double resX = Acts::clampValue<float>(thisPosition.x()-bestTruth.x());
    double resY = Acts::clampValue<float>(thisPosition.y()-bestTruth.y());
    double resZ = Acts::clampValue<float>(thisPosition.z()-bestTruth.z());
    double resT = Acts::clampValue<float>(thisPosition.w()-bestTruth.w());

    // if this is the hard scatter vertedx we found earlier
    // store in the hs resos
    if (vtxIndex == hs_idx) {
      // resolutions
      m_acts_hs_resX.push_back(resX);
      m_acts_hs_resY.push_back(resY);
      m_acts_hs_resZ.push_back(resZ);
      m_acts_hs_resT.push_back(resT);

      for(auto trk : thisVtx.tracks()) {
	auto trackParams = trk.originalParams.as<Acts::BoundTrackParameters>();
	// spatialImpactParameters returns <d0,z0>
	// fill Hard Scatter track z0
	double z0 = trackParams->spatialImpactParameters().y();
	ACTS_DEBUG("Hard Scatter Track z0: " << z0);
	// fill Hard Scatter track z0-vertex z
	double z0mpvz = z0-thisPosition.z();
	ACTS_DEBUG("Hard Scatter Track z0 minus PV z: " << z0mpvz);
	// fill Hard Scatter Track eta
	double eta = -std::log(std::tan(trackParams->theta()/2));
        ACTS_DEBUG("Hard Scatter Track eta: " << eta);
	// fill Hard Scatter Track pT
	double pt = trackParams->transverseMomentum();
	ACTS_DEBUG("Hard Scatter Track pT: " << pt);
	
	m_acts_track_z0.push_back(z0);
	m_acts_z0_pvz.push_back(z0mpvz);
	m_acts_track_eta.push_back(eta);
	m_acts_track_pt.push_back(pt);

      }
      m_acts_nTracksHS = thisVtx.tracks().size();
    }

    m_acts_resX.push_back(resX);
    m_acts_resY.push_back(resY);
    m_acts_resZ.push_back(resZ);
    m_acts_resT.push_back(resT);

  }
    
  // same thing with athena reco vertices
  for (const auto& [vtxIndex, thisVtx] : Acts::enumerate(recoVertices)) {
    // position
    m_athe_vx.push_back(Acts::clampValue<float>(thisVtx.x()));
    m_athe_vy.push_back(Acts::clampValue<float>(thisVtx.y()));
    m_athe_vz.push_back(Acts::clampValue<float>(thisVtx.z()));
    m_athe_vt.push_back(Acts::clampValue<float>(thisVtx.w()));

    //find closest (non-same) reco vertex in z
    double nearest_dz=1000000, abs_dz=1000000;
    for(const auto& [idx, thatVtx] : Acts::enumerate(recoVertices)) {
      if(vtxIndex == idx) {
	continue; // skip if same
      }
      auto diff     = thisVtx.z()-thatVtx.z();
      auto abs_diff = std::abs(diff);
      if(abs_diff < abs_dz) {
	nearest_dz=diff;
	abs_dz = abs_diff;
      }
    }
    m_athe_dzNearest.push_back(nearest_dz);
    m_athe_dzNearestAbs.push_back(abs_dz);

    // find index of truth vertex with min(reco_z-truth_z[idx])
    auto bestTruth     = truthVertices[0];
    double bestAbsDiff = std::abs(thisVtx.z()-bestTruth.z());
    for(const auto& [idx, truthVtx] : Acts::enumerate(truthVertices)) {
      auto absDiff = std::abs(thisVtx.z()-truthVtx.z());
      if(absDiff < bestAbsDiff) {
	bestAbsDiff = absDiff;
	bestTruth = truthVtx;
      }
    }

    // Now use this to fill resolutions
    double resX = Acts::clampValue<float>(thisVtx.x()-bestTruth.x());
    double resY = Acts::clampValue<float>(thisVtx.y()-bestTruth.y());
    double resZ = Acts::clampValue<float>(thisVtx.z()-bestTruth.z());
    double resT = Acts::clampValue<float>(thisVtx.w()-bestTruth.w());

    if(vtxIndex == 0) { // first vertex is HS reco vertex
      // resolutions
      m_athe_hs_resX.push_back(resX);
      m_athe_hs_resY.push_back(resY);
      m_athe_hs_resZ.push_back(resZ);
      m_athe_hs_resT.push_back(resT);

    }
    
    m_athe_resX.push_back(resX);
    m_athe_resY.push_back(resY);
    m_athe_resZ.push_back(resZ);
    m_athe_resT.push_back(resT);
  }

  for (const auto& [vtxIndex, thisVtx] : Acts::enumerate(truthVertices)) {
    m_truthX.push_back(Acts::clampValue<float>(thisVtx.x()));
    m_truthY.push_back(Acts::clampValue<float>(thisVtx.y()));
    m_truthZ.push_back(Acts::clampValue<float>(thisVtx.z()));
    m_truthT.push_back(Acts::clampValue<float>(thisVtx.w()));

    // just for shits and giggles
    double nearest_dz=1000000, abs_dz=1000000;
    for(const auto& [idx, thatVtx] : Acts::enumerate(truthVertices)) {
      if (vtxIndex == idx) {
	continue;
      }
      auto diff     = thisVtx.z()-thatVtx.z();
      auto abs_diff = std::abs(diff);
      if(abs_diff < abs_dz) {
	nearest_dz=diff;
	abs_dz = abs_diff;
      }
    }
    m_truth_dzNearest.push_back(nearest_dz);
    m_truth_dzNearestAbs.push_back(abs_dz);
  }

  for (auto [pt, eta, z0, dpvz]: Acts::zip(truthHSTrackInfo.trackPt,
					   truthHSTrackInfo.trackEta,
					   truthHSTrackInfo.trackZ0,
					   truthHSTrackInfo.dZ0PVZ)) {
    m_truth_track_pt.push_back(pt);
    m_truth_track_eta.push_back(eta);
    m_truth_track_z0.push_back(z0);
    m_truth_z0_pvz.push_back(dpvz);
  }

  for (auto [pt, eta, z0, dpvz]: Acts::zip(recoHSTrackInfo.trackPt,
					   recoHSTrackInfo.trackEta,
					   recoHSTrackInfo.trackZ0,
					   recoHSTrackInfo.dZ0PVZ)) {
    m_athe_track_pt.push_back(pt);
    m_athe_track_eta.push_back(eta);
    m_athe_track_z0.push_back(z0);
    m_athe_z0_pvz.push_back(dpvz);
  }
  
  m_athe_nTracksHS = recoHSTrackInfo.nHSTracks;
  m_true_nTracksHS = truthHSTrackInfo.nHSTracks;

  m_outputTree->Fill();
  // Clear vector-like objects
  m_acts_vx.clear();
  m_acts_vy.clear();
  m_acts_vz.clear();
  m_acts_vt.clear();
  m_acts_resX.clear();
  m_acts_resY.clear();
  m_acts_resZ.clear();
  m_acts_resT.clear();
  m_acts_hs_resX.clear();
  m_acts_hs_resY.clear();
  m_acts_hs_resZ.clear();
  m_acts_hs_resT.clear();
  m_acts_dzNearest.clear();
  m_acts_dzNearestAbs.clear();
  m_acts_track_z0.clear();
  m_acts_z0_pvz.clear();
  m_acts_track_eta.clear();
  m_acts_track_pt.clear();

  m_athe_vx.clear();
  m_athe_vy.clear();
  m_athe_vz.clear();
  m_athe_vt.clear();
  m_athe_resX.clear();
  m_athe_resY.clear();
  m_athe_resZ.clear();
  m_athe_resT.clear();
  m_athe_hs_resX.clear();
  m_athe_hs_resY.clear();
  m_athe_hs_resZ.clear();
  m_athe_hs_resT.clear();
  m_athe_dzNearest.clear();
  m_athe_dzNearestAbs.clear();
  m_athe_z0_pvz.clear();
  m_athe_track_eta.clear();
  m_athe_track_z0.clear();
  m_athe_track_pt.clear();

  
  m_truthX.clear();
  m_truthY.clear();
  m_truthZ.clear();
  m_truthT.clear();
  m_truth_dzNearest.clear();
  m_truth_dzNearestAbs.clear();
  m_truth_z0_pvz.clear();
  m_truth_track_eta.clear();
  m_truth_track_z0.clear();
  m_truth_track_pt.clear();

  m_nDoF.clear();
  m_chiSquared.clear();
  
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
