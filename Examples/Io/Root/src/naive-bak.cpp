// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthenaNTupleReaderNaive.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Root/RootUtility.hpp"

#include <cstdint>
#include <iostream>
#include <optional>
#include <stdexcept>

#include <TChain.h>

ActsExamples::RootAthenaNTupleReaderNaive::RootAthenaNTupleReaderNaive(
    const ActsExamples::RootAthenaNTupleReaderNaive::Config& config,
    Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.inputFilePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.inputTreeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Track Parameters
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);

  // Truth Vertex Parameters
  m_outputTruthVtxParameters.initialize(m_cfg.outputTruthVtxParameters);
  m_outputTruthVtxNTracks.initialize(m_cfg.outputTruthVtxNTracks);

  // Reco Vertex Parameters
  m_outputRecoVtxParameters.initialize(m_cfg.outputRecoVtxParameters);
  m_outputRecoVtxNTracks.initialize(m_cfg.outputRecoVtxNTracks);

  // Athena Vertex Classification
  m_outputTruthVtxHS.initialize(m_cfg.outputTruthVertexIsHS);
  m_outputRecoVtxHS.initialize(m_cfg.outputRecoVertexIsHS);

  m_inputChain = std::make_unique<TChain>(m_cfg.inputTreeName.c_str());

  // unused event identifier
  int eventNumber = 0;
  // ACTS_INFO("uh ohhh, stinkyyyy");

  // Set the branches
  m_inputChain->SetBranchAddress("eventNumber", &eventNumber);
  m_inputChain->SetBranchAddress("Track_d0", &m_branches.track_d0);
  m_inputChain->SetBranchAddress("Track_z0", &m_branches.track_z0);
  m_inputChain->SetBranchAddress("Track_theta", &m_branches.track_theta);
  m_inputChain->SetBranchAddress("Track_phi", &m_branches.track_phi);
  m_inputChain->SetBranchAddress("Track_qOverP", &m_branches.track_qOverP);
  m_inputChain->SetBranchAddress("Track_time", &m_branches.track_t);
  // m_inputChain->SetBranchAddress("Track_z", &m_branches.track_z);
  ACTS_INFO("Set track vars, event num");

  // Covariance stuff
  m_inputChain->SetBranchAddress("Track_var_d0",          &m_branches.track_var_d0);
  m_inputChain->SetBranchAddress("Track_var_z0",          &m_branches.track_var_z0);
  m_inputChain->SetBranchAddress("Track_var_phi0",        &m_branches.track_var_phi);
  m_inputChain->SetBranchAddress("Track_var_theta",       &m_branches.track_var_theta);
  m_inputChain->SetBranchAddress("Track_var_qOverP",      &m_branches.track_var_qOverP);
  m_inputChain->SetBranchAddress("Track_cov_d0z0",        &m_branches.track_cov_d0z0);
  m_inputChain->SetBranchAddress("Track_cov_d0phi0",      &m_branches.track_cov_d0phi);
  m_inputChain->SetBranchAddress("Track_cov_d0theta",     &m_branches.track_cov_d0theta);
  m_inputChain->SetBranchAddress("Track_cov_d0qOverP",    &m_branches.track_cov_d0qOverP);
  m_inputChain->SetBranchAddress("Track_cov_z0phi0",      &m_branches.track_cov_z0phi);
  m_inputChain->SetBranchAddress("Track_cov_z0theta",     &m_branches.track_cov_z0theta);
  m_inputChain->SetBranchAddress("Track_cov_z0qOverP",    &m_branches.track_cov_z0qOverP);
  m_inputChain->SetBranchAddress("Track_cov_phi0theta",   &m_branches.track_cov_phitheta);
  m_inputChain->SetBranchAddress("Track_cov_phi0qOverP",  &m_branches.track_cov_phiqOverP);
  m_inputChain->SetBranchAddress("Track_cov_thetaqOverP", &m_branches.track_cov_tehtaqOverP);
  ACTS_INFO("Set COV matrix");

  // Truth vertex
  m_inputChain->SetBranchAddress("TruthVtx_x", &m_branches.truthvertex_x);
  m_inputChain->SetBranchAddress("TruthVtx_y", &m_branches.truthvertex_y);
  m_inputChain->SetBranchAddress("TruthVtx_z", &m_branches.truthvertex_z);
  m_inputChain->SetBranchAddress("TruthVtx_time", &m_branches.truthvertex_t);
  
  m_inputChain->SetBranchAddress("TruthVtx_isHS", &m_branches.truthvertex_isHS);
  m_inputChain->SetBranchAddress("TruthVtx_track_idx", &m_branches.truthvertex_tracks_idx);

  // Reco vertex
  m_inputChain->SetBranchAddress("RecoVtx_x", &m_branches.recovertex_x);
  m_inputChain->SetBranchAddress("RecoVtx_y", &m_branches.recovertex_y);
  m_inputChain->SetBranchAddress("RecoVtx_z", &m_branches.recovertex_z);
  m_inputChain->SetBranchAddress("RecoVtx_isHS", &m_branches.recovertex_isHS);
  m_inputChain->SetBranchAddress("RecoVtx_track_idx", &m_branches.recovertex_tracks_idx);

  ACTS_INFO("Set vertex info");
  
  auto path = m_cfg.inputFilePath;

  // add file to the input chain
  m_inputChain->Add(path.c_str());
  ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.inputTreeName
                            << "'.");

  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

ActsExamples::RootAthenaNTupleReaderNaive::~RootAthenaNTupleReaderNaive() = default;

ActsExamples::ProcessCode ActsExamples::RootAthenaNTupleReaderNaive::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read track parameters from ntuple.");

  Acts::Vector3 pos(0, 0, 0);
  std::shared_ptr<Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(pos);

  if (context.eventNumber >= m_events) {
    ACTS_ERROR("event out of bounds");
    return ProcessCode::ABORT;
  }

  std::lock_guard<std::mutex> lock(m_read_mutex);

  int entry = static_cast<int>(context.eventNumber);
  
  // std::cout << "m_inputChain pointer: " << m_inputChain << std::endl;
  // std::cout << "Total entries in chain: " << m_inputChain->GetEntries() << std::endl;
  // std::cout << "Trying to access entry: " << entry << std::endl;
  // std::cout << "Tree inside chain: " << m_inputChain->GetTree() << std::endl;
  // std::cout << "Attempting GetEntry now..." << std::endl;

  m_inputChain->GetEntry(entry);
  ACTS_INFO("Reading event: " << context.eventNumber
                              << " stored as entry: " << entry);

  const unsigned int nTracks = m_branches.track_d0->size();
  const unsigned int nTruthVtx = m_branches.truthvertex_z->size();
  const unsigned int nRecoVtx = m_branches.recovertex_z->size();

  ACTS_DEBUG("nTracks = " << nTracks);
  ACTS_DEBUG("nTruthVtx = " << nTruthVtx);
  ACTS_DEBUG("nRecoVtx = " << nRecoVtx);

  TrackParametersContainer trackContainer;
  trackContainer.reserve(nTracks);

  for (unsigned int i = 0; i < nTracks; i++) {
    Acts::BoundVector params;

    params[Acts::BoundIndices::eBoundLoc0] = (*m_branches.track_d0)[i];
    params[Acts::BoundIndices::eBoundLoc1] = (*m_branches.track_z0)[i];
    params[Acts::BoundIndices::eBoundPhi] = (*m_branches.track_phi)[i];
    params[Acts::BoundIndices::eBoundTheta] = (*m_branches.track_theta)[i];
    params[Acts::BoundIndices::eBoundQOverP] = (*m_branches.track_qOverP)[i];
    params[Acts::BoundIndices::eBoundTime] = (*m_branches.track_t)[i];

    // Construct and fill covariance matrix
    Acts::BoundSquareMatrix cov;

    // Variances
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.track_var_d0)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.track_var_z0)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.track_var_phi)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.track_var_theta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.track_var_qOverP)[i];

    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.track_cov_d0z0)[i];
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.track_cov_d0phi)[i];
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.track_cov_d0theta)[i];
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.track_cov_d0qOverP)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.track_cov_z0phi)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.track_cov_z0theta)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.track_cov_z0qOverP)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.track_cov_phitheta)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.track_cov_phiqOverP)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.track_cov_tehtaqOverP)[i];

    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.track_cov_d0z0)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.track_cov_d0phi)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.track_cov_d0theta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.track_cov_d0qOverP)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.track_cov_z0phi)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.track_cov_z0theta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.track_cov_z0qOverP)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.track_cov_phitheta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.track_cov_phiqOverP)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.track_cov_tehtaqOverP)[i];

    // TODO we do not have a hypothesis at hand here. defaulting to pion
    Acts::BoundTrackParameters tc(surface, params, cov,
                                  Acts::ParticleHypothesis::pion());
    trackContainer.push_back(tc);
  }

  std::vector<Acts::Vector4> truthVertexContainer;
  std::vector<int> nTracksTruth;
  std::vector<bool> truthIsHS;
  for (unsigned int i = 0; i < nTruthVtx; i++) {
    Acts::Vector4 vtx((*m_branches.truthvertex_x)[i], (*m_branches.truthvertex_y)[i],
                      (*m_branches.truthvertex_z)[i], (*m_branches.truthvertex_t)[i]);
    
    truthVertexContainer.push_back(vtx);
    // nTracksTruth.push_back((*m_branches.truthvertex_tracks_idx)[i].size());
    // truthIsHS.push_back((*m_branches.truthvertex_isHS)[i]);
  }
  std::vector<Acts::Vector4> recoVertexContainer;
  std::vector<int> nTracksReco;
  std::vector<bool> recoIsHS;
  for (unsigned int i = 0; i < nRecoVtx; i++) {
    Acts::Vector4 vtx((*m_branches.recovertex_x)[i], (*m_branches.recovertex_y)[i],
                      (*m_branches.recovertex_z)[i], 0);
    // FIXME: add time if there is a valid time for the vertex maybe
    
    recoVertexContainer.push_back(vtx);
    // nTracksReco.push_back((*m_branches.recovertex_tracks_idx)[i].size());
    // recoIsHS.push_back((*m_branches.recovertex_isHS)[i]);
  }

  // Track Parameters
  m_outputTrackParameters(context, std::move(trackContainer));

  // Truth Vertex Parameters
  m_outputTruthVtxParameters(context, std::move(truthVertexContainer));
  m_outputRecoVtxParameters(context, std::move(recoVertexContainer));

  // Reco Vertex Parameters
  // m_outputTruthVtxNTracks(context, std::move(nTracksTruth));
  // m_outputRecoVtxNTracks(context, std::move(nTracksReco));

  // Athena Vertex Classification
  // m_outputTruthVtxHS(context, std::move(truthIsHS));
  // m_outputRecoVtxHS(context, std::move(recoIsHS));

  // Return success flag
  return ProcessCode::SUCCESS;
}
