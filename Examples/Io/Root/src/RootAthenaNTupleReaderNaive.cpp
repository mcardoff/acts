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

#include <cmath>
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

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_outputTruthVtxParameters.initialize(m_cfg.outputTruthVtxParameters);
  m_outputRecoVtxParameters.initialize(m_cfg.outputRecoVtxParameters);
  
  m_outputRecoHSTrackInfo.initialize(m_cfg.outputRecoHSTrackInfo);
  m_outputTruthHSTrackInfo.initialize(m_cfg.outputTruthHSTrackInfo);

  m_inputChain = std::make_unique<TChain>(m_cfg.inputTreeName.c_str());

  // unused event identifier
  std::int32_t eventNumber = 0;

  // Set the branches
  m_inputChain->SetBranchAddress("eventNumber", &eventNumber);
  m_inputChain->SetBranchAddress("Track_d0", &m_branches.track_d0);
  m_inputChain->SetBranchAddress("Track_z0", &m_branches.track_z0);
  m_inputChain->SetBranchAddress("Track_theta", &m_branches.track_theta);
  m_inputChain->SetBranchAddress("Track_phi", &m_branches.track_phi);
  m_inputChain->SetBranchAddress("Track_qOverP", &m_branches.track_qOverP);
  m_inputChain->SetBranchAddress("Track_time", &m_branches.track_t);

  m_inputChain->SetBranchAddress("Track_eta", &m_branches.track_eta);
  m_inputChain->SetBranchAddress("Track_pt", &m_branches.track_pt);
  // m_inputChain->SetBranchAddress("Track_z", &m_branches.track_z);
  ACTS_INFO("Set track vars, event num");

  // Covariance stuff
  m_inputChain->SetBranchAddress("Track_var_d0", &m_branches.track_var_d0);
  m_inputChain->SetBranchAddress("Track_var_z0", &m_branches.track_var_z0);
  m_inputChain->SetBranchAddress("Track_var_phi0", &m_branches.track_var_phi);
  m_inputChain->SetBranchAddress("Track_var_theta",
                                 &m_branches.track_var_theta);
  m_inputChain->SetBranchAddress("Track_var_qOverP",
                                 &m_branches.track_var_qOverP);
  m_inputChain->SetBranchAddress("Track_cov_d0z0", &m_branches.track_cov_d0z0);
  m_inputChain->SetBranchAddress("Track_cov_d0phi0",
                                 &m_branches.track_cov_d0phi);
  m_inputChain->SetBranchAddress("Track_cov_d0theta",
                                 &m_branches.track_cov_d0theta);
  m_inputChain->SetBranchAddress("Track_cov_d0qOverP",
                                 &m_branches.track_cov_d0qOverP);
  m_inputChain->SetBranchAddress("Track_cov_z0phi0",
                                 &m_branches.track_cov_z0phi);
  m_inputChain->SetBranchAddress("Track_cov_z0theta",
                                 &m_branches.track_cov_z0theta);
  m_inputChain->SetBranchAddress("Track_cov_z0qOverP",
                                 &m_branches.track_cov_z0qOverP);
  m_inputChain->SetBranchAddress("Track_cov_phi0theta",
                                 &m_branches.track_cov_phitheta);
  m_inputChain->SetBranchAddress("Track_cov_phi0qOverP",
                                 &m_branches.track_cov_phiqOverP);
  m_inputChain->SetBranchAddress("Track_cov_thetaqOverP",
                                 &m_branches.track_cov_tehtaqOverP);
  ACTS_INFO("Set COV matrix");

  // Truth vertex
  m_inputChain->SetBranchAddress("TruthVtx_x", &m_branches.truthvertex_x);
  m_inputChain->SetBranchAddress("TruthVtx_y", &m_branches.truthvertex_y);
  m_inputChain->SetBranchAddress("TruthVtx_z", &m_branches.truthvertex_z);
  m_inputChain->SetBranchAddress("TruthVtx_time", &m_branches.truthvertex_t);

  m_inputChain->SetBranchAddress("RecoVtx_x", &m_branches.recovertex_x);
  m_inputChain->SetBranchAddress("RecoVtx_y", &m_branches.recovertex_y);
  m_inputChain->SetBranchAddress("RecoVtx_z", &m_branches.recovertex_z);

  m_inputChain->SetBranchAddress("Track_truthVtx_idx", &m_branches.track_truthvtx_idx);
  m_inputChain->SetBranchAddress("Track_recoVtx_idx",  &m_branches.track_recovtx_idx);

  ACTS_DEBUG("!!!!!!!!SET NEW :3");

  ACTS_INFO("Set vertex info");
  
  auto path = m_cfg.inputFilePath;

  // add file to the input chain
  m_inputChain->Add(path.c_str());
  // ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.inputTreeName
                            // << "'.");

  m_events = m_inputChain->GetEntries();
  // ACTS_DEBUG("The full chain has " << m_events << " entries.");
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

  auto entry = context.eventNumber;
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
  HardScatterTrackInfo truthHSInfo;
  HardScatterTrackInfo recoHSInfo;

  for (unsigned int i = 0; i < nTracks; i++) {
    Acts::BoundVector params;

    params[Acts::BoundIndices::eBoundLoc0] = (*m_branches.track_d0)[i];
    params[Acts::BoundIndices::eBoundLoc1] = (*m_branches.track_z0)[i];
    params[Acts::BoundIndices::eBoundPhi] = (*m_branches.track_phi)[i];
    params[Acts::BoundIndices::eBoundTheta] = (*m_branches.track_theta)[i];
    params[Acts::BoundIndices::eBoundQOverP] = (*m_branches.track_qOverP)[i];
    params[Acts::BoundIndices::eBoundTime] = (*m_branches.track_t)[i];
    ACTS_DEBUG("Set Parameters");

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
    ACTS_DEBUG("Set Variances");

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
    ACTS_DEBUG("Set Covariances");

    // Need to implement the following track selections
    // Pixel Hits >= 3
    // Pixel Holes <= 1
    // Strip Hits >= 0
    // Tot. Si Hits >= 7
    // pT > 900 MeV
    // |d_0| < 1.0 mm
    // sigma(d_0) < 0.35 mm
    // sigma(z_0sin(theta)) < 2.5 mm


    double theta = (*m_branches.track_var_theta)[i];
    double sinsqth = std::sin(theta)*std::sin(theta);
    double varz0 = (*m_branches.track_var_z0)[i]*(*m_branches.track_var_z0)[i];
    double varth = (*m_branches.track_var_theta)[i]*(*m_branches.track_var_theta)[i];
    int nPixelHits = (*m_branches.track_hits_pixel)[i];
    int nStripHits (*m_branches.track_hits_strip)[i];
    int nSiHits = nPixelHits + nStripHits;
    int nPixelHole = (*m_branches.track_hole_strip)[i];
    double pT = (*m_branches.track_pT)[i];
    double d0 = (*m_branches.track_d0)[i];
    double z0 = (*m_branches.track_z0)[i];
    double sd0 = std::sqrt((*m_branches.track_var_z0)[i]);
    double sz0sin = std::sqrt(sinsqth*varz0 + z0*z0*(1-sinsqth)*varth);

    bool pixel_hits_cut = nPixelHits > 3;
    bool pixel_hole_cut = nPixelHole <= 1;
    bool strip_hits_cut = nStripHits >= 0;
    bool si_hits_cut = nSiHits >= 7;
    bool pt_cut = pT > 0.9;
    bool d0_cut = std::abs(d0) < 1.0;
    bool sd0_cut = sd0 < 0.35;
    bool sz0_cut = sz0sin < 2.5;
    bool all_cuts = pixel_hits_cut && pixel_hole_cut && strip_hits_cut && si_hits_cut && pt_cut && d0_cut && sd0_cut && sz0_cut;
    
      
    if (all_cuts) {
      // TODO we do not have a hypothesis at hand here. defaulting to pion
      Acts::BoundTrackParameters tc(surface, params, cov,
				    Acts::ParticleHypothesis::pion());
      trackContainer.push_back(tc);
    }
    // already looping over tracks, might as well extract information

    // if this track ~ truth HS Vertex
    if((*m_branches.track_truthvtx_idx)[i] == 0) { 
      truthHSInfo.trackPt.push_back((*m_branches.track_pt)[i]);
      truthHSInfo.trackEta.push_back((*m_branches.track_eta)[i]);
      truthHSInfo.trackZ0.push_back((*m_branches.track_z0)[i]);
      truthHSInfo.dZ0PVZ.push_back((*m_branches.track_z0)[i]
				   -(*m_branches.truthvertex_z)[0]);
      truthHSInfo.nHSTracks++;
    }
    ACTS_DEBUG("Set Truth HS Info");

    // if this track ~ reco HS Vertex
    if((*m_branches.track_recovtx_idx)[i] == 0) { 
      recoHSInfo.trackPt.push_back((*m_branches.track_pt)[i]);
      recoHSInfo.trackEta.push_back((*m_branches.track_eta)[i]);
      recoHSInfo.trackZ0.push_back((*m_branches.track_z0)[i]);
      recoHSInfo.dZ0PVZ.push_back((*m_branches.track_z0)[i]
				  -(*m_branches.recovertex_z)[0]);
      recoHSInfo.nHSTracks++;
    }
    ACTS_DEBUG("Set Reco HS Info");
    
  }

  // truth information
  std::vector<Acts::Vector4> truthVertexContainer;
  for (unsigned int i = 0; i < nTruthVtx; i++) {
    Acts::Vector4 vtx((*m_branches.truthvertex_x)[i], (*m_branches.truthvertex_y)[i],
                      (*m_branches.truthvertex_z)[i], (*m_branches.truthvertex_t)[i]);
    truthVertexContainer.push_back(vtx);
  }
  ACTS_DEBUG("Set truth vertex positions");
  
  std::vector<Acts::Vector4> recoVertexContainer;
  for (unsigned int i = 0; i < nRecoVtx; i++) {
    Acts::Vector4 vtx((*m_branches.recovertex_x)[i], (*m_branches.recovertex_y)[i],
                      (*m_branches.recovertex_z)[i], 0);
    recoVertexContainer.push_back(vtx);
  }
  ACTS_DEBUG("Set reco vertex positions");

  // see if we can calculate it this way
  ACTS_DEBUG("Found " << recoHSInfo.nHSTracks << " tracks in reco HS Vertex");
  ACTS_DEBUG("Found " << truthHSInfo.nHSTracks << " tracks in truth HS Vertex");

  m_outputTrackParameters(context, std::move(trackContainer));
  m_outputTruthVtxParameters(context, std::move(truthVertexContainer));
  m_outputTruthHSTrackInfo(context, std::move(truthHSInfo));
  m_outputRecoVtxParameters(context, std::move(recoVertexContainer));
  m_outputRecoHSTrackInfo(context, std::move(recoHSInfo));
  
  // m_outputNRecoTracksHS(context, std::move(nRecoTracksHS));
  // m_outputNTruthTracksHS(context, std::move(nTruthTracksHS));

  // Return success flag
  return ProcessCode::SUCCESS;
}
