// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "DetectorBuild.hpp"

namespace Acts {
namespace Test {

  using id = unsigned long int;

  ///
  /// @brief This struct collects surfaces which are hit by the propagator and
  /// which carries at least one measurement
  ///
  struct SurfaceCollector
  {
    // Collection of measurements sorted by their surfaces
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    /// @brief Constructor
    SurfaceCollector() = default;

    using result_type = std::vector<Surface const*>;

    /// @brief Operater that is callable by an ActionList. The function collects
    /// the surfaces
    ///
    /// @tparam propagator_state_t Type of the propagator state
    /// @param [in] state State of the propagator
    /// @param [out] result Vector of matching surfaces
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, result_type& result) const
    {
      if (measurements.find(state.navigation.currentSurface)
          != measurements.end())
        result.push_back(state.navigation.currentSurface);
    }
  };

  ///
  /// @brief Unit test for Kalman fitter with measurements along the x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_initialization)
  {
    // Build detector
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

    // Construct measurements
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    ActsSymMatrixD<2> cov2D;
    cov2D << 1. * units::_mm, 0., 0., 1. * units::_mm;

    Vector3D       pos(-2. * units::_m, 0., 0.);
    Surface const* sur = detector->lowestTrackingVolume(pos)
                             ->associatedLayer(pos)
                             ->surfaceArray()
                             ->at(pos)[0];
    measurements[sur].push_back(
        std::move(Measurement<id, eLOC_0, eLOC_1>(*sur, 0, cov2D, 0., 0.)));
    pos = {-1. * units::_m, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 1, cov2D, 0., 0.));

    ActsSymMatrixD<1> cov1D;
    cov1D << 1. * units::_mm;

    pos = {1. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 2, cov1D, 0.));

    pos = {1. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 3, cov1D, 0.));

    pos = {2. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 4, cov1D, 0.));

    pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 5, cov1D, 0.));

    // Build navigator
    Navigator navi(detector);
    navi.resolvePassive   = true;
    navi.resolveMaterial  = false;
    navi.resolveSensitive = false;
    // Use default stepper
    StraightLineStepper sls;
    // Build navigator
    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(-3. * units::_m, 0., 0.),
        startMom(1. * units::_GeV, 0., 0);

    SingleCurvilinearTrackParameters<NeutralPolicy> sbtp(
        std::move(covPtr), startParams, startMom);

    // Create surface collection
    ActionList<SurfaceCollector> aList;
    aList.get<SurfaceCollector>().measurements = measurements;
    // Set options for propagator
    Propagator<StraightLineStepper,
               Navigator>::Options<ActionList<SurfaceCollector>>
        propOpts;
    propOpts.actionList = aList;

    // Launch and collect
    const auto& result = prop.propagate(sbtp, *sur, propOpts);
    const std::vector<Surface const*>& surResult
        = result.get<typename SurfaceCollector::result_type>();

    BOOST_TEST(surResult.size() == 6);
  }
}  // namespace Test
}  // namespace Acts
