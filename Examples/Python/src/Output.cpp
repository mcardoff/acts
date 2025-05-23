// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "ActsExamples/Io/Csv/CsvBFieldWriter.hpp"
#include "ActsExamples/Io/Csv/CsvExaTrkXGraphWriter.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Csv/CsvProtoTrackWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSeedWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointsBucketWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackParameterWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ActsExamples/Io/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Obj/ObjSimHitWriter.hpp"
#include "ActsExamples/Io/Obj/ObjTrackingGeometryWriter.hpp"
#include "ActsExamples/Io/Root/RootBFieldWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Io/Root/RootNuclearInteractionParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationSummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootSeedWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootSpacepointWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootVertexWriter.hpp"
// #include "ActsExamples/Io/Root/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/TrackFinderNTupleWriter.hpp"
#include "ActsExamples/Io/Root/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/VertexNTupleWriter.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupReader.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupWriter.hpp"

#include "ActsExamples/Io/Root/RootVertexWriterNaive.hpp"

#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace Acts {
class TrackingGeometry;
namespace detail {
struct Step;
}  // namespace detail
}  // namespace Acts
namespace ActsExamples {
class IWriter;
struct AlgorithmContext;
}  // namespace ActsExamples

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace {
template <ActsExamples::CsvBFieldWriter::CoordinateType CType, bool Grid>
void register_csv_bfield_writer_binding(
    pybind11::class_<ActsExamples::CsvBFieldWriter>& w) {
  std::string name =
      std::string(CType == ActsExamples::CsvBFieldWriter::CoordinateType::XYZ
                      ? "Xyz"
                      : "Rz") +
      std::string(Grid ? "Grid" : "Gridless");

  using Config = ActsExamples::CsvBFieldWriter::Config<CType, Grid>;
  w.def_static((std::string("run") + name).c_str(),
               [](const Config& config, Acts::Logging::Level level) {
                 ActsExamples::CsvBFieldWriter::run(config, level);
               },
               py::arg("config"), py::arg("level"));
  auto c = py::class_<Config>(w, (std::string("Config") + name).c_str())
               .def(py::init<>());
  ACTS_PYTHON_STRUCT(c, fileName, bField, range, bins);
}
}  // namespace

namespace Acts::Python {

void addOutput(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::ObjPropagationStepsWriter, mex,
                             "ObjPropagationStepsWriter", collection, outputDir,
                             outputScalor, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::ObjSimHitWriter, mex,
                             "ObjSimHitWriter", inputSimHits, outputDir,
                             outputStem, outputPrecision, drawConnections,
                             momentumThreshold, momentumThresholdTraj,
                             nInterpolatedPoints, keepOriginalHits);

  {
    auto c = py::class_<ViewConfig>(m, "ViewConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, visible, color, offset, lineThickness,
                       surfaceThickness, quarterSegments, triangulate,
                       outputName);

    patchKwargsConstructor(c);

    py::class_<Color>(m, "Color")
        .def(py::init<>())
        .def(py::init<int, int, int>())
        .def(py::init<double, double, double>())
        .def(py::init<std::string_view>())
        .def_readonly("rgb", &Color::rgb);
  }

  py::class_<IVisualization3D>(m, "IVisualization3D")
      .def("write", py::overload_cast<const std::filesystem::path&>(
                        &IVisualization3D::write, py::const_));

  {
    using Writer = ActsExamples::ObjTrackingGeometryWriter;
    auto w = py::class_<Writer, std::shared_ptr<Writer>>(
                 mex, "ObjTrackingGeometryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const AlgorithmContext&,
                                                 const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, outputScalor, outputPrecision, outputDir,
                       containerView, volumeView, sensitiveView, passiveView,
                       gridView);
  }

  // Bindings for the binning in e.g., TrackFinderPerformanceWriter
  {
    py::class_<PlotHelpers::Binning>(mex, "Binning")
        .def(py::init<std::string, int, double, double>(), "title"_a, "bins"_a,
             "bMin"_a, "bMax"_a)
        .def(py::init<std::string, std::vector<double>>(), "title"_a, "bins"_a);

    py::class_<EffPlotTool::Config>(mex, "EffPlotToolConfig")
        .def(py::init<std::map<std::string, PlotHelpers::Binning>>(),
             "varBinning"_a);

    py::class_<FakeRatePlotTool::Config>(mex, "FakeRatePlotToolConfig")
        .def(py::init<std::map<std::string, PlotHelpers::Binning>>(),
             "varBinning"_a);

    py::class_<DuplicationPlotTool::Config>(mex, "DuplicationPlotToolConfig")
        .def(py::init<std::map<std::string, PlotHelpers::Binning>>(),
             "varBinning"_a);
  }

  // ROOT WRITERS
  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootPropagationStepsWriter, mex,
                             "RootPropagationStepsWriter", collection, filePath,
                             fileMode);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootPropagationSummaryWriter, mex,
                             "RootPropagationSummaryWriter",
                             inputSummaryCollection, filePath, fileMode);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootParticleWriter, mex,
                             "RootParticleWriter", inputParticles, filePath,
                             fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootVertexWriter, mex,
                             "RootVertexWriter", inputVertices, filePath,
                             fileMode, treeName);
  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootVertexWriterNaive, mex,
                             "RootVertexWriterNaive", inputVertices,
			     inputTruthVertices, inputRecoVertices,
			     inputRecoHSTrackInfo, inputTruthHSTrackInfo,
			     filePath, fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::TrackFinderNTupleWriter, mex,
                             "TrackFinderNTupleWriter", inputTracks,
                             inputParticles, inputParticleMeasurementsMap,
                             inputTrackParticleMatching, filePath, fileMode,
                             treeNameTracks, treeNameParticles);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::TrackFitterPerformanceWriter, mex,
                             "TrackFitterPerformanceWriter", inputTracks,
                             inputParticles, inputTrackParticleMatching,
                             filePath, resPlotToolConfig, effPlotToolConfig,
                             trackSummaryPlotToolConfig);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::RootTrackParameterWriter, mex, "RootTrackParameterWriter",
      inputTrackParameters, inputProtoTracks, inputParticles, inputSimHits,
      inputMeasurementParticlesMap, inputMeasurementSimHitsMap, filePath,
      treeName, fileMode);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::RootMaterialTrackWriter, mex, "RootMaterialTrackWriter",
      inputMaterialTracks, filePath, fileMode, treeName, recalculateTotals,
      prePostStep, storeSurface, storeVolume, collapseInteractions);

  {
    using Writer = ActsExamples::RootBFieldWriter;
    auto w =
        py::class_<Writer>(mex, "RootBFieldWriter")
            .def_static(
                "run",
                [](const Writer::Config& config, Acts::Logging::Level level) {
                  Writer::run(config, Acts::getDefaultLogger("RootBFieldWriter",
                                                             level));
                },
                py::arg("config"), py::arg("level"));

    py::enum_<Writer::GridType>(w, "GridType")
        .value("rz", Writer::GridType::rz)
        .value("xyz", Writer::GridType::xyz);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, treeName, fileName, fileMode, bField, gridType,
                       rBounds, zBounds, rBins, zBins, phiBins);
  }

  {
    using Writer = ActsExamples::RootMeasurementWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootMeasurementWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, inputMeasurements, inputClusters, inputSimHits,
                       inputMeasurementSimHitsMap, filePath, fileMode,
                       surfaceByIdentifier);
  }

  py::class_<IMaterialWriter, std::shared_ptr<IMaterialWriter>>(
      mex, "IMaterialWriter");

  py::class_<ActsExamples::ITrackParamsLookupWriter,
             std::shared_ptr<ActsExamples::ITrackParamsLookupWriter>>(
      mex, "ITrackParamsLookupWriter");

  py::class_<ActsExamples::ITrackParamsLookupReader,
             std::shared_ptr<ActsExamples::ITrackParamsLookupReader>>(
      mex, "ITrackParamsLookupReader");

  {
    using Writer = ActsExamples::RootMaterialWriter;
    auto w = py::class_<Writer, IMaterialWriter, std::shared_ptr<Writer>>(
                 mex, "RootMaterialWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, processSensitives, processApproaches,
                       processRepresenting, processBoundaries, processVolumes,
                       folderSurfaceNameBase, folderVolumeNameBase, voltag,
                       boutag, laytag, apptag, sentag, ntag, vtag, otag, mintag,
                       maxtag, ttag, x0tag, l0tag, atag, ztag, rhotag, filePath,
                       fileMode);
  }

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootSeedWriter, mex,
                             "RootSeedWriter", inputSeeds, writingMode,
                             filePath, fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootSimHitWriter, mex,
                             "RootSimHitWriter", inputSimHits, filePath,
                             fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::RootSpacepointWriter, mex,
                             "RootSpacepointWriter", inputSpacepoints, filePath,
                             fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::RootTrackStatesWriter, mex, "RootTrackStatesWriter",
      inputTracks, inputParticles, inputTrackParticleMatching, inputSimHits,
      inputMeasurementSimHitsMap, filePath, treeName, fileMode);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::RootTrackSummaryWriter, mex, "RootTrackSummaryWriter",
      inputTracks, inputParticles, inputTrackParticleMatching, filePath,
      treeName, fileMode, writeCovMat, writeGsfSpecific, writeGx2fSpecific);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::VertexNTupleWriter, mex, "VertexNTupleWriter",
      inputVertices, inputTracks, inputTruthVertices, inputParticles,
      inputSelectedParticles, inputTrackParticleMatching, bField, filePath,
      treeName, fileMode, vertexMatchThreshold, trackMatchThreshold,
      writeTrackInfo);

  // CSV WRITERS
  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvParticleWriter, mex,
                             "CsvParticleWriter", inputParticles, outputDir,
                             outputStem, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvMeasurementWriter, mex,
                             "CsvMeasurementWriter", inputMeasurements,
                             inputClusters, inputMeasurementSimHitsMap,
                             outputDir, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvSimHitWriter, mex,
                             "CsvSimHitWriter", inputSimHits, outputDir,
                             outputStem, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvSpacePointWriter, mex,
                             "CsvSpacePointWriter", inputSpacepoints, outputDir,
                             outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvSpacePointsBucketWriter, mex,
                             "CsvSpacePointsBucketWriter", inputBuckets,
                             outputDir, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvTrackWriter, mex,
                             "CsvTrackWriter", inputTracks, outputDir, fileName,
                             inputMeasurementParticlesMap, outputPrecision,
                             nMeasurementsMin, truthMatchProbMin, ptMin);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvSeedWriter, mex, "CsvSeedWriter",
                             inputTrackParameters, inputSimSeeds, inputSimHits,
                             inputMeasurementParticlesMap,
                             inputMeasurementSimHitsMap, fileName, outputDir);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::CsvTrackingGeometryWriter, mex, "CsvTrackingGeometryWriter",
      trackingGeometry, outputDir, outputPrecision, writeSensitive,
      writeBoundary, writeSurfaceGrid, writeLayerVolume, writePerEvent);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::TrackFinderPerformanceWriter, mex,
      "TrackFinderPerformanceWriter", inputTracks, inputParticles,
      inputTrackParticleMatching, inputParticleTrackMatching, filePath,
      fileMode, effPlotToolConfig, fakeRatePlotToolConfig,
      duplicationPlotToolConfig, trackSummaryPlotToolConfig,
      subDetectorTrackSummaryVolumes, writeMatchingDetails);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::RootNuclearInteractionParametersWriter, mex,
      "RootNuclearInteractionParametersWriter", inputSimulationProcesses,
      filePath, fileMode, interactionProbabilityBins, momentumBins,
      invariantMassBins, multiplicityMax, writeOptionalHistograms,
      nSimulatedEvents);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvTrackParameterWriter, mex,
                             "CsvTrackParameterWriter", inputTrackParameters,
                             inputTracks, outputDir, outputStem,
                             outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvProtoTrackWriter, mex,
                             "CsvProtoTrackWriter", inputSpacepoints,
                             inputPrototracks, outputDir);

  {
    using Writer = ActsExamples::CsvBFieldWriter;

    auto w = py::class_<Writer>(mex, "CsvBFieldWriter");

    py::enum_<Writer::CoordinateType>(w, "CoordinateType")
        .value("rz", Writer::CoordinateType::RZ)
        .value("xyz", Writer::CoordinateType::XYZ);

    register_csv_bfield_writer_binding<Writer::CoordinateType::XYZ, true>(w);
    register_csv_bfield_writer_binding<Writer::CoordinateType::XYZ, false>(w);
    register_csv_bfield_writer_binding<Writer::CoordinateType::RZ, true>(w);
    register_csv_bfield_writer_binding<Writer::CoordinateType::RZ, false>(w);
  }

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::CsvExaTrkXGraphWriter, mex,
                             "CsvExaTrkXGraphWriter", inputGraph, outputDir,
                             outputStem);
}

}  // namespace Acts::Python
