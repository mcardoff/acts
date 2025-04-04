#!/usr/bin/env python3
import os
from collections import namedtuple

import acts
import acts.examples

from acts.examples import (
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    AdaptiveMultiVertexFinderAlgorithm,
    VertexNTupleWriter,
    RootVertexWriterNaive,
)

from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

# acts.examples.dump_args_calls(locals())  # show python binding calls

s = acts.examples.Sequencer(numThreads=1)
u = acts.UnitConstants
field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
global_logging_level = acts.logging.INFO

# beam_pos = acts.Vector3([50*u.um, 50*u.um, 1e2])
# beam_cov = acts.SquareMatrix3([[50*u.um, 0, 0],
                               # [0, 50*u.um, 0],
                               # [0, 0, 100.0]])

NTupleReaderConfig = namedtuple(
    "Config",
    [
        "inputTreeName",            # name of the input tree
        "inputFilePath",            # name of the input file
        "outputTrackParameters",    # track parameters
        "outputTruthVtxParameters", # truth vtx 4-position
        "outputTruthHSTrackInfo",
        "outputRecoVtxParameters",  # reco vtx 4-position (t always 0)
        "outputRecoHSTrackInfo",
    ],
    defaults = [None] * 7, 
)

ntuple_config = NTupleReaderConfig(
    # inputs
    inputTreeName="ntuple",
    inputFilePath="~/Projects/ntuple/user.scheong.42871997.Output._000031.SuperNtuple.root", # currently this is the VBF H -> inv sample
    # track parameters
    outputTrackParameters="trackParameters",
    # truth vertex info
    outputTruthVtxParameters="truthVertices",
    outputTruthHSTrackInfo="truthTracksHSInfo",
    # reco vertex info
    outputRecoVtxParameters="recoVertices",
    outputRecoHSTrackInfo="recoTracksHSInfo",
)

TrackSelectorConfig = namedtuple(
    "SelectorConfig",
    [
        
    ]
)

AMVFConfig = namedtuple(
    "Config",
    [
        "inputTrackParameters",
        "inputTruthParticles", # only used in truthSeeder
        "inputTruthVertices",  # only used in truthSeeder
        "outputProtoVertices",
        "outputVertices",
        "bField",
        "seedFinder",
        "maxIterations",
        "initialVariances",
        "useTime",
        "spatialBinExtent",   # only used in AdaptiveGridSeeder
        "temportalBinExtent", # only used in AdaptiveGridSeeder
    ],
    defaults = [None] * 12,
)

amvf_config = AMVFConfig(
    inputTrackParameters=ntuple_config.outputTrackParameters,
    outputProtoVertices="outputProtoVertices",
    outputVertices="outputVertices",
    bField=field,
    seedFinder=acts.VertexSeedFinder.GaussianSeeder,
    useTime=False,
)

RootVertexWriterConfig = namedtuple(
    "Config",
    [
        "inputVertices",         # ACTS vertices
        "inputTruthVertices",    # Athena Truth vertices
        "inputRecoVertices",     # Athena Reco vertices
        "inputTruthHSTrackInfo", # Athena Truth HS Vertex track info
        "inputRecoHSTrackInfo",  # Athena Reco HS Vertex track info
        "filePath",
        "fileMode",
        "treeName",
    ],
    defaults = [None] * 8,
)

writer_config = RootVertexWriterConfig(
    inputVertices=amvf_config.outputVertices,
    inputTruthVertices=ntuple_config.outputTruthVtxParameters,
    inputTruthHSTrackInfo=ntuple_config.outputTruthHSTrackInfo,
    inputRecoVertices=ntuple_config.outputRecoVtxParameters,
    inputRecoHSTrackInfo=ntuple_config.outputRecoHSTrackInfo,
    filePath='./amvf_test.root', # be sure to include extension
    treeName='tree',
)

# inputs to RootAthenaNTupleReader: config, logginglevel
# NTupleReader = acts.examples.RootAthenaNTupleReader(
NTupleReader = acts.examples.RootAthenaNTupleReaderNaive(
    level=global_logging_level,
    **acts.examples.defaultKWArgs(
        # inputs
        inputTreeName=ntuple_config.inputTreeName,
        inputFilePath=ntuple_config.inputFilePath,
        # tracks 
        outputTrackParameters=ntuple_config.outputTrackParameters,
        # truth
        outputTruthVtxParameters=ntuple_config.outputTruthVtxParameters,
        outputTruthHSTrackInfo=ntuple_config.outputTruthHSTrackInfo,
        # reco
        outputRecoVtxParameters=ntuple_config.outputRecoVtxParameters,
        outputRecoHSTrackInfo=ntuple_config.outputRecoHSTrackInfo,
    )
)

# add root reader
s.addReader(NTupleReader)

findVertices = AdaptiveMultiVertexFinderAlgorithm(
    level=global_logging_level,
    inputTrackParameters=amvf_config.inputTrackParameters,
    outputProtoVertices=amvf_config.outputProtoVertices,
    outputVertices=amvf_config.outputVertices,
    bField=amvf_config.bField,
    seedFinder=amvf_config.seedFinder,
    **acts.examples.defaultKWArgs(
        initialVariances=amvf_config.initialVariances,
        useTime=amvf_config.useTime,
    ),
)

# s.addAlgorithm(findVertices)

rootVertexWriter = RootVertexWriterNaive(
    level=global_logging_level,
    **acts.examples.defaultKWArgs(
        inputVertices         = writer_config.inputVertices,
        inputTruthVertices    = writer_config.inputTruthVertices,
        inputRecoVertices     = writer_config.inputRecoVertices,
        inputTruthHSTrackInfo = writer_config.inputTruthHSTrackInfo,
        inputRecoHSTrackInfo  = writer_config.inputRecoHSTrackInfo,
        filePath              = writer_config.filePath,
        fileMode              = writer_config.fileMode,
        treeName              = writer_config.treeName,
    )
)


# add vertexwriter
# s.addWriter(rootVertexWriter)

s.run()
