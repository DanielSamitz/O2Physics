# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

add_subdirectory(PID)

o2physics_add_dpl_workflow(trackextension
                    SOURCES trackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                                          O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2::DetectorsBase
                                          O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(trackselection
                    SOURCES trackselection.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-selection
                    SOURCES eventSelection.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase O2::CCDB O2::CommonConstants
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multiplicity-table
                    SOURCES multiplicityTable.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(centrality-table
                    SOURCES centralityTable.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2::CCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(timestamp
                    SOURCES timestamp.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::DetectorsRaw O2Physics::AnalysisCore O2::CommonDataFormat O2::CCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(weak-decay-indices
                    SOURCES weakDecayIndices.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ft0-corrected-table
                    SOURCES ft0CorrectedTable.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore  O2Physics::AnalysisCCDB O2::DetectorsBase O2::CCDB O2::CommonConstants
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-propagation
                    SOURCES trackPropagation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-converter
                    SOURCES mcConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fdd-converter
                    SOURCES fddConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(calo-clusters
                    SOURCES caloClusterProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DataFormatsPHOS O2::PHOSBase O2::PHOSReconstruction
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fwdtrackextension
                    SOURCES fwdtrackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                                          O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2::DetectorsBase
                                          O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)
