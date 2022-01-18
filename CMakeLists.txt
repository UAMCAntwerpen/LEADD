cmake_minimum_required(VERSION 3.15)

# Define the project.
project(LEADD
  LANGUAGES CXX
  VERSION 2.0.0
  DESCRIPTION "Lamarckian Evolutionary Algorithm for de novo Drug Design"
  HOMEPAGE_URL "https://github.com/UAMCAntwerpen/LEADD"
)

# Create an option to allow the user to specify whether they want to build
# Python bindings or not. The default is ON.
option(BUILD_PYTHON_BINDINGS "Build LEADD Python bindings using pybind11" ON)

# Specify the C++ standard.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define the directories where the source code can be found and where the
# compiled libraries/executables should be stored.
set(LEADD_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LEADD_SOURCE_DIR ${LEADD_ROOT_DIR}/source)
set(LEADD_LIB_DIR ${LEADD_ROOT_DIR}/lib)
set(LEADD_BIN_DIR ${LEADD_ROOT_DIR}/bin)

# Include the LEADD library directory in the default RPATH of the targets.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${LEADD_LIB_DIR})

# Check that the RDKit was found. It's the responsibility of the user to specify its location.
if(NOT DEFINED RDKit_INCLUDE_DIRS OR NOT IS_DIRECTORY ${RDKit_INCLUDE_DIRS})
  message(FATAL_ERROR "Variable RDKit_INCLUDE_DIRS (${RDKit_INCLUDE_DIRS}) doesn't specify a valid directory.")
elseif(NOT DEFINED RDKit_LIBRARY_DIRS OR NOT IS_DIRECTORY ${RDKit_LIBRARY_DIRS})
  message(FATAL_ERROR "Variable RDKit_LIBRARY_DIRS (${RDKit_LIBRARY_DIRS}) doesn't specify a valid directory.")
else()
  set(RDKit_FOUND True)
endif()

# Search for the rest of LEADD's dependencies.
find_package(Boost REQUIRED COMPONENTS serialization program_options)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(SQLite3 REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Add the dependencies' directories to the include and library paths.
include_directories(
  ${RDKit_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${HDF5_CXX_INCLUDE_DIRS}
  ${SQLite3_INCLUDE_DIRS}
  ${OpenMP_CXX_INCLUDE_DIRS}
)

link_directories(
  ${RDKit_LIBRARY_DIRS}
  ${Boost_LIBRARY_DIRS}
)

message("-- Extra compilation flags: ${CMAKE_CXX_FLAGS}")

# Define all libraries included within the LEADD distribution and their dependencies.
add_library(Graph SHARED ${LEADD_SOURCE_DIR}/Graph.cpp)

add_library(Connection SHARED ${LEADD_SOURCE_DIR}/Connection.cpp)
target_link_libraries(Connection PUBLIC
  boost_serialization
)

add_library(Pseudofragment SHARED ${LEADD_SOURCE_DIR}/Pseudofragment.cpp)
target_link_libraries(Pseudofragment PUBLIC
  Connection
  RDKitRDGeneral
  RDKitGraphMol
  RDKitSmilesParse
  RDKitFingerprints
)

add_library(PseudofragmentDB SHARED ${LEADD_SOURCE_DIR}/PseudofragmentDB.cpp)
target_link_libraries(PseudofragmentDB PUBLIC
  Pseudofragment
  OpenMP::OpenMP_CXX
  ${SQLite3_LIBRARIES}
)

add_library(FragmentationSettings SHARED ${LEADD_SOURCE_DIR}/FragmentationSettings.cpp)

add_library(Fragmentation SHARED ${LEADD_SOURCE_DIR}/Fragmentation.cpp)
target_link_libraries(Fragmentation PUBLIC
  Pseudofragment
  FragmentationSettings
  RDKitChemTransforms
  RDKitForceFieldHelpers
)

add_library(ReconstructionSettings SHARED ${LEADD_SOURCE_DIR}/ReconstructionSettings.cpp)

add_library(ConnectionCompatibilities SHARED ${LEADD_SOURCE_DIR}/ConnectionCompatibilities.cpp)
target_link_libraries(ConnectionCompatibilities PUBLIC
  PseudofragmentDB
)

add_library(ConnectionQueryResults SHARED ${LEADD_SOURCE_DIR}/ConnectionQueryResults.cpp)
target_link_libraries(ConnectionQueryResults PUBLIC
  ConnectionCompatibilities
)

add_library(Reconstruction SHARED ${LEADD_SOURCE_DIR}/Reconstruction.cpp)
target_link_libraries(Reconstruction PUBLIC
  Graph
  ReconstructionSettings
  ConnectionQueryResults
  Fragmentation
  RDKitMolDraw2D
  RDKitDepictor
  ${HDF5_CXX_LIBRARIES}
)

add_library(SAScore SHARED ${LEADD_SOURCE_DIR}/SAScore.cpp)
target_link_libraries(SAScore PUBLIC
  RDKitGraphMol
  RDKitFingerprints
  RDKitDescriptors
  boost_serialization
)

add_library(LEADD SHARED ${LEADD_SOURCE_DIR}/LEADD.cpp)
target_link_libraries(LEADD PUBLIC
  Reconstruction
  SAScore
)

# Define the executables included wihin the LEADD distribution and their dependencies.
add_executable(DecodeConnection ${LEADD_SOURCE_DIR}/DecodeConnection.cpp)

add_executable(EncodeConnection ${LEADD_SOURCE_DIR}/EncodeConnection.cpp)

add_executable(Fragment ${LEADD_SOURCE_DIR}/Fragment.cpp)
target_link_libraries(Fragment
  PseudofragmentDB
  Fragmentation
  RDKitFileParsers
  boost_program_options
)

add_executable(MakeFragmentSimilarityMatrix ${LEADD_SOURCE_DIR}/MakeFragmentSimilarityMatrix.cpp)
target_link_libraries(MakeFragmentSimilarityMatrix
  PseudofragmentDB
  boost_program_options
  ${HDF5_CXX_LIBRARIES}
)

add_executable(PrecalculateConnectionQueryResults ${LEADD_SOURCE_DIR}/PrecalculateConnectionQueryResults.cpp)
target_link_libraries(PrecalculateConnectionQueryResults
  ReconstructionSettings
  ConnectionQueryResults
  boost_program_options
)

add_executable(InspectConnectionQueryResults ${LEADD_SOURCE_DIR}/InspectConnectionQueryResults.cpp)
target_link_libraries(InspectConnectionQueryResults
  ConnectionQueryResults
  boost_program_options
)

add_executable(BinAtomPairFrequencies ${LEADD_SOURCE_DIR}/BinAtomPairFrequencies.cpp)
target_link_libraries(BinAtomPairFrequencies
  ConnectionQueryResults
  boost_program_options
)

add_executable(MakeFeatureLibrary ${LEADD_SOURCE_DIR}/MakeFeatureLibrary.cpp)
target_link_libraries(MakeFeatureLibrary
  SAScore
  RDKitFileParsers
  boost_program_options
)

add_executable(CalculateSAScore ${LEADD_SOURCE_DIR}/CalculateSAScore.cpp)
target_link_libraries(CalculateSAScore
  SAScore
  RDKitSmilesParse
  RDKitFileParsers
  boost_program_options
)

add_executable(ConvertToReconstruction ${LEADD_SOURCE_DIR}/ConvertToReconstruction.cpp)
target_link_libraries(ConvertToReconstruction
  Reconstruction
  RDKitFileParsers
  boost_program_options
)

add_executable(StandaloneLEADD ${LEADD_SOURCE_DIR}/StandaloneLEADD.cpp)
target_link_libraries(StandaloneLEADD
  LEADD
  boost_program_options
)

# Specify where the libraries should be installed.
install(
  TARGETS
    Graph
    Connection
    Pseudofragment
    PseudofragmentDB
    FragmentationSettings
    Fragmentation
    ReconstructionSettings
    ConnectionCompatibilities
    ConnectionQueryResults
    Reconstruction
    SAScore
    LEADD
  DESTINATION ${LEADD_LIB_DIR}
)

# Specify where the executables should be installed.
install(
  TARGETS
    DecodeConnection
    EncodeConnection
    Fragment
    MakeFragmentSimilarityMatrix
    PrecalculateConnectionQueryResults
    InspectConnectionQueryResults
    BinAtomPairFrequencies
    MakeFeatureLibrary
    CalculateSAScore
    StandaloneLEADD
    ConvertToReconstruction
  DESTINATION ${LEADD_BIN_DIR}
)

# If requested, build the Python bindings.
if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(${LEADD_ROOT_DIR}/pybind11)
  find_package(Python REQUIRED COMPONENTS Interpreter Development)

  pybind11_add_module(pyLEADD ${LEADD_SOURCE_DIR}/pyLEADD.cpp)
  target_link_libraries(pyLEADD PUBLIC
    LEADD
  )

  install(
    TARGETS pyLEADD
    DESTINATION ${LEADD_LIB_DIR}
  )
endif()
