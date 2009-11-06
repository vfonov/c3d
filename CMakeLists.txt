PROJECT(CONVERT3D)

# Nonsense for cmake 2.6 compatibility
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Versioning information
SET(C3D_VERSION_MAJOR 0)
SET(C3D_VERSION_MINOR 8)
SET(C3D_VERSION_PATCH 2)
SET(C3D_VERSION_FULL "${C3D_VERSION_MAJOR}.${C3D_VERSION_MINOR}.${C3D_VERSION_PATCH}")

# Include directories
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/adapters)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/VoxBoIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/PovRayIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/N3ITK)

# Get ITK
FIND_PACKAGE(ITK 3.14 REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Get Boost (for regex and whatever else)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_DEBUG TRUE)
FIND_PACKAGE(Boost 1.33 COMPONENTS regex)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
  LINK_LIBRARIES(${Boost_LIBRARIES})
ELSE(Boost_FOUND)
  MESSAGE(FATAL_ERROR
    "Could not find Boost libraries. Try adding a BOOST_ROOT entry to the cache")
ENDIF(Boost_FOUND)

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
  SOURCE_GROUP("Adapter Sources" REGULAR_EXPRESSION "adapters/*cxx")
  SOURCE_GROUP("Adapter Headers" REGULAR_EXPRESSION "adapters/*h")
ENDIF(WIN32)

SET(SOURCES
  adapters/AddImages.cxx
  adapters/AntiAliasImage.cxx
  adapters/ApplyMetric.cxx
  adapters/BiasFieldCorrection.cxx
  adapters/BinaryImageCentroid.cxx
  adapters/ClipImageIntensity.cxx
  adapters/ComputeFFT.cxx
  adapters/ComputeOverlaps.cxx
  adapters/ConnectedComponents.cxx
  adapters/CopyTransform.cxx
  adapters/CreateImage.cxx
  adapters/CreateInterpolator.cxx
  adapters/ExtractRegion.cxx
  adapters/HistogramMatch.cxx
  adapters/ImageERF.cxx
  adapters/ImageLaplacian.cxx
  adapters/GeneralLinearModel.cxx
  adapters/LabelStatistics.cxx
  adapters/LevelSetSegmentation.cxx
  adapters/MathematicalMorphology.cxx
  adapters/MixtureModel.cxx
  adapters/MultiplyImages.cxx
  adapters/NormalizedCrossCorrelation.cxx
  adapters/PadImage.cxx
  adapters/PeronaMalik.cxx
  adapters/PrintImageInfo.cxx
  adapters/Rank.cxx
  adapters/ReadImage.cxx
  adapters/ReciprocalImage.cxx
  adapters/ReplaceIntensities.cxx
  adapters/ResampleImage.cxx
  adapters/ResliceImage.cxx
  adapters/SampleImage.cxx
  adapters/ScaleShiftImage.cxx
  adapters/SignedDistanceTransform.cxx
  adapters/SmoothImage.cxx
  adapters/SplitMultilabelImage.cxx
  adapters/StapleAlgorithm.cxx
  adapters/ThresholdImage.cxx
  adapters/TrimImage.cxx
  adapters/UnaryMathOperation.cxx
  adapters/UpdateMetadataKey.cxx
  adapters/Vote.cxx
  adapters/VoxelwiseRegression.cxx
  adapters/WarpImage.cxx
  adapters/WarpLabelImage.cxx
  adapters/WeightedSum.cxx
  adapters/WriteImage.cxx
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx)

# Pull in FFTW
FIND_LIBRARY(FFTW_LIB fftw3f)

# Define the ITK libraries
SET(ITK_LIBRARIES
  ITKNumerics ITKIO ITKCommon ITKBasicFilters ITKStatistics ITKAlgorithms)

# Configure the version file
CONFIGURE_FILE(
  ${CONVERT3D_SOURCE_DIR}/ConvertImageVersion.cxx.in
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx @ONLY IMMEDIATE)

# Get the extra stuff compiled
SUBDIRS(itkextras)
SUBDIRS(utilities)

ADD_LIBRARY(cnd_adapters ${SOURCES})
ADD_LIBRARY(cnd_driver ConvertImageND.cxx)

LINK_LIBRARIES(
  cnd_adapters cnd_driver ${ITK_LIBRARIES} ITKVoxBoIO ITKPovRayIO)

ADD_EXECUTABLE(c3d Convert3DMain.cxx)

ADD_EXECUTABLE(c2d Convert2DMain.cxx)

INSTALL_TARGETS("/bin" c3d c2d)

INCLUDE(cmake/Package.cmake)

# CDash Support
ENABLE_TESTING()
INCLUDE(CTest)
