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
ENABLE_TESTING()
INCLUDE(CTest)

# Include directories
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/adapters)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/VoxBoIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/PovRayIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/N3ITK)

IF(NOT MINC_TOOLKIT_BUILD)

  OPTION(USE_MINC4ITK "Use ITK MINC1/2 interface" ON)

  # Get ITK
  FIND_PACKAGE(ITK REQUIRED) # not compatible with ITK v 4 

  IF(USE_MINC4ITK)
    FIND_PACKAGE(MINC4ITK REQUIRED)
    FIND_PACKAGE(LIBMINC REQUIRED)
  ENDIF(USE_MINC4ITK)
  
  # Pull in FFTW
  FIND_LIBRARY(FFTW_LIB fftw3f)
  
ELSE(NOT MINC_TOOLKIT_BUILD)
  SET(USE_MINC4ITK ON)
ENDIF(NOT MINC_TOOLKIT_BUILD)

INCLUDE(${ITK_USE_FILE})

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
  adapters/BiasFieldCorrectionN4.cxx
  adapters/BinaryImageCentroid.cxx
  adapters/ClipImageIntensity.cxx
  adapters/ComputeFFT.cxx
  adapters/ComputeOverlaps.cxx
  adapters/ConnectedComponents.cxx
  adapters/CoordinateMap.cxx
  adapters/CopyTransform.cxx
  adapters/CreateImage.cxx
  adapters/CreateInterpolator.cxx
  adapters/ExtractRegion.cxx
  adapters/ExtractSlice.cxx
  adapters/FlipImage.cxx
  adapters/HistogramMatch.cxx
  adapters/ImageERF.cxx
  adapters/ImageLaplacian.cxx
  adapters/GeneralLinearModel.cxx
  adapters/LabelOverlapMeasures.cxx
  adapters/LabelStatistics.cxx
  adapters/LandmarksToSpheres.cxx
  adapters/LevelSetSegmentation.cxx
  adapters/MathematicalMorphology.cxx
  adapters/MixtureModel.cxx
#  adapters/MRFVote.cxx
  adapters/MultiplyImages.cxx
  adapters/NormalizedCrossCorrelation.cxx
  adapters/OverlayLabelImage.cxx
  adapters/PadImage.cxx
  adapters/PeronaMalik.cxx
  adapters/PrintImageInfo.cxx
  adapters/Rank.cxx
  adapters/ReadImage.cxx
  adapters/ReciprocalImage.cxx
  adapters/ReorderStack.cxx
  adapters/ReplaceIntensities.cxx
  adapters/ResampleImage.cxx
  adapters/ResliceImage.cxx
  adapters/SampleImage.cxx
  adapters/ScaleShiftImage.cxx
  adapters/SetSform.cxx
  adapters/SetOrientation.cxx
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
  adapters/WeightedSumVoxelwise.cxx
  adapters/WrapDimension.cxx
  adapters/WriteImage.cxx
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx)


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

IF(USE_MINC4ITK)
  LINK_DIRECTORIES(${LIBMINC_LIBRARY_DIRS} ${MINC4ITK_LIBRARY_DIRS})
  LINK_LIBRARIES(${MINC4ITK_LIBRARIES} ${EZMINC_LIBRARIES} ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES})

  IF(HAVE_MINC2)
   LINK_DIRECTORIES(${HDF5_LIBRARY_DIRS})
   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
   LINK_LIBRARIES(${HDF5_LIBRARY} ${ZLIB_LIBRARY})
  ENDIF(HAVE_MINC2)

  IF(HAVE_MINC1)
   LINK_DIRECTORIES(${NETCDF_LIBRARY_DIRS})
   INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIR})
   LINK_LIBRARIES(${NETCDF_LIBRARY})
  ENDIF(HAVE_MINC1)

  INCLUDE_DIRECTORIES(${LIBMINC_INCLUDE_DIRS} ${MINC4ITK_INCLUDE_DIRS})
  ADD_DEFINITIONS(-DHAVE_MINC4ITK)
ENDIF(USE_MINC4ITK)

ADD_EXECUTABLE(c3d Convert3DMain.cxx)

ADD_EXECUTABLE(c2d Convert2DMain.cxx)

INSTALL(TARGETS 
          c3d c2d
          RUNTIME DESTINATION bin )

INCLUDE(cmake/Package.cmake)

