# ----------------------------------------------------------------
# INSTALLATION AND PACKAGING with CPack
# ----------------------------------------------------------------

# On Win32, we must include the redistributable
IF(MSVC80 OR MSVC90)
  FIND_PROGRAM(VCREDIST_X86 vcredist_x86.exe)
  IF(VCREDIST_X86)
    INSTALL(FILES ${VCREDIST_X86} DESTINATION bin)
    SET(CPACK_NSIS_EXTRA_INSTALL_COMMANDS 
      "ExecWait '\\\"$INSTDIR\\\\bin\\\\vcredist_x86.exe\\\" /q:a'"
      "CreateShortCut 'c:\\Windows\\system32\\cmd.exe' '$INSTDIR\\c3dshell.lnk'")
  ENDIF(VCREDIST_X86)
ENDIF(MSVC80 OR MSVC90)

# Allow package generation
SET(CPACK_PACKAGE_NAME "c3d")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "C3D Medical Image Processing Tool")
SET(CPACK_PACKAGE_VENDOR "itksnap.org")
SET(CPACK_PACKAGE_VERSION_MAJOR "${C3D_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${C3D_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${C3D_VERSION_PATCH}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "c3d-${C3D_VERSION_FULL}")
SET(CPACK_NSIS_MODIFY_PATH ON)


# Shamelessly stolen from ParaView_
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "c3d-${C3D_VERSION_FULL}")
IF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
  EXEC_PROGRAM(uname ARGS "-m" OUTPUT_VARIABLE CMAKE_SYSTEM_PROCESSOR)
ENDIF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
IF(NOT DEFINED CPACK_SYSTEM_NAME)
  SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
  IF(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win64-${CMAKE_SYSTEM_PROCESSOR})
  ELSE(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win32-${CMAKE_SYSTEM_PROCESSOR})
  ENDIF(CMAKE_CL_64)
ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)
IF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
  SET(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
ENDIF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)

# Show GPL license
SET(CPACK_RESOURCE_FILE_LICENSE "${CONVERT3D_SOURCE_DIR}/COPYING.txt")

IF(WIN32 AND NOT UNIX)

  SET(CPACK_GENERATOR "ZIP" "NSIS")

ELSE(WIN32 AND NOT UNIX)

  # Set the generator to either STGZ or Apple
  IF(NOT APPLE)
    SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "pauly2@mail.med.upenn.edu")
    SET(CPACK_GENERATOR "RPM" "TGZ" "STGZ" "DEB")
  ELSE(NOT APPLE)
    SET(CPACK_GENERATOR "PackageMaker")
  ENDIF(NOT APPLE)

ENDIF(WIN32 AND NOT UNIX)

INCLUDE(CPack)
