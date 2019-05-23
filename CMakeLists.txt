macro (POST_COPY SRC_TARGET SRC_PATH SRC_MATCH DEST_DIR)
  # Make sure dir exists
  MAKE_DIRECTORY ( "${DEST_DIR}" )  
  # Locate all files matching pattern
  message ( STATUS "Packaging: ${SRC_PATH}/${SRC_MATCH} -> ${DEST_DIR}" )    
  if ( "${SRC_MATCH}" MATCHES "(.*)\\*(.*)" )  
     # wildcard given, find all matching files (MUST EXIST NOW, AT CMAKE TIME, BEFORE BUILD)
     FILE ( GLOB src_list ABSOLUTE "${SRC_PATH}/${SRC_MATCH}" )
  else ()
     # specific file, just copy it
     SET ( src_list "${SRC_PATH}/${SRC_MATCH}")
  endif ()  
  # Add post-build copy
  foreach ( src_file ${src_list})
    ADD_CUSTOM_COMMAND ( TARGET ${SRC_TARGET}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${src_file} ${DEST_DIR}
    )
  endforeach () 
endmacro ()

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Create solution (cmake calls it a 'project')
PROJECT ( ProteinLevelSetSurface )
GET_FILENAME_COMPONENT ( BUILD_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../build/" REALPATH )
GET_FILENAME_COMPONENT ( SOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../" REALPATH )
GET_FILENAME_COMPONENT ( CMAKE_PATH  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake_modules/" REALPATH )
SET ( BASE_BUILD_PATH  ${BUILD_PATH} CACHE PATH "Base build output path.." )
SET ( BASE_SOURCE_PATH  ${SOURCESPATH} CACHE PATH "Base source code path..")
SET ( CMAKE_MODULE_PATH ${CMAKE_PATH} CACHE PATH "Path to Find.cmake modules" )

# Better search hints for Boost 
#  FILE ( GLOB children RELATIVE ${BASE_BUILD_PATH} ${BASE_BUILD_PATH}/*)
#  UNSET ( BOOST_ROOT )
#  FOREACH(subdir ${children})
#      IF ( "${subdir}" MATCHES "boost(.*)" )
#        SET ( BOOST_ROOT ${BASE_BUILD_PATH}/${subdir}/ CACHE PATH "Path to Boost")
#        SET ( BOOST_LIBRARYDIR "${BOOST_ROOT}/lib64-msvc-14.0" )
#     ENDIF()
#  ENDFOREACH()    


# Find required libraries
FIND_PACKAGE ( OpenVDB REQUIRED )
FIND_PACKAGE ( OpenEXR REQUIRED )  
FIND_PACKAGE ( Boost REQUIRED COMPONENTS iostreams system thread )
FIND_PACKAGE ( ILMBase REQUIRED )
FIND_PACKAGE ( TBB REQUIRED )
FIND_PACKAGE ( ZLIB REQUIRED )
FIND_PACKAGE ( OpenGL REQUIRED )
FIND_PACKAGE ( GLFW3 REQUIRED)
IF (WIN32)
  FIND_PACKAGE ( GLEW REQUIRED )
ENDIF()
FIND_PACKAGE ( CGAL REQUIRED)

INCLUDE(${CGAL_USE_FILE})

IF ( GLFW_FOUND )
  SET ( OPENVDB_USE_GLFW_FLAG "-DOPENVDB_USE_GLFW_3" CACHE STRING "Use GLFW3 for OpenVDB")
ENDIF ()

IF (APPLE)
  FIND_LIBRARY ( COCOA_LIBRARY Cocoa )
  FIND_LIBRARY ( IOKIT_LIBRARY IOKit )
  FIND_LIBRARY ( COREVIDEO_LIBRARY CoreVideo )
ENDIF ()
IF ( UNIX AND NOT APPLE )
  MESSAGE ( "UNIX" )
  SET ( XI_LIBRARY Xi )
  SET ( XRANDR_LIBRARY Xrandr )
ENDIF ()
IF (WIN32)
  # Because of implicit linking!  
  LINK_DIRECTORIES ( ${Boost_LIBRARY_DIR_RELEASE} )
ENDIF ()

# Project source files
SET ( ProteinLevelSetSurface_SOURCE_FILES
  main.cpp
  Density.h
  Density.cpp
 )
# Create project 
ADD_EXECUTABLE ( ProteinLevelSetSurface
  ${ProteinLevelSetSurface_SOURCE_FILES}
  )

SET_SOURCE_FILES_PROPERTIES ( ${ProteinLevelSetSurface_SOURCE_FILES}
  PROPERTIES
  COMPILE_FLAGS "${OPENVDB_USE_GLFW_FLAG} -DGL_GLEXT_PROTOTYPES=1 /bigobj"
)

# Project libraries (will become Project -> Linker -> Input)
TARGET_LINK_LIBRARIES ( ProteinLevelSetSurface debug
  ${OPENVDB_LIBRARY}
  ${OPENEXR_ILMIMF_LIBRARY_DEBUG}
  ${ILMBASE_HALF_LIBRARY_DEBUG}
  ${Boost_SYSTEM_LIBRARY_DEBUG}
  ${Boost_THREAD_LIBRARY_DEBUG}
  ${TBB_LIBRARY_DEBUG}  
  ${ZLIB_LIBRARY}
  ${GLFW_LIBRARY}
  ${GLEW_LIBRARY}
  ${CGAL_LIBRARIES}
  ${CGAL_3RD_PARTY_LIBRARIES}
  )
TARGET_LINK_LIBRARIES ( ProteinLevelSetSurface optimized
  ${OPENVDB_LIBRARY}
  ${OPENGL_gl_LIBRARY}
  ${OPENGL_glu_LIBRARY}
  ${COCOA_LIBRARY}
  ${IOKIT_LIBRARY}
  ${COREVIDEO_LIBRARY}
  ${OPENEXR_ILMIMF_LIBRARY_RELEASE}
  ${ILMBASE_HALF_LIBRARY_RELEASE}
  ${Boost_SYSTEM_LIBRARY_RELEASE}
  ${Boost_THREAD_LIBRARY_RELEASE}
  ${TBB_LIBRARY}    
  ${ZLIB_LIBRARY}
  ${GLFW_LIBRARY}
  ${GLEW_LIBRARY}
  ${XI_LIBRARY}
  ${XRANDR_LIBRARY}
  ${CGAL_LIBRARIES}
  ${CGAL_3RD_PARTY_LIBRARIES}
  )

# Project includes (will become Project -> Additional Includes)
INCLUDE_DIRECTORIES ( ${OPENVDB_INCLUDE_DIR} )
INCLUDE_DIRECTORIES ( ${OPENEXR_INCLUDE_DIR} )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIR} )
INCLUDE_DIRECTORIES ( ${ILMBASE_INCLUDE_DIR} )
INCLUDE_DIRECTORIES ( ${TBB_INCLUDE_DIR} )
INCLUDE_DIRECTORIES ( ${OPENGL_INCLUDE_DIR} )

INCLUDE_DIRECTORIES ( ${ZLIB_INCLUDE_DIR} )
INCLUDE_DIRECTORIES ( ${GLFW_INCLUDE_DIRS} )

IF (WIN32)
  INCLUDE_DIRECTORIES ( ${GLEW_INCLUDE_PATH} )
ENDIF ()

INCLUDE_DIRECTORIES ( .. ) # So they can simulate referencing installed headers
INCLUDE_DIRECTORIES ( . ) # So they can simulate referencing installed headers

IF (WIN32)
  ADD_DEFINITIONS ( -D_ITERATOR_DEBUG_LEVEL=0 -D_WIN32 -DNOMINMAX -DOPENVDB_STATICLIB -DOPENEXR_DLL )
  ADD_DEFINITIONS ( ${OPENVDB_USE_GLFW_FLAG} )
ELSE ()
  ADD_DEFINITIONS ( -pthread -fPIC )
ENDIF ()

# Copy required dlls to output path
# (Will go to Project -> Build Events -> Post Build )
POST_COPY ( ProteinLevelSetSurface "${GLEW_ROOT_DIR}/bin/Release/x64" glew32.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${OPENEXR_HOME}/lib IlmImf.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${ILMBASE_HOME}/lib Half.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${ILMBASE_HOME}/lib Iex.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${ILMBASE_HOME}/lib Imath.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${ILMBASE_HOME}/lib IlmThread.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${ZLIB_ROOT_DIR}/lib zlib.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${TBB_BINARY_DIR} *.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${Boost_LIBRARY_DIR_RELEASE} boost_system*.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${Boost_LIBRARY_DIR_RELEASE} boost_thread*.dll "${CMAKE_BINARY_DIR}/Release" )
POST_COPY ( ProteinLevelSetSurface ${Boost_LIBRARY_DIR_RELEASE} boost_chrono*.dll "${CMAKE_BINARY_DIR}/Release" )
