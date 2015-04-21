#-----------------------------------------------------------------------------
# Name of the project
PROJECT(LogDomainDemonsRegistration)
find_package(ITK REQUIRED)

#-----------------------------------------------------------------------------
# Requires a recent version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


#-----------------------------------------------------------------------------
# Include some useful macros
INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)
INCLUDE(${CMAKE_SOURCE_DIR}/SDCMakeMacros.txt)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

#find_package(ITK REQUIRED)
#find_package(Matlab)

#-----------------------------------------------------------------------------
#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists


SET(Required_Packages
ITK
Matlab
)


#-----------------------------------------------------------------------------
#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

include(${ITK_USE_FILE})
#-----------------------------------------------------------------------------
#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc

SET(Libraries
ITKCommon
 #ITKIO
ITKStatistics
)


#-----------------------------------------------------------------------------
ADD_MEX_FILE(velocityfieldexp mex_velocityfieldexp.cpp)
#`TARGET_LINK_LIBRARIES(velocityfieldexp ${Libraries}  ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(velocityfieldexp  ${ITK_LIBRARIES})


ADD_MEX_FILE(invcondemonsforces mex_invcondemonsforces.cpp)
TARGET_LINK_LIBRARIES(invcondemonsforces   ${ITK_LIBRARIES})

ADD_MEX_FILE(weightedfwdemonsforces mex_weightedfwdemonsforces.cpp)
TARGET_LINK_LIBRARIES(weightedfwdemonsforces   ${ITK_LIBRARIES})

ADD_MEX_FILE(warpimage mex_warpimage.cpp)
TARGET_LINK_LIBRARIES(warpimage  ${ITK_LIBRARIES})

ADD_MEX_FILE(warplabelimage mex_warplabelimage.cpp)
TARGET_LINK_LIBRARIES(warplabelimage  ${ITK_LIBRARIES})

ADD_MEX_FILE(deffieldharmonicenergy mex_deffieldharmonicenergy.cpp)
TARGET_LINK_LIBRARIES(deffieldharmonicenergy  ${ITK_LIBRARIES})

ADD_MEX_FILE(deffieldjacobiandeterminant mex_deffieldjacobiandeterminant.cpp)
TARGET_LINK_LIBRARIES(deffieldjacobiandeterminant  ${ITK_LIBRARIES})

ADD_MEX_FILE(deffieldjacobiandist mex_deffieldjacobiandist.cpp)
TARGET_LINK_LIBRARIES(deffieldjacobiandist  ${ITK_LIBRARIES})

ADD_MEX_FILE(writewarpfile mex_writewarpfile.cpp)
TARGET_LINK_LIBRARIES(writewarpfile  ${Libraries} ${ITK_LIBRARIES})

ADD_MEX_FILE(readwarpfile mex_readwarpfile.cpp)
TARGET_LINK_LIBRARIES(readwarpfile  ${Libraries} ${ITK_LIBRARIES})
#ADD_MEX_FILE(writewarpfile mex_writewarpfile.cpp)
#TARGET_LINK_LIBRARIES(writewarpfile  ${Libraries} ${ITK_LIBRARIES})

#ADD_MEX_FILE(readwarpfile mex_readwarpfile.cpp)
#TARGET_LINK_LIBRARIES(readwarpfile  ${Libraries} ${ITK_LIBRARIES})




