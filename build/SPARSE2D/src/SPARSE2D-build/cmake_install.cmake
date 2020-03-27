# Install script for directory: /local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/extern")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sparse2d" TYPE FILE FILES
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Array.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Atrou3D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Border.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/CErf.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/CMem.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/DefFunc.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/DefMath.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/FCur.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/FFTN.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/FFTN_1D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/FFTN_2D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Filter.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Fista.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/FloatTrans.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/GMCA.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/GenMedian.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/GetLongOptions.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/GlobalInc.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM1D_ALDCT.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM1D_Block.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM1D_Dct.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM1D_IO.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM3D_IO.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_BIO.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Block2D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Color.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_DCT.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Deconv.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Edge.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Graphics.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_IO.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_IOCol.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_IOTools.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Lut.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Math.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Morpho.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Noise.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Obj.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Prob.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Regul.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Rot.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_Sigma.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/IM_VisTool.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Licence.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/LineCol.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MGA_Inc.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D1D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D_Deconv.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D_Filter.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D_NoiseModel.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D_Obj.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D_Regul.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D_Segment.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR1D_Sigma.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR2D1D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR3D_Obj.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Abaque.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_CorrNoise.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Deconv.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Filter.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_HaarPoisson.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_HaarPoissonFilter.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_ListObj.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_MVM.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Noise.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_NoiseModel.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Obj.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Psupport.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Rayleigh.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Sigma.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_SoftRegul.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Support.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_Threshold.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_VisElem.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MR_VisTree.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MW1D_Filter.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MW_Deconv.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MW_Filter.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MatrixOper.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Memory.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MeyerWT.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/MeyerWT1D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Mr_FE.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Mr_FewEvent.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Mr_FewEvent1d.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Mr_FewEvent2d.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/NR.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/NR_util.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Nesterov.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/OptMedian.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/SB_Filter.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/SB_Filter1D.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/SB_Filter_float.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/SoftInfo.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/TempArray.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/TempMemory.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/Usage.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/VMS.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/WPackets.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/WT2D_CF.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/WT_Bord.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/WT_Mirror.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/fractal.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/io.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/longnam.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/macro.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/mr_FiltUsage.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/mr_com.h"
    "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D/src/libsparse2d/writefits3d.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/libsparse2d.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im1d_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im1d_deconv")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im1d_deconv"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/im1d_deconv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im1d_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im1d_deconv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im1d_deconv")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_deconv")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_deconv"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/im_deconv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_deconv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_deconv")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_mirror" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_mirror")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_mirror"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/im_mirror")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_mirror" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_mirror")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im_mirror")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_filter")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_filter"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr1d_filter")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_filter")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_filter")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_gmca")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_gmca"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr1d_gmca")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_gmca")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_gmca")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_recons")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_recons"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr1d_recons")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_recons")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_recons")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_trans" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_trans")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_trans"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr1d_trans")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_trans" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_trans")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr1d_trans")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr2d1d_trans" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr2d1d_trans")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr2d1d_trans"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr2d1d_trans")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr2d1d_trans" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr2d1d_trans")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr2d1d_trans")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr_deconv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr_filter")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr_gmca")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_info" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_info")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_info"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr_info")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_info" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_info")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_info")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr_recons")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_sigma" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_sigma")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_sigma"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr_sigma")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_sigma" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_sigma")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_sigma")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mr_transform")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/mw_deconv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/local/home/themelis/Development/Cpp/glimpse-wiener-beta/build/SPARSE2D/src/SPARSE2D-build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
