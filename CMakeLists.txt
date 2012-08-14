cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

PROJECT(PolyaffineLogDomainDemonsRegistration)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

SET(PolyaffineLogDomainDemonsRegistrationHeaders
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkESMDemonsRegistrationMaskFunction.h
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkESMDemonsRegistrationMaskFunction.txx
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineLogDomainDeformableRegistrationFilter.h
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineLogDomainDeformableRegistrationFilter.txx
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineLogDomainDemonsRegistrationFilter.h
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineLogDomainDemonsRegistrationFilter.txx
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineSymmetricLogDomainDemonsRegistrationFilter.h
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineSymmetricLogDomainDemonsRegistrationFilter.txx
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineStationaryVelocityFieldTransform.h
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkPolyaffineStationaryVelocityFieldTransform.txx
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkManualPolyaffineStationaryVelocityFieldTransform.h
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/itkManualPolyaffineStationaryVelocityFieldTransform.txx  
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/WeightedSVDTree.h
)

SET(PolyaffineLogDomainDemonsRegistrationSrcs
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/Main.cxx
  ${PolyaffineLogDomainDemonsRegistration_SOURCE_DIR}/WeightedSVDTree.cxx
)

ADD_EXECUTABLE( PolyaffineLogDomainDemonsRegistration ${PolyaffineLogDomainDemonsRegistrationHeaders} ${PolyaffineLogDomainDemonsRegistrationSrcs})

TARGET_LINK_LIBRARIES( PolyaffineLogDomainDemonsRegistration
  ${ITK_LIBRARIES}
)
