cmake_minimum_required(VERSION 3.20 FATAL_ERROR) #minimum version of CMAKE required
project(UserActions)

# Load ActRoot custom script
set(ACTROOT $ENV{ACTROOT})
include(${ACTROOT}/Scripts.cmake)

#And call function
# First user action
add_userlibrary(NAME RecRANSAC SOURCES RecRANSAC.h RecRANSAC.cxx LINK ActAlgorithm)
add_userlibrary(NAME FilterDecay SOURCES FilterDecay.h FilterDecay.cxx LINK ActAlgorithm)
add_userlibrary(NAME CorrelateImplantDecay SOURCES CorrelateImplantDecay.h CorrelateImplantDecay.cxx LINK ActAlgorithm)
