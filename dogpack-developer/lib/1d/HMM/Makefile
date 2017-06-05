#
#               Makefile for the dogpack code:
#
#
#
FFLAGS =  -c -O4
LFLAGS =
LINK    = g++
.cpp.o: ; $(LINK) $(FFLAGS) $*.cpp -o $*.o

HMMOBJECTS = \
  $(DOGPACK)/1d/lib/HMM/RunDogpackHMM.o \
  $(DOGPACK)/1d/lib/HMM/GetParamsHMM.o \
  $(DOGPACK)/1d/lib/HMM/SnapToMaster.o \
  $(DOGPACK)/1d/lib/HMM/GridPatch.o \
  $(DOGPACK)/1d/lib/HMM/BeforeStep_sub.o \
  $(DOGPACK)/1d/lib/HMM/AfterStep_sub.o \
  $(DOGPACK)/1d/lib/HMM/DogSolveRK_HMM.o \
  $(DOGPACK)/1d/lib/HMM/AdvanceOneStage_RK.o \
  $(DOGPACK)/1d/lib/HMM/SetBC_SubGrid.o \
  $(DOGPACK)/1d/lib/HMM/SetValues_MainGrid.o \
  $(DOGPACK)/1d/lib/HMM/ConstructL_HMM.o \
  $(DOGPACK)/1d/lib/HMM/ConstructL_HMM_short.o \
  $(DOGPACK)/1d/lib/HMM/UpdateSoln_HMM.o \
  $(DOGPACK)/1d/lib/HMM/GetHermiteCoeffs.o

SOURCES = \
  $(DOGPACK)/1d/lib/HMM/RunDogpackHMM.cpp \
  $(DOGPACK)/1d/lib/HMM/GetParamsHMM.cpp \
  $(DOGPACK)/1d/lib/HMM/SnapToMaster.cpp \
  $(DOGPACK)/1d/lib/HMM/GridPatch.cpp \
  $(DOGPACK)/1d/lib/HMM/BeforeStep_sub.cpp \
  $(DOGPACK)/1d/lib/HMM/AfterStep_sub.cpp \
  $(DOGPACK)/1d/lib/HMM/DogSolveRK_HMM.cpp \
  $(DOGPACK)/1d/lib/HMM/AdvanceOneStage_RK.cpp \
  $(DOGPACK)/1d/lib/HMM/SetBC_SubGrid.cpp \
  $(DOGPACK)/1d/lib/HMM/SetValues_MainGrid.cpp \
  $(DOGPACK)/1d/lib/HMM/ConstructL_HMM.cpp \
  $(DOGPACK)/1d/lib/HMM/ConstructL_HMM_short.cpp \
  $(DOGPACK)/1d/lib/HMM/UpdateSoln_HMM.cpp \
  $(DOGPACK)/1d/lib/HMM/GetHermiteCoeffs.cpp

hmm: $(HMMOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > dogpack1hmm.`date +%y%m%d`

### DO NOT remove this line - make depends on it ###
