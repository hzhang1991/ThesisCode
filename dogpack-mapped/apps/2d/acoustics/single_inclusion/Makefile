
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# compiler options -- see file $(DOGPACK)/config/Makefile.config

# overridden library objects
IntegrateBasis     = $(DOGPACK)/lib/2d/cart/IntegrateBasis
gluInvertMatrix     = $(DOGPACK)/lib/2d/cart/gluInvertMatrix
RiemannSolve  = RiemannSolve
IntegrateBasis = IntegrateBasisBackup
# place to add object modules
#
ALL_OBJECTS = \
$(COMMON_OBJECTS) \
$(SL_LIB_OBJECTS) \
AcousticParams.o \
AppSolver.o \
Parameters.o \
mapc2p.o\
#$(IntegrateBasis).o\
$(gluInvertMatrix).o\

ALL_SOURCES = \
$(COMMON_SOURCES) \
$(SL_LIB_SOURCES) \
AcousticParams.cpp \
AppSolver.cpp \
Parameters.cpp \
mapc2p.cpp\
#$(IntegrateBasis).cpp\
$(gluInvertMatrix).cpp\

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
