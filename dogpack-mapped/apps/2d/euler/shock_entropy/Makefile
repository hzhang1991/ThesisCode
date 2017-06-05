
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# overridden library objects
#
AfterQinit  =  AfterQinit
gluInvertMatrix     = $(DOGPACK)/lib/2d/cart/gluInvertMatrix
IntegrateBasis = IntegrateBasisBackup
RiemannSolve = $(DOGPACK)/lib/2d/RiemannSolveLLF
ApplyLimiter = ApplyPrimitiveLimiter

ALL_OBJECTS = $(COMMON_OBJECTS) \
AppSolver.o \
mapc2p.o\
EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
AppSolver.cpp \
mapc2p.cpp\
EulerParams.cpp \
*.h \

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
