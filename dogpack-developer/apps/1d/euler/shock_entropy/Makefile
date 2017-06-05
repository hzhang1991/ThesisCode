### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(DOGPACK)/apps/1d/euler/lib/Makefile.defs


# overridden library objects
#
InitApp       =  InitApp
DogSolveUser  = $(DOGPACK)/lib/1d/DogSolveTDRK
ApplyLimiter  = ApplyPrimitivePointLimiter
#ApplyLimiter = ApplyPointLimiter_Backup

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) \
EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
EulerParams.cpp \

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
