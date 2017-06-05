### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(DOGPACK)/apps/2d/euler/lib/Makefile.defs

# overridden library objects
#
AfterQinit          = AfterQinit
RiemannSolve        = $(DOGPACK)/lib/2d/RiemannSolveLLF
DogSolveUser        = $(DOGPACK)/lib/2d/cart/DogSolveTD
#ApplyLimiter        = Apptest3
ApplyLimiter        = $(DOGPACK)/lib/2d/cart/ApplyCharBarLimiter2
ApplyPosLimiter     = ApplyPosLimiter
ConstructL          = ConstructLGood
SetBndFlux          = SetBndFlux
SetBndValues        = SetBndValues2



ALL_OBJECTS = $(COMMON_OBJECTS) \
AppSolver.o \
EulerParams.o \


ALL_SOURCES = $(COMMON_SOURCES) \
AppSolver.cpp \
EulerParams.cpp \
*.h \

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
