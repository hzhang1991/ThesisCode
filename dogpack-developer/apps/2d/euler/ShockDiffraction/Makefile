
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(DOGPACK)/apps/2d/euler/lib/Makefile.defs


# overridden library objects
#
# For this problem, we need to override how the ConstructL's are called in
# order to set the correct boundary conditions.  The problem is how the corner
# element is handled, because it is assumed that when calling ConstructL, that
# a single ghost cell is capable of handling all of the data for interior
# points.  In these overridden routines, we introduce a SetBndValuesX and
# SetBndValuesY in order to get the correct reflective boundary conditions for
# the corner element.
AfterQinit          = AfterQinit
LaxWendroff         = LaxWendroff
ConstructL_LLF      = ConstructL_LLF
ApplyPosMPPLimiter  = ApplyPosMPPLimiter
ConstructL          = ConstructL_globalLF
AfterFullTimeStep   = AfterFullTimeStep
ApplyLimiter        = ApptestBest

ALL_OBJECTS = $(COMMON_OBJECTS) \
AppSolver.o \
EulerParams.o \
RiemannSolveGLF.o \

ALL_SOURCES = $(COMMON_SOURCES) \
AppSolver.cpp \
EulerParams.cpp \
RiemannSolveGLF.cpp \
*.h \

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
