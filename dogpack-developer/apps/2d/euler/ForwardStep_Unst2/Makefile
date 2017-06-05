
### definitions ###

# default/library definitions
COMPILE_WITH_UNST = 1
include $(DOGPACK)/lib/2d/Makefile.defs


# Default (application level) definintions
include $(DOGPACK)/apps/2d/euler/lib/Makefile.defs

# overridden library objects
#
#SetBndValues_Unst = SetBndValues_Unst
AfterQinit_Unst   = AfterQinit_Unst
ApplyShockLimiter_Unst = Apptest3
DFluxFunc         = DFluxFunc
D2FluxFunc         = D2FluxFunc
FluxFunc         = FluxFunc
SetWaveSpd       = SetWaveSpd
SetBndValues_Unst= SetBndValues_Unst2
#ConstructL_Unst  = ConstructL_Unst_glf2
#ConstructL_Unst  = ConstructL_Unst_glf
#ConstructL_Unst  = ConstructL_Unst_glf_good
ConstructL_Unst  = ConstructL_Unst_good
#ConstructL_Unst  = ConstructL_Unst
SetBndFlux       = SetBndFlux
LaxWendroff_Unst = LaxWendroff_Unst

# Eigen decomposition of the Flux function
ProjectLeftEig   = ProjectLeftEig
ProjectRightEig  = ProjectRightEig



## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
## ------------------------------------------------------------------------- ##

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
