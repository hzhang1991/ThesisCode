
### definitions ###

# default/library definitions
COMPILE_WITH_UNST = 1
include $(DOGPACK)/scalar/lib/hybrid/Makefile.defs

include $(DOGPACK)/scalar/apps/hybrid/VlasovPoisson/lib/Makefile.defs

# overridden library objects
#
#ComputeElectricField = ComputeElectricField
SetBndValues_Unst     = SetBndValues_Unst

ALL_OBJECTS = $(COMMON_OBJECTS) $(VLASOV_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(VLASOV_LIB_SOURCES)

### targets ###

include $(DOGPACK)/scalar/lib/hybrid/Makefile.targets

### DO NOT remove this line - make depends on it ###
