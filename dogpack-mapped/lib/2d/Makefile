### definitions ###

# default/library definitions
COMPILE_WITH_UNST = 1
include $(DOGPACK)/lib/2d/Makefile.defs
ALL_OBJECTS = $(LIB2D_OBJECTS)
ALL_SOURCES = $(LIB2D_SOURCES)

### targets ###

lib: $(LIB2D_OBJECTS)

include $(DOGPACK)/lib/2d/Makefile.targets

cleanlib:
	rm -f *.o cart/*.o unst/*.o

### DO NOT remove this line - make depends on it ###
