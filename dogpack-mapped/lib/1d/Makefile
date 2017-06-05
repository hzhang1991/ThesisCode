### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs
ALL_OBJECTS = $(LIB1D_OBJECTS)
ALL_SOURCES = $(LIB1D_SOURCES)

### targets ###

lib: $(LIB1D_OBJECTS)

$(Output).o: $(Output).cpp
	$(COMPILE) $(SET_COMPILE_WITH_HDF5) \
	$(Output).cpp -o $(Output).o

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
