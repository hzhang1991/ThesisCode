#
#               Makefile for the lib code:
#
#       To make an executable, type:            make
#       (from the application sub-directory)

### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/asynch/Makefile.defs

# compiler options
#
#FFLAGS =  -c -O4
#LFLAGS =  
# -show: display g++ commands that would be executed
# -echo: show *all* commands and execute them
#LINK    = /usr/local/hdf5/bin/h5c++ -show # -E

TAGSOURCES = $(ALL_SOURCES)

### targets ###

.cpp.o: ; $(LINK) $(FFLAGS) $*.cpp -o $*.o

lib: $(LIBOBJECTS_COMMON)

$(Output).o: $(Output).cpp
	$(COMPILE) $(SET_COMPILE_WITH_HDF5) \
	$(Output).cpp -o $(Output).o

program : $(ALL_SOURCES)
	cat $(ALL_SOURCES) > dog1program.cpp

clean:
	rm -f *.o */*.o

retags:
	ctags $(TAGSOURCES)
	$(DOGPACK)/scripts/makefiletags $(TAGSOURCES) >> tags
	sort tags -o tags

tags: retags

### DO NOT remove this line - make depends on it ###
