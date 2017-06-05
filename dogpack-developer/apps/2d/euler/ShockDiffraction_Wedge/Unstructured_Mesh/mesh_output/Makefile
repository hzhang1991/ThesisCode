#
#       MeshGenC++ Makefile
#       To make an executable, type: make
#       (from the application sub-directory)

### definitions ###

# default/library definitions
include $(MESHGENCPP)/lib/2d/Makefile.defs

# compiler options
#FFLAGS =  -c -O4
#LFLAGS =  
#LINK    = g++

# overridden library objects
#
#MeshPreProcess   = MeshPreProcess
#MeshPostProcess1 = MeshPostProcess1
#MeshPostProcess2 = MeshPostProcess2

# place to add object modules
#
#ALL_OBJECTS = $(COMMON_OBJECTS)
#ALL_SOURCES = $(COMMON_SOURCES)

### targets ###

.cpp.o: ; $(LINK) $(FFLAGS) $*.cpp -o $*.o

mesh.exe: $(ALL_OBJECTS)
	$(LINK) $(LFLAGS) $(ALL_OBJECTS) -o mesh.exe

include $(MESHGENCPP)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###