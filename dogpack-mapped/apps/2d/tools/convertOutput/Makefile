
### definitions ###

include $(DOGPACK)/lib/2d/Makefile.defs

DAT2H5_SOURCES = \
  $(DOGPACK)/lib/assert.h \
  $(DOGPACK)/lib/assert.cpp \
  $(DOGPACK)/lib/getdelim.h \
  $(DOGPACK)/lib/getdelim.c \
  $(DOGPACK)/lib/debug.h \
  $(DOGPACK)/lib/debug.cpp \
  $(DOGPACK)/lib/IniDocument.h \
  $(DOGPACK)/lib/IniDocument.cpp \
  $(DOGPACK)/lib/dog_ini.h \
  $(DOGPACK)/lib/dog_ini.cpp \
  $(DOGPACK)/lib/DogParams.h \
  $(DOGPACK)/lib/DogParams.cpp \
  $(DOGPACK)/lib/tensors.h \
  $(DOGPACK)/lib/tensors.cpp \
  $(DOGPACK)/lib/2d/cart/DogParamsCart2.h \
  $(DOGPACK)/lib/2d/cart/DogParamsCart2.cpp \
  $(DOGPACK)/lib/2d/cart/Qinit_restart.cpp \
  $(DOGPACK)/lib/2d/cart/Output.cpp \
  $(DOGPACK)/lib/2d/cart/HDF5.cpp \
  dat2h5.cpp


DAT2H5_OBJECTS = \
  $(DOGPACK)/lib/assert.o \
  $(DOGPACK)/lib/getdelim.o \
  $(DOGPACK)/lib/debug.o \
  $(DOGPACK)/lib/IniDocument.o \
  $(DOGPACK)/lib/dog_ini.o \
  $(DOGPACK)/lib/DogParams.o \
  $(DOGPACK)/lib/tensors.o \
  $(DOGPACK)/lib/2d/cart/DogParamsCart2.o \
  $(DOGPACK)/lib/2d/cart/Qinit_restart.o \
  $(DOGPACK)/lib/2d/cart/Output.o \
  $(DOGPACK)/lib/2d/cart/HDF5.o \

### targets ###

.cpp.o: ; $(COMPILE) $*.cpp -o $*.o

dat2h5: dat2h5.cpp $(DAT2H5_OBJECTS)
	$(LINK) -c dat2h5.cpp \
            -I$(DOGPACK)/lib/2d/cart \
            -I$(DOGPACK)/lib/2d \
            -I$(DOGPACK)/lib \
	    -I$(HDF5_INCLUDES) 
	$(LINK) $(LFLAGS) $(DAT2H5_OBJECTS) dat2h5.o -o dat2h5

cleanallo:
	rm -f $(DAT2H5_OBJECTS) dat2h5

retags:

tags: retags
	ctags $(DAT2H5_SOURCES)
	$(DOGPACK)/scripts/makefiletags $(DAT2H5_SOURCES) >> tags
	sort tags -o tags

