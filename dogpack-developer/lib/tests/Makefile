include $(DOGPACK)/lib/Makefile.defs

DOGDEFS_SOURCES = \
 $(TOPLIB)/constants.h \
 $(LIB_SYS)/assert.h \
 $(LIB_SYS)/assert.cpp \
 $(LIB_SYS)/fcmp.c \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \
 $(TOPLIB)/tensors.h \
 $(TOPLIB)/tensors.cpp \

DOGDEFS_OBJECTS = \
 $(LIB_SYS)/assert.o \
 $(LIB_SYS)/fcmp.o \
 $(LIB_SYS)/debug.o \
 $(TOPLIB)/tensors.o \

DOGDEFS_DEPENDENCIES = \
 $(DOGDEFS_SOURCES) \
 $(DOGDEFS_OBJECTS) \

TEST_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test.o

TEST_CONSTANTS_OBJECTS = \
  $(CONSTANTS_OBJECTS) \
  test_constants.o

TEST_DOG_MATH_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test_dog_math.o

TEST_TENSOR_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test_tensors.o

TEST_DTENSORBC6_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test_dtensorbc6.o

TIME_SOURCES = \
 $(DOGDEFS_SOURCES) \
 $(LIB_SYS)/ext_time.h \
 $(LIB_SYS)/ext_time.cpp \
 time_test.cpp \

TIME_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  $(LIB_SYS)/ext_time.o \
 #time_test.o \

%.o : %.cpp
	$(LINK) $(FFLAGS) $< -o $@

#.cpp.o: ; $(LINK) $(FFLAGS) $*.cpp -o $*.o

EXECUTABLES = \
  test \
  test_tensors \
  test_dtensorbc6 \
  test_constants \
  test_dog_math \
  time_test \

LFLAGS = -gdb

all: $(EXECUTABLES)

test: $(DOGDEFS_DEPENDENCIES) test.cpp test.o
	$(LINK) $(LFLAGS) $(TEST_OBJECTS) -o test

test_constants: $(DOGDEFS_DEPENDENCIES) test_constants.cpp test_constants.o
	$(LINK) $(LFLAGS) $(TEST_CONSTANTS_OBJECTS) -o test_constants

test_dog_math: $(DOGDEFS_DEPENDENCIES) test_dog_math.cpp test_dog_math.o
	$(LINK) $(LFLAGS) $(TEST_DOG_MATH_OBJECTS) -o test_dog_math

test_tensors:  $(DOGDEFS_DEPENDENCIES) test_tensors.cpp test_tensors.o
	$(LINK) $(LFLAGS) $(TEST_TENSOR_OBJECTS) -o test_tensors

test_dtensorbc6:  $(DOGDEFS_DEPENDENCIES) test_dtensorbc6.cpp test_dtensorbc6.o
	$(LINK) $(LFLAGS) $(TEST_DTENSORBC6_OBJECTS) -o test_dtensorbc6

time_test: time_test.o $(TIME_OBJECTS) $(TIME_SOURCES)
	$(COMPILE) $(TIME_FLAGS) time_test.cpp -o time_test.o
	$(LINK) time_test.o $(TIME_OBJECTS) -o time_test

CLEAN_FILES = *.o $(EXECUTABLES)

include $(DOGPACK)/lib/Makefile.targets
