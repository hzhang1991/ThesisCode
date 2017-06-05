include $(DOGPACK)/lib/sys/Makefile.defs

TEST_SOURCES = \
 $(LIB_SYS)/assert.h \
 $(LIB_SYS)/assert.cpp \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \

TEST_OBJECTS = \
 $(LIB_SYS)/assert.o \
 $(LIB_SYS)/debug.o \

INTERVAL_SOURCES = \
 $(LIB_SYS)/Interval.h \
 $(LIB_SYS)/assert.h \
 $(LIB_SYS)/assert.cpp \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \
 test_Interval.cpp

INTERVAL_OBJECTS = \
 $(LIB_SYS)/assert.o \
 $(LIB_SYS)/debug.o \
 test_Interval.o

ASSERT_SOURCES = \
 $(LIB_SYS)/AssertStream.h \
 $(LIB_SYS)/assert.h \
 $(LIB_SYS)/assert.cpp \
 $(LIB_SYS)/fcmp.h \
 $(LIB_SYS)/fcmp.c \
 test_assert.cpp \

ASSERT_OBJECTS = \
 $(LIB_SYS)/assert.o \
 $(LIB_SYS)/fcmp.o \
 test_assert.o \

INI_SOURCES = \
 $(LIB_SYS)/getdelim.h \
 $(LIB_SYS)/getdelim.c \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \
 $(LIB_SYS)/ErrorStream.h \
 $(LIB_SYS)/IniDocument.h \
 $(LIB_SYS)/IniDocument.cpp \
 test_ini.cpp \

INI_OBJECTS = \
 $(LIB_SYS)/getdelim.o \
 $(LIB_SYS)/debug.o \
 $(LIB_SYS)/IniDocument.o \
 test_ini.o \

DEBUG_SOURCES = \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \
 $(LIB_SYS)/DebugStream.h \
 $(LIB_SYS)/ErrorStream.h \
 test_debug.cpp \

DEBUG_OBJECTS = \
 $(LIB_SYS)/debug.o \
 test_debug.o \

TIME_SOURCES = \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \
 $(LIB_SYS)/ext_time.h \
 $(LIB_SYS)/ext_time.cpp \
 time_test.cpp \

TIME_OBJECTS = \
 $(LIB_SYS)/debug.o \
 $(LIB_SYS)/ext_time.o \
 #time_test.o \

ALL_SOURCES = \
  $(INTERVAL_SOURCES) \
  $(ASSERT_SOURCES) \
  $(INI_SOURCES) \
  $(DEBUG_SOURCES) \

ALL_OBJECTS = \
  $(INTERVAL_OBJECTS) \
  $(ASSERT_OBJECTS) \
  $(INI_OBJECTS) \
  $(DEBUG_OBJECTS) \

CLEAN_TARGETS = \
  clean_assert \
  clean_debug \
  clean_ini \
  clean_Interval \

ALL_TARGETS = \
  test_assert \
  test_debug \
  test_ini \
  test_Interval \
  time_test \

all: $(ALL_TARGETS)

test: test.o $(TEST_OBJECTS) test.cpp $(TEST_SOURCES)
	$(LINK) $(TEST_OBJECTS) -o test

test_assert: $(ASSERT_OBJECTS) $(ASSERT_SOURCES)
	$(LINK) -I.. $(ASSERT_OBJECTS) -o test_assert

test_debug: $(DEBUG_OBJECTS) $(DEBUG_SOURCES)
	$(LINK) $(DEBUG_OBJECTS) -o test_debug

test_ini: $(INI_OBJECTS) $(INI_SOURCES)
	$(COMPILE) -c test_ini.cpp
	$(LINK) $(INI_OBJECTS) -o test_ini

test_Interval:  $(INTERVAL_OBJECTS) $(INTERVAL_SOURCES)
	$(LINK) $(INTERVAL_OBJECTS) -o test_Interval

# To compile with fast math do:
#TIME_FLAGS=-ffast-math make time_test
# or merely
#TIME_FLAGS=-funsafe-math-optimizations make time_test
# which seems to give all the benefit
time_test: time_test.o $(TIME_OBJECTS) $(TIME_SOURCES)
	$(COMPILE) $(TIME_FLAGS) time_test.cpp -o time_test.o
	$(LINK) time_test.o $(TIME_OBJECTS) -o time_test

clean_assert:
	rm -f $(ASSERT_OBJECTS) test_assert

clean_debug:
	rm -f $(DEBUG_OBJECTS) test_debug

clean_ini:
	rm -f $(INI_OBJECTS) test_ini

clean_Interval:
	rm -f $(INTERVAL_OBJECTS) test_Interval

clean_time:
	rm -f $(TIME_OBJECTS) time_test

include $(DOGPACK)/lib/sys/Makefile.targets

