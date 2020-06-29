CC = $(HOME)/hdf5-1.10.1-linux-centos7-x86_64-gcc485-shared/bin/h5cc
CFLAGS = -std=gnu99 -Wall -O3 -DNDEBUG
TARGET = precomp_distance
TARGET_SUB = precomp_distance_module
TARGET_ALL = $(TARGET) $(TARGET_SUB)
TEST = test
CXX = $(HOME)/hdf5-1.10.1-linux-centos7-x86_64-gcc485-shared/bin/h5c++
CXXFLAGS = -std=c++11 -Wall

# For memory leak detection
#CPPUTEST_HOME = $(HOME)/cpputest_home
#CPPFLAGS += -I$(CPPUTEST_HOME)/include
#CXXFLAGS += -include $(CPPUTEST_HOME)/include/CppUTest/MemoryLeakDetectorNewMacros.h
#CFLAGS += -include $(CPPUTEST_HOME)/include/CppUTest/MemoryLeakDetectorMallocMacros.h
#LD_LIBRARIES = -L$(CPPUTEST_HOME)/lib -lCppUTest -lCppUTestExt

all: $(TARGET)

$(TEST): CPPUTEST_HOME = $(HOME)/cpputest_home
$(TEST).o: CPPFLAGS += -I$(CPPUTEST_HOME)/include
$(TEST).o: $(TARGET_SUB).h
$(TEST): LD_LIBRARIES = -L$(CPPUTEST_HOME)/lib -lCppUTest -lCppUTestExt
$(TEST): $(TEST).o $(TARGET_SUB).o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LD_LIBRARIES)

$(TARGET_SUB).o: $(TARGET_SUB).h

$(TARGET): $(TARGET).o $(TARGET_SUB).o

.PHONY: clean
clean:
	$(RM) *.o $(TARGET_ALL) $(TEST) $(TEST).tmp.*
