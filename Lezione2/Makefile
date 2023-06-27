# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -O3 --std=c++11

# Libraries
LIBS = -lm

# Source files
SRCS21 = Esercizio021.cpp random.cpp
SRCS22 = Esercizio022.cpp random.cpp

# Object files
OBJS21 = $(SRCS21:.cpp=.o)
OBJS22 = $(SRCS22:.cpp=.o)

# Executables
EXEC21 = Esercizio021
EXEC22 = Esercizio022

# Default target
all: $(EXEC21) $(EXEC22)

# Esercizio021 target
$(EXEC21): $(OBJS21) error.h function.h integralMC.h distribution.h random.h 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Esercizio022 target error.h randomwalk.h random.h 
$(EXEC22): $(OBJS22)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Object files targets
%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean target
clean:
	rm -rf $(OBJS21) $(OBJS22) $(EXEC21) $(EXEC22) seed.out *.dat