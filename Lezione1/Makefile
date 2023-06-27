# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -O3 --std=c++11

# Libraries
LIBS = -lm

# Source files
SRCS11 = Esercizio011.cpp random.cpp
SRCS12 = Esercizio012.cpp random.cpp
SRCS13 = Esercizio013.cpp random.cpp 

# Object files
OBJS11 = $(SRCS11:.cpp=.o)
OBJS12 = $(SRCS12:.cpp=.o)
OBJS13 = $(SRCS13:.cpp=.o)

# Executables
EXEC11 = Esercizio011
EXEC12 = Esercizio012
EXEC13 = Esercizio013

# Default target
all: $(EXEC11) $(EXEC12) $(EXEC13)

# Esercizio011 target
$(EXEC11): $(OBJS11) error.h random.h
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Esercizio012 target
$(EXEC12): $(OBJS12) error.h random.h
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Esercizio013 target
$(EXEC13): $(OBJS13) error.h random.h
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Object files targets
%.o: %.cpp  %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean target
clean:
	rm -rf $(OBJS11) $(OBJS12) $(OBJS13) $(EXEC11) $(EXEC12) $(EXEC13) seed.out *.dat