# Compiler and flags
CXX = g++
CXXFLAGS = `root-config --cflags --libs` -O2 -Wall
# Target executable
TARGET = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA
# Source files
SRCS = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA.C
# Object files
OBJS = $(SRCS:.C=.o)
# Default target
all: $(TARGET)
# Rule for creating the target
$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS)
# Rule for creating object files
%.o: %.C
	$(CXX) -c $< -o $@ $(CXXFLAGS)
# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)
# Run the program
run: $(TARGET)
	./$(TARGET) myfile.root
