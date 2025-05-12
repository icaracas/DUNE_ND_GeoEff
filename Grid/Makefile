# Compiler and flags
CXX = g++
# includes and library paths
INCLUDES = -I./Calcs -I$(OSCLIB_INC) -I $(EIGEN_INC)
LIBS = -L$(OSCLIB_LIB) -lOscLib
#compile and link flags
CXXFLAGS = `root-config --cflags --libs` -O2 -Wall $(INCLUDES)
LDFLAGS = `root-config --libs` $(LIBS)
# Target executable
#TARGET = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA
#TARGET = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_NOOscSpectrumNOCAFAna_FDRateAtND
TARGET = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_WithOscSpectrumNOCAFAna_FDRateAtND 
# Source files
#SRCS = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA.C
#SRCS = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_NOOscSpectrumNOCAFAna_FDRateAtND.C
#SRCS = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_WithOscSpectrumNOCAFAna_FDRateAtND.C 
SRCS = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_WithOscSpectrumNOCAFAna_FDRateAtND.C \
       Calcs/Calcs.cxx

# Object files
OBJS = $(SRCS:.C=.o)
OBJS := $(OBJS:.cxx=.o)
# Default target
all: $(TARGET)
# Rule for creating the target
$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)
# Rule for creating object files
%.o: %.C
	$(CXX) -c $< -o $@ $(CXXFLAGS)
%.o: %.cxx
	$(CXX) -c $< -o $@ $(CXXFLAGS)
# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)
# Run the program
run: $(TARGET)
	./$(TARGET) myfile.root
