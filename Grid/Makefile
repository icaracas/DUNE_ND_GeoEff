# Compiler and flags
CXX = g++

# includes and library paths
INCLUDES = -I./Calcs -I$(OSCLIB_INC) -I $(EIGEN_INC)
LIBS = -L$(OSCLIB_LIB) -lOscLib

#compile and link flags
CXXFLAGS = `root-config --cflags --libs` -O2 -Wall $(INCLUDES)
LDFLAGS = `root-config --libs` $(LIBS)

# Target executable
## First target No Osc:
TARGET1 = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_NOOscSpectrumNOCAFAna_FDRateAtND
SRCS1 = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_NOOscSpectrumNOCAFAna_FDRateAtND.C
OBJS1 = $(SRCS1:.C=.o)

#Second target With Osc:
TARGET2 = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_WithOscSpectrumNOCAFAna_FDRateAtND 
# Source files
SRCS2 = NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA_WithOscSpectrumNOCAFAna_FDRateAtND.C \
       Calcs/Calcs.cxx

# Object files
OBJS2 = $(SRCS2:.C=.o)
OBJS2 := $(OBJS2:.cxx=.o)
# Default target builts both
all: $(TARGET1) $(TARGET2) 

# Rule for creating the target
$(TARGET1): $(OBJS1)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(TARGET2): $(OBJS2)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Rule for creating object files
%.o: %.C
	$(CXX) -c $< -o $@ $(CXXFLAGS)
%.o: %.cxx
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Clean up build files
clean:
	rm -f $(OBJS1) $(OBJS2) $(TARGET1) $(TARGET2)
# Run target 1
run1: $(TARGET1)
	./$(TARGET1) myfile1.root
run2: $(TARGET2)
	./$(TARGET2) myfile2.root
