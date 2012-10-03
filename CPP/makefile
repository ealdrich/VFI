# Use g++ to compile .cpp files
CPP  = g++
SDIR = .

# Paths to headers and libraries
IGSL = /opt/local/include
LGSL = /opt/local/lib

# Add CUDA libraries to the link line
LFLAGS += -L$(LGSL) -lgslcblas

# Include standard optimization flags
CPPFLAGS = -O2 -g -c -I$(IGSL) -I$(SDIR)

# List of all the objects you need
OBJECTS  = ar1.o kGrid.o vfInit.o binaryVal.o vfStep.o gridMax.o binaryMax.o global.o

# Rule that tells make how to make the program from the objects
main :	main.o $(OBJECTS)
	$(CPP) -o main main.o $(OBJECTS) $(LFLAGS) 

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f main
