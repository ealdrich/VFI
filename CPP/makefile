# Use icpc to compile .cpp files
CPP  = icpc

# Add CUDA libraries to the link line
LFLAGS += -lgslcblas

# Include standard optimization flags
CPPFLAGS = -O1 -g -c

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
