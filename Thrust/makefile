# Compiler for .cpp files
CPP  = g++

# Use nvcc to compile .cu files
NVCC = nvcc
NVCCFLAGS = -arch sm_20 # For fermi's in keeneland

# Add CUDA Paths
ICUDA = /usr/local/cuda/include
LCUDA = /usr/local/cuda/lib64

# Add CUDA libraries to the link line
LFLAGS += -lcuda -lcudart -L$(LCUDA) -lgomp

# Include standard optimization flags
CPPFLAGS = -O3 -c -I$(ICUDA) -Xcompiler -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP

# List of all the objects you need
OBJECTS  = timer.o ar1.o kGrid.o vfInit.o parameters.o

# Rule that tells make how to make the program from the objects
main :	main.o $(OBJECTS)
	$(CPP) -o main main.o $(OBJECTS) $(LFLAGS) 

# Rule that tells make how to turn a .cu file into a .o
%.o: %.cu
		$(NVCC) ${NVCCFLAGS} $(CPPFLAGS) -c $<

# How does make know how to turn a .cpp into a .o?  It's built-in!
# but if you wanted to type it out it would look like:
# %.o: %.cpp
# 	$(CPP) $(CPPFLAGS) -c $<

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f main
