# Warnings
WFLAGS	:= -Wall -Wextra -Wsign-conversion -Wsign-compare

# Optimization and architecture
OPT  := -O3
ARCH	:= -march=native

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++11
CXXXSTD := -std=c++14

# Linker options
LDOPT	:= $(OPT)
LDFLAGS := -fopenmp
BIN = "/usr/local/gcc/6.4.0/bin/gcc"
.DEFAULT_GOAL := all

.PHONY: debug
debug : OPT  := -O0 -g -fno-omit-frame-pointer -fsanitize=addres
debug : LDFLAGS := -fsanitize=address
debug : ARCH :=
debug : $(EXEC)

all : gals_advection galsMPIOMP galsMPI

gals_advection: GALS_Advection.cpp
	@ echo Compiling $<...
	$(CXX) GALS_Advection.cpp -o GALS_Advection $(WFLAGS) $(OPT) $(LDFLAGS) $(CFLAGS) $(CXXXSTD) #-pg -fprofile-arcs -ftest-coverage

galsMPIOMP: GALS_Advection_MPI_OpenMP.cpp
	@ echo Compiling $<...
	module load openmpi/2.1.1;mpicxx GALS_Advection_MPI_OpenMP.cpp -o GALS_Advection_MPI_OMP $(WFLAGS) $(OPT) $(LDFLAGS) $(CFLAGS) $(CXXXSTD) #-pg -fprofile-arcs -ftest-coverage

galsMPI: GALS_Advection_MPI.cpp
	@ echo Compiling $<...
	module load openmpi/2.1.1;mpicxx GALS_Advection_MPI.cpp -o GALS_Advection_MPI $(WFLAGS) $(OPT) $(LDFLAGS) $(CFLAGS) $(CXXXSTD) #-pg -fprofile-arcs -ftest-coverage

# TODO: add targets for building executables


.PHONY: clean
clean:
	rm -f *.o *.exe *.dat *.out *.txt
	rm GALS_Advection GALS_Advection_MPI_OMP GALS_Advection_MPI
