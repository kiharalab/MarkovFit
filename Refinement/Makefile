CPP=g++
CFLAGS=-O3 -Wall -Werror
MPICPP=mpicxx
# mpi.h throws warnings so we can't treat them as errors
MPICFLAGS=-O3 -Wall
FINALFLAGS=
LD=$(CPP)
UNAME=$(shell uname -s)

ANN_LIBS=-L../libraries/ann/lib
ANN_INCLUDE=-I../libraries/ann/include
ANN_LDFLAGS=-lANN


REFINE_MPI_OBJS=main_refine_mpi.o refine_mpi.o mrf_features.o pdb.o utils.o transformable_pdb.o point_transformation.o rmsd.o contact.o scoring.o
REFINE_OBJS=main_refine.o refine.o mrf_features.o pdb.o utils.o transformable_pdb.o point_transformation.o rmsd.o contact.o scoring.o

REFINE_LIBS=$(ANN_LIBS)
REFINE_INCLUDE=$(ANN_INCLUDE)
REFINE_LDFLAGS=$(ANN_LDFLAGS)

ifeq ($(UNAME),Linux)
BINARIES= refine refine_mpi
else
BINARIES= refine refine_mpi
endif

all: $(BINARIES)

pdb.o:
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) -c $(subst .o,.cc,$@)

refine: $(REFINE_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(REFINE_INCLUDE) $(REFINE_LIBS) $(REFINE_OBJS) $(REFINE_LDFLAGS) -o $@  

refine_mpi: $(REFINE_MPI_OBJS)
	$(MPICPP) -DWITH_MPI $(CFLAGS) $(FINALFLAGS) $(REFINE_INCLUDE) $(REFINE_LIBS) $(REFINE_MPI_OBJS) $(REFINE_LDFLAGS) -o $@  


contact.o scoring.o rmsd.o refine.o:
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) -c $(subst .o,.cc,$@)

refine_mpi.o: refine.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(CFLAGS) $(ANN_INCLUDE) -c $(subst _mpi.o,.cc,$@) -o $@

main_refine_mpi.o: main_refine.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(CFLAGS) $(ANN_INCLUDE) -c $(subst _mpi.o,.cc,$@) -o $@

install:
	cp $(BINARIES) ~/bin

clean:
	rm -f $(REFINE_OBJS) deps $(BINARIES)

.SUFFIXES: .cc .cpp .o

.cc.o:
	$(CPP) $(CFLAGS) -c $< -o $@

.cpp.o:
	$(CPP) $(CFLAGS) -c $< -o $@

deps:
	$(CPP) -MM *.cc $(ANN_INCLUDE) > deps

include deps
