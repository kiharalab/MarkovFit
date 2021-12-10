CPP=g++
CFLAGS=-O3 -Wall -Werror
MPICPP=mpicxx
# mpi.h throws warnings so we can't treat them as errors
MPICFLAGS=-O3 -Wall
FINALFLAGS=
LD=$(CPP)

ANN_LIBS=-L../libraries/ann/lib
ANN_INCLUDE=-I../libraries/ann/include
ANN_LDFLAGS=-lANN

PAIRWISE_SCORES_MPI_OBJS=main_pairwise_scores_mpi.o pairwise_scores_mpi.o pairwise_features.o pdb.o utils.o transformable_pdb.o point_transformation.o rmsd.o contact.o scoring.o
PAIRWISE_SCORES_OBJS=main_pairwise_scores.o pairwise_scores.o pairwise_features.o pdb.o utils.o transformable_pdb.o point_transformation.o rmsd.o contact.o scoring.o

PAIRWISE_SCORES_LIBS=$(ANN_LIBS) 
PAIRWISE_SCORES_INCLUDE=$(ANN_INCLUDE) 
PAIRWISE_SCORES_LDFLAGS=$(ANN_LDFLAGS) 

BINARIES= pairwise_scores pairwise_scores_mpi

all: $(BINARIES)

pdb.o:
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) -c $(subst .o,.cc,$@)

pairwise_scores: $(PAIRWISE_SCORES_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(PAIRWISE_SCORES_INCLUDE) $(PAIRWISE_SCORES_LIBS) $(PAIRWISE_SCORES_OBJS) $(PAIRWISE_SCORES_LDFLAGS) -o $@  

pairwise_scores_mpi: $(PAIRWISE_SCORES_MPI_OBJS)
	$(MPICPP) -DWITH_MPI $(CFLAGS) $(FINALFLAGS) $(PAIRWISE_SCORES_INCLUDE) $(PAIRWISE_SCORES_LIBS) $(PAIRWISE_SCORES_MPI_OBJS) $(PAIRWISE_SCORES_LDFLAGS) -o $@  

contact.o scoring.o rmsd.o pairwise_scores.o:
	$(CPP) $(CFLAGS) $(ANN_INCLUDE)  -c $(subst .o,.cc,$@)

pairwise_scores_mpi.o: pairwise_scores.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(CFLAGS) $(ANN_INCLUDE)  -c $(subst _mpi.o,.cc,$@) -o $@

main_pairwise_scores_mpi.o: main_pairwise_scores.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(CFLAGS) $(ANN_INCLUDE)  -c $(subst _mpi.o,.cc,$@) -o $@

install:
	cp $(BINARIES) ~/bin

clean:
	rm -f $(PAIRWISE_SCORES_OBJS) deps $(BINARIES)

.SUFFIXES: .cc .cpp .o

.cc.o:
	$(CPP) $(CFLAGS) -c $< -o $@

.cpp.o:
	$(CPP) $(CFLAGS) -c $< -o $@

deps:
	$(CPP) -MM *.cc  $(ANN_INCLUDE) > deps

include deps
