# Project_3 Makefile

#First Rule
CCSD: hartreefock.o CCSD.o molecule.o #scf.o
	c++ $^ -g -o $@

#scf.o: scf.cc
#	c++ -g -c $^
molecule.o: molecule.cc
	c++ -g -c $^
hartreefock.o: hartreefock.cc
	c++ -g -c $^
CCSD.o: CCSD.cc
	c++ -g -c $^
clean:
	/bin/rm -f CCSD *.o core *~

#PROG = scf
#CXX = g++
#CXXFLAGS = -g -Wall -std=c++11
#CPPFLAGS = -I$(HOME)/src/eigen -I$(HOME)/src/eigen/unsupported
#CXXLIBS =
#
#$(PROG): $(PROG).o
#    $(CXX) $^ $(CXXLIBS) -o $@
#
#clean:
#    rm -f *.o $(PROG)
