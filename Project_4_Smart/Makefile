# Project_3 Makefile

#First Rule
MP2: hartreefock.o MP2.o molecule.o
	c++ $^ -o $@
molecule.o: molecule.cc
	c++ -c $^
hartreefock.o: hartreefock.cc
	c++ -c $^
MP2.o: MP2.cc
	c++ -c $^
clean:
	/bin/rm -f MP2 *.o core *~
