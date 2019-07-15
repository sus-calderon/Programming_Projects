#This is the title of my executable program
HVF_calculator: molecule.o hessian.o harm_vib_calc.o
	c++ $^ -o $@

molecule.o: molecule.cc
	c++ -c $^

hessian.o: hessian.cc
	c++ -c $^

harm_vib_calc.o: harm_vib_calc.cc
	c++ -c $^

clean:
	/bin/rm -f HVF_calculator *.o core *~

