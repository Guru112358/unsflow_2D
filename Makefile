CXX := g++
ACCFLAGS_1 := -Ofast  -march=native  -fopenmp 


unsflow: main.cpp 
	${CXX} ${ACCFLAGS_1} -o unsflow main.cpp 


clean:
	rm -f *.o  unsflow 
