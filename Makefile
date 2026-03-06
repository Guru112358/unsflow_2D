CXX := g++
ACCFLAGS_1 := -Ofast  -march=native   -pg -g 


unsflow: main.cpp 
	${CXX} ${ACCFLAGS_1} -o unsflow main.cpp 


clean:
	rm -f *.o  unsflow 
