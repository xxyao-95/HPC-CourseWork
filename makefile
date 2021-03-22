CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O3
LIBS = -lblas
default: SPH_run

main.o: SPH_main.cpp
	$(CXX) $(CXXFLAGS) -o main.o -c SPH_main.cpp

SPH_class.o: SPH_parallel.cpp SPH_parallel.h
	$(CXX) $(CXXFLAGS) -o SPH_class.o -c SPH_parallel.cpp

DataStructure.o: datastructure/DataStructure.cpp datastructure/DataStructure.h
	$(CXX) $(CXXFLAGS) -o DataStructure.o -c datastructure/DataStructure.cpp

gen_case.o: cases/generate_cases.cpp cases/generate_cases.h
	$(CXX) $(CXXFLAGS) -o gen_case.o -c cases/generate_cases.cpp

SPH_run: main.o SPH_class.o DataStructure.o gen_case.o
	$(CXX) $(CXXFLAGS) -o SPH_run main.o SPH_class.o DataStructure.o gen_case.o $(LIBS)

clean:
	rm *.o SPH_run
