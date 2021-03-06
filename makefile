CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O3 -g
LIBS = -lblas -lboost_program_options
TARGET = SPH


default: $(TARGET)

main.o: SPH_main.cpp
	$(CXX) $(CXXFLAGS) -o main.o -c SPH_main.cpp

SPH_class.o: SPH_parallel.cpp SPH_parallel.h
	$(CXX) $(CXXFLAGS) -o SPH_class.o -c SPH_parallel.cpp

DataStructure.o: datastructure/DataStructure.cpp datastructure/DataStructure.h
	$(CXX) $(CXXFLAGS) -o DataStructure.o -c datastructure/DataStructure.cpp

gen_case.o: cases/generate_cases.cpp cases/generate_cases.h
	$(CXX) $(CXXFLAGS) -o gen_case.o -c cases/generate_cases.cpp

parse.o: parse_argument/parse_argument.cpp parse_argument/parse_argument.h
	$(CXX) $(CXXFLAGS) -o parse.o -c parse_argument/parse_argument.cpp

algo.o: algorithm/Findneighbour.cpp algorithm/Findneighbour.h
	$(CXX) $(CXXFLAGS) -o algo.o -c algorithm/Findneighbour.cpp

$(TARGET): main.o SPH_class.o DataStructure.o gen_case.o parse.o algo.o
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.o SPH_class.o DataStructure.o gen_case.o parse.o algo.o $(LIBS)


.PHONY: clean droplet block-drop dam-break 

droplet: $(TARGET)
	mpiexec -np 6 ./$(TARGET) --ic-droplet

block-drop: $(TARGET)
	mpiexec -np 6 ./$(TARGET) --ic-block-drop

dam-break: $(TARGET)
	mpiexec -np 6 ./$(TARGET) --ic-dam-break

clean:
	rm *.o SPH
