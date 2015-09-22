CXX = g++
CXXFLAGS = -Wall -g -std=gnu++11

main: main.o virus.o host.o hostpopulation.o
	$(CXX) $(CXXFLAGS) main.o virus.o host.o hostpopulation.o -o main

main.o: main.cpp virus.hpp host.hpp hostpopulation.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

virus.o: virus.cpp virus.hpp
	$(CXX) $(CXXFLAGS) -c virus.cpp -o virus.o

host.o: host.cpp host.hpp
	$(CXX) $(CXXFLAGS) -c host.cpp -o host.o

hostpopulation.o: hostpopulation.cpp hostpopulation.hpp
	$(CXX) $(CXXFLAGS) -c hostpopulation.cpp -o hostpopulation.o

clean: 
	rm -f *.o
