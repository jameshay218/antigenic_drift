CXX = g++
CXXFLAGS = -Wall -g

main: main.o virus.o
	$(CXX) $(CXXFLAGS) main.o virus.o -o main

main.o: main.cpp virus.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

Virus.o: virus.cpp virus.hpp
	$(CXX) $(CXXFLAGS) -c virus.cpp -o virus.o

clean: 
	rm -f *.o
