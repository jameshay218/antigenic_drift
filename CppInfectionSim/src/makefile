CXX = g++
CXXFLAGS = -Wall -g

main: main.o virus.o host.o hostpopulation.o viruspopulation.o
	$(CXX) $(CXXFLAGS) main.o virus.o host.o hostpopulation.o viruspopulation.o -o main

main.o: main.cpp virus.hpp host.hpp hostpopulation.hpp viruspopulation.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

virus.o: virus.cpp virus.hpp
	$(CXX) $(CXXFLAGS) -c virus.cpp -o virus.o

host.o: host.cpp host.hpp
	$(CXX) $(CXXFLAGS) -c host.cpp -o host.o

hostpopulation.o: hostpopulation.cpp hostpopulation.hpp
	$(CXX) $(CXXFLAGS) -c hostpopulation.cpp -o hostpopulation.o

viruspopulation.o: viruspopulation.cpp viruspopulation.hpp
	$(CXX) $(CXXFLAGS) -c viruspopulation.cpp -o viruspopulation.o

clean: 
	rm -f *.o
