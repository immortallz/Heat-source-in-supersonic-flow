CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2 -Iinclude

OBJS = cone.o main.o

all: run

project: $(OBJS)
	$(CXX) $(CXXFLAGS) -o project $(OBJS)

run: project
	./project
	rm -f $(OBJS) project

cone.o: cone.cpp r.h
	$(CXX) $(CXXFLAGS) -c cone.cpp

main.o: main.cpp r.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f $(OBJS) project