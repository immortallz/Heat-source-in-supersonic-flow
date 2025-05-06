CXX = g++
CXXFLAGS = -std=c++11 -fopenmp -Wall -O2

OBJS = BaseArray.o EFGR.o main.o

all: run

project: $(OBJS)
	$(CXX) $(CXXFLAGS) -o project $(OBJS)

run: project
	./project

BaseArray.o: BaseArray.cpp r.h
	$(CXX) $(CXXFLAGS) -c BaseArray.cpp

EFGR.o: EFGR.cpp r.h
	$(CXX) $(CXXFLAGS) -c EFGR.cpp

main.o: main.cpp r.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f $(OBJS) project