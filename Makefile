CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2 -Iinclude

OBJS = BaseArray.o EFGR.o heatSourceConfiguration.o bodySurface.o normalizedXiConfiguration.o solver.o main.o

all: run

project: $(OBJS)
	$(CXX) $(CXXFLAGS) -o project $(OBJS)

run: project
	./project

BaseArray.o: BaseArray.cpp r.h
	$(CXX) $(CXXFLAGS) -c BaseArray.cpp

EFGR.o: EFGR.cpp r.h
	$(CXX) $(CXXFLAGS) -c EFGR.cpp

heatSourceConfiguration.o: heatSourceConfiguration.cpp r.h
	$(CXX) $(CXXFLAGS) -c heatSourceConfiguration.cpp

bodySurface.o: bodySurface.cpp r.h
	$(CXX) $(CXXFLAGS) -c bodySurface.cpp

normalizedXiConfiguration.o: normalizedXiConfiguration.cpp r.h
	$(CXX) $(CXXFLAGS) -c normalizedXiConfiguration.cpp

solver.o: solver.cpp r.h
	$(CXX) $(CXXFLAGS) -c solver.cpp

main.o: main.cpp r.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f $(OBJS) project