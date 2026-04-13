CXX = g++-15
CXXFLAGS = -std=c++11 -Wall -O2 -g

OBJS = volterra.o main.o

all: run

volterra: $(OBJS)
	$(CXX) $(CXXFLAGS) -o volterra $(OBJS)

run: volterra
	./volterra

volterra.o: volterra.cpp volterra.h
	$(CXX) $(CXXFLAGS) -c volterra.cpp

main.o: main.cpp volterra.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f $(OBJS) volterra

.PHONY: all run clean
