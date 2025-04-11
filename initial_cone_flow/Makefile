prog: cone.o main.o
	g++ main.o cone.o -o prog.exe && del main.o cone.o && prog.exe && del prog.exe
cone.o: cone.cpp
	g++ -c cone.cpp
main.o: main.cpp
	g++ -c main.cpp