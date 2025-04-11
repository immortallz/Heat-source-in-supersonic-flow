prog: Runge-Kutta.o test_RK.o
	g++ Runge-Kutta.o test_RK.o -o prog.exe && del test_RK.o Runge-Kutta.o && prog.exe && del prog.exe
Runge-Kutta.o: Runge-Kutta.cpp
	g++ -c Runge-Kutta.cpp
test_RK.o: test_RK.cpp
	g++ -c test_RK.cpp