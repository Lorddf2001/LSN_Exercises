CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

Esercizio091.exe : Esercizio091.o random.o map.h route.h population.h
	$(CC) random.o Esercizio091.o -o Esercizio091.exe
Esercizio091.o : Esercizio091.cpp map.h route.h population.h
	$(CC) -c Esercizio091.cpp -o Esercizio091.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

run:
	mpiexec -np 6 Esercizio091.exe

clean :
	rm *.o Esercizio091.exe
removedata :
	rm *.dat
