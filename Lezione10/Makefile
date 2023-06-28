CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

Esercizio10.exe : Esercizio10.o random.o map.h route.h population.h
	$(CC) random.o Esercizio10.o -o Esercizio10.exe
Esercizio10.o : Esercizio10.cpp map.h route.h population.h
	$(CC) -c Esercizio10.cpp -o Esercizio10.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

run:
	mpiexec -np 6 Esercizio10.exe

clean :
	rm *.o Esercizio10.exe
removedata :
	rm *.dat
