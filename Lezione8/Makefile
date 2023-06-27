CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio082.exe : Esercizio082.o random.o Hamiltonian_metropolis.h
	$(CC) random.o Esercizio082.o -o Esercizio082.exe
Esercizio082.o : Esercizio082.cpp Hamiltonian_metropolis.h
	$(CC) -c Esercizio082.cpp -o Esercizio082.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o Esercizio082.exe
removedata :
	rm *.dat
