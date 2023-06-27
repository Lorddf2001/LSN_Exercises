CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio031 : Esercizio031.o random.o
	$(CC) Esercizio031.o random.o -o Esercizio031
Esercizio031.o : Esercizio031.cpp error.h
	$(CC) -c Esercizio031.cpp -o Esercizio031.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio031 seed.out *.dat