OFLAGS = -Wall -g
CFLAGS = -c
CC = g++
exec = test_sym

make: test.o SymmetricAirfoil.o
	$(CC) $(OFLAGS) test.o SymmetricAirfoil.o -o $(exec)

test.o:
	$(CC) $(CFLAGS) test.cpp

SymmetricAirfoil.o:
	$(CC) $(CFLAGS) SymmetricAirfoil.cpp

clean:
	rm -f core *.o
