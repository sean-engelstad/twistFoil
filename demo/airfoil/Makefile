CFLAGS = -O
CC = g++

AIRFOIL_DIR=${TWIST_FOIL}/src/airfoil

make: symmetric_test.o ${AIRFOIL_DIR}/SymmetricAirfoil.o
	$(CC) $(CFLAGS) -o sym_test symmetric_test.o ${AIRFOIL_DIR}/SymmetricAirfoil.o

symmetric_test.o: symmetric_test.cpp ${AIRFOIL_DIR}/SymmetricAirfoil.cpp
	$(CC) $(CFLAGS) -I${AIRFOIL_DIR} -c symmetric_test.cpp

clean:
	rm -f *.o
