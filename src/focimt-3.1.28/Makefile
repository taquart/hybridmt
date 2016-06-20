CC = g++
CFLAGS = -O3 -Wall 
OBJ = getopts.o faultsolution.o focimtaux.o inputdata.o timedist.o traveltime.o usmtcore.o trinity_library.o

all: focimt

focimt: $(OBJ) moment_tensor.cpp 
	$(CC) $(CFLAGS) moment_tensor.cpp -o focimt $(OBJ) -lcairo

faultsolution.o: faultsolution.cpp 
	$(CC) -c $(CFLAGS) faultsolution.cpp

focimtaux.o: focimtaux.cpp
	$(CC) -c $(CFLAGS) focimtaux.cpp

getopts.o: getopts.cpp
	$(CC) -c $(CFLAGS) getopts.cpp

inputdata.o: inputdata.cpp
	$(CC) -c $(CFLAGS) inputdata.cpp 

timedist.o: timedist.cpp
	$(CC) -c $(CFLAGS) timedist.cpp

usmtcore.o: usmtcore.cpp
	$(CC) -c $(CFLAGS) usmtcore.cpp

traveltime.o: traveltime.cpp
	$(CC) -c $(CFLAGS) traveltime.cpp 

trinity_library.o: trinity_library.cpp
	$(CC) -c $(CFLAGS) trinity_library.cpp
