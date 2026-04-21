TARGET = genScan
SRC    = $(TARGET).cpp

CC     = g++
CFLAGS = -O2 -Wall -ansi -pedantic -std=c++11
LFLAGS = -o $(TARGET) -fopenmp

$(TARGET): $(SRC)
	$(CC) $(SRC) $(CFLAGS) $(LFLAGS)

run: $(TARGET)
	./$(TARGET) ../tmp/Felis_catus_9.dna.plain CAT

clean:
	rm -f $(TARGET) a.out *.o *~ *#
